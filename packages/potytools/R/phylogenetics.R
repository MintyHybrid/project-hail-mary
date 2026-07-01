# Phylogenetic Analysis and Host-Virus Co-phylogeny
# Compare viral gene phylogenies with host taxonomy

# ==============================================================================
# SEQUENCE ALIGNMENT AND DISTANCE CALCULATION
# ==============================================================================

#' Build phylogenetic tree from DNA sequences
#'
#' @param sequences DNAStringSet
#' @param method Distance method ("K80", "JC69", "raw", "TN93")
#' @param tree_method Tree building ("nj", "upgma")
#' @param bootstrap Number of bootstrap replicates (0 = no bootstrap)
#' @return phylo object
build_phylo_tree <- function(sequences,
                             method = "K80",
                             tree_method = "nj",
                             bootstrap = 0) {

  if (inherits(sequences, "DNAStringSet")) {
    seq_names <- names(sequences)
  } else {
    stop("sequences must be DNAStringSet")
  }

  message(sprintf("Building phylogenetic tree from %d sequences...\n", length(sequences)))

  # Convert to phangorn format
  phyDat_seq <- phangorn::as.phyDat(sequences)

  # Calculate distance matrix
  message(sprintf("  Calculating %s distances...\n", method))
  if (method == "K80") {
    dist_matrix <- ape::dist.dna(ape::as.DNAbin(sequences), model = "K80")
  } else if (method == "JC69") {
    dist_matrix <- ape::dist.dna(ape::as.DNAbin(sequences), model = "JC69")
  } else if (method == "TN93") {
    dist_matrix <- ape::dist.dna(ape::as.DNAbin(sequences), model = "TN93")
  } else if (method == "raw") {
    dist_matrix <- ape::dist.dna(ape::as.DNAbin(sequences), model = "raw")
  } else {
    dist_matrix <- phangorn::dist.ml(phyDat_seq)
  }

  # Build tree
  message(sprintf("  Building tree using %s...\n", tree_method))
  if (tree_method == "nj") {
    tree <- ape::njs(dist_matrix)
  } else if (tree_method == "upgma") {
    tree <- phangorn::upgma(dist_matrix)
  } else {
    stop("tree_method must be 'nj' or 'upgma'")
  }

  # Bootstrap support
  if (bootstrap > 0) {
    message(sprintf("  Calculating %d bootstrap replicates...\n", bootstrap))
    bs <- phangorn::bootstrap.phyDat(phyDat_seq,
      FUN = function(x) {
        if (tree_method == "nj") {
          ape::njs(phangorn::dist.ml(x))
        } else {
          phangorn::upgma(phangorn::dist.ml(x))
        }
      },
      bs = bootstrap)

    # Add bootstrap values to tree
    tree$node.label <- ape::prop.clades(tree, bs) / bootstrap * 100
  }

  message("  Done!\n")

  return(tree)
}

#' Build phylogenies from multiple genomic regions
#'
#' @param genome_alignment DNAStringSet of full genomes (aligned)
#' @param regions Named list of (start, end) coordinates
#' @param method Distance method
#' @param bootstrap Number of bootstrap replicates
#' @return Named list of phylo objects
build_multi_region_trees <- function(genome_alignment,
                                     regions,
                                     method = "K80",
                                     bootstrap = 100) {

  trees <- list()

  for (region_name in names(regions)) {
    message(sprintf("\n=== Building tree for %s ===\n", region_name))

    coords <- regions[[region_name]]
    start_pos <- coords[1]
    end_pos <- coords[2]

    # Extract region
    region_seqs <- Biostrings::subseq(genome_alignment, start = start_pos, end = end_pos)

    # Build tree
    tree <- build_phylo_tree(region_seqs,
      method = method,
      tree_method = "nj",
      bootstrap = bootstrap)

    trees[[region_name]] <- tree
  }

  return(trees)
}




# ==============================================================================
# TREE COMPARISON METRICS
# ==============================================================================

#' Compare two phylogenetic trees
#'
#' @param tree1 First phylo object
#' @param tree2 Second phylo object
#' @param labels Common tip labels (if NULL, uses intersection)
#' @return List with comparison metrics
compare_trees <- function(tree1, tree2, labels = NULL) {

  if (is.null(labels)) {
    # Find common labels
    labels <- intersect(tree1$tip.label, tree2$tip.label)
  }

  if (length(labels) < 4) {
    message("Not enough common labels to compare trees\n")
    return(NULL)
  }

  # Prune trees to common labels
  tree1_pruned <- ape::keep.tip(tree1, labels)
  tree2_pruned <- ape::keep.tip(tree2, labels)

  # Ensure trees are binary
  tree1_pruned <- ape::multi2di(tree1_pruned)
  tree2_pruned <- ape::multi2di(tree2_pruned)

  # Calculate Robinson-Foulds distance
  rf_dist <- phangorn::RF.dist(tree1_pruned, tree2_pruned, normalize = TRUE)

  # Path length difference
  path_diff <- phangorn::path.dist(tree1_pruned, tree2_pruned)

  # Cophenetic correlation
  coph_cor <- stats::cor(stats::cophenetic(tree1_pruned), stats::cophenetic(tree2_pruned))

  results <- list(
    n_taxa = length(labels),
    RF_distance = rf_dist,
    path_difference = mean(path_diff),
    cophenetic_correlation = coph_cor
  )

  message("\nTree Comparison Metrics:\n")
  message(sprintf("  Common taxa: %d\n", results$n_taxa))
  message(sprintf("  Robinson-Foulds distance: %.3f\n", results$RF_distance))
  message(sprintf("  Mean path difference: %.3f\n", results$path_difference))
  message(sprintf("  Cophenetic correlation: %.3f\n", results$cophenetic_correlation))
  message("\nInterpretation:\n")
  message("  RF = 0: Identical topologies\n")
  message("  RF = 1: Completely different topologies\n")
  message("  Coph. corr. > 0.8: Similar branch lengths\n\n")

  return(results)
}

# ==============================================================================
# TANGLEGRAM VISUALIZATION
# ==============================================================================

#' Create tanglegram comparing two trees
#'
#' @param tree1 First phylo object
#' @param tree2 Second phylo object
#' @param main_title Plot title
#' @return dendextend tanglegram object
plot_tanglegram <- function(tree1, tree2,
                            main_title = "Tanglegram Comparison") {
  # Convert to dendrograms
  dend1 <- stats::as.dendrogram(ape::as.hclust.phylo(tree1))
  dend2 <- stats::as.dendrogram(ape::as.hclust.phylo(tree2))

  # Create dendlist
  dend_list <- dendextend::dendlist(dend1, dend2)

  # Plot tanglegram
  dendextend::tanglegram(dend1, dend2,
    main = main_title,
    sub = sprintf("Entanglement: %.3f",
      dendextend::entanglement(dend_list)),
    margin_inner = 5,
    lwd = 1,
    edge.lwd = 1,
    type = "rectangle",
    lab.cex = 0.6)
}

# ==============================================================================
# COMPREHENSIVE PHYLOGENETIC ANALYSIS
# ==============================================================================

#' Complete phylogenetic analysis workflow
#'
#' @param genome_alignment DNAStringSet
#' @param regions Named list of genomic regions
#' @param bootstrap Number of bootstrap replicates
#' @return List with all trees and comparisons
phylogenetic_analysis_workflow <- function(genome_alignment,
                                           regions,
                                           bootstrap = 100) {

  message("\n")
  message(paste(rep("=", 80), collapse = ""), "\n")
  message("COMPREHENSIVE PHYLOGENETIC ANALYSIS\n")
  message(paste(rep("=", 80), collapse = ""), "\n\n")

  # Build trees for each region
  message("Step 1: Building phylogenetic trees for each region\n")
  message(paste(rep("-", 60), collapse = ""), "\n")

  viral_trees <- build_multi_region_trees(
    genome_alignment,
    regions,
    method = "K80",
    bootstrap = bootstrap
  )

  # Compare viral trees
  message("\n\nStep 2: Comparing viral phylogenies\n")
  message(paste(rep("-", 60), collapse = ""), "\n")

  tree_comparisons <- list()
  region_names <- names(viral_trees)

  for (i in 1:(length(region_names) - 1)) {
    for (j in (i + 1):length(region_names)) {
      region1 <- region_names[i]
      region2 <- region_names[j]

      message(sprintf("\nComparing %s vs %s:\n", region1, region2))

      comparison <- compare_trees(viral_trees[[region1]],
        viral_trees[[region2]])

      tree_comparisons[[paste(region1, region2, sep = "_vs_")]] <- comparison
    }
  }

  message("\n")
  message(paste(rep("=", 80), collapse = ""), "\n")
  message("PHYLOGENETIC ANALYSIS COMPLETE\n")
  message(paste(rep("=", 80), collapse = ""), "\n\n")

  return(list(
    viral_trees = viral_trees,
    tree_comparisons = tree_comparisons
  ))
}

