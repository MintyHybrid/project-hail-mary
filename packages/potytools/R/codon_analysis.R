# Viral Motif Codon Usage Analysis for R
# Works with DNAStringSet and AAStringSet from Biostrings

# ==============================================================================
# GENETIC CODE DEFINITIONS
# ==============================================================================

# Standard genetic code table
GENETIC_CODE <- list(
  TTT = "F", TTC = "F", TTA = "L", TTG = "L",
  TCT = "S", TCC = "S", TCA = "S", TCG = "S",
  TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
  TGT = "C", TGC = "C", TGA = "*", TGG = "W",
  CTT = "L", CTC = "L", CTA = "L", CTG = "L",
  CCT = "P", CCC = "P", CCA = "P", CCG = "P",
  CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
  CGT = "R", CGC = "R", CGA = "R", CGG = "R",
  ATT = "I", ATC = "I", ATA = "I", ATG = "M",
  ACT = "T", ACC = "T", ACA = "T", ACG = "T",
  AAT = "N", AAC = "N", AAA = "K", AAG = "K",
  AGT = "S", AGC = "S", AGA = "R", AGG = "R",
  GTT = "V", GTC = "V", GTA = "V", GTG = "V",
  GCT = "A", GCC = "A", GCA = "A", GCG = "A",
  GAT = "D", GAC = "D", GAA = "E", GAG = "E",
  GGT = "G", GGC = "G", GGA = "G", GGG = "G"
)

# Amino acid to codons mapping
get_aa_to_codons <- function() {
  aa_to_codons <- list()
  for (codon in names(GENETIC_CODE)) {
    aa <- GENETIC_CODE[[codon]]
    if (is.null(aa_to_codons[[aa]])) {
      aa_to_codons[[aa]] <- c()
    }
    aa_to_codons[[aa]] <- c(aa_to_codons[[aa]], codon)
  }
  return(aa_to_codons)
}

AA_TO_CODONS <- get_aa_to_codons()

# ==============================================================================
# CORE CODON EXTRACTION FUNCTIONS
# ==============================================================================

#' Extract codons from DNA sequence
#'
#' @param seq Character string or DNAString
#' @param frame Starting frame (0, 1, or 2)
#' @return Character vector of codons
#' @export
extract_codons <- function(seq, frame = 0) {
  if (inherits(seq, "DNAString") || inherits(seq, "RNAString")) {
    seq <- as.character(seq)
  }
  seq <- toupper(gsub("U", "T", seq))
  starts <- seq(frame + 1, nchar(seq) - 2, by = 3)
  codons <- substring(seq, starts, starts + 2)
  codons[nchar(codons) == 3 & !grepl("N", codons)]
}

#' Extract codons from DNAStringSet
#'
#' @param dna_set DNAStringSet object
#' @return List of character vectors of codons
#' @export
extract_codons_from_set <- function(dna_set) {
  lapply(as.character(dna_set), extract_codons)
}

#' Get all codons from multiple sequences
#'
#' @param sequences Character vector or DNAStringSet
#' @return Character vector of all codons
get_all_codons <- function(sequences) {
  if (inherits(sequences, "DNAStringSet")) {
    sequences <- as.character(sequences)
  }
  unlist(lapply(sequences, extract_codons), use.names = FALSE)
}

# ==============================================================================
# CODON FREQUENCY CALCULATIONS
# ==============================================================================

#' Calculate codon frequencies
#'
#' @param sequences Character vector or DNAStringSet
#' @return Named vector of codon frequencies
#' @export
calculate_codon_frequencies <- function(sequences) {
  all_codons <- get_all_codons(sequences)
  codon_counts <- table(all_codons)
  codon_freq <- codon_counts / sum(codon_counts)
  return(codon_freq)
}

#' Calculate RSCU (Relative Synonymous Codon Usage)
#'
#' @param sequences Character vector or DNAStringSet
#' @return Data frame with codon, amino acid, and RSCU values
#' @export
calculate_rscu <- function(sequences) {
  if (inherits(sequences, "DNAStringSet")) {
    sequences <- as.character(sequences)
  }
  ### manually checked whether width(ref_seqs) %% 3 all == 0!!! 20260423
  # Concatenate all sequences
  all_seq <- paste(sequences, collapse = "")
  all_seq <- tolower(all_seq) # seqinr expects lowercase

  # Convert to seqinr format
  seq_vector <- s2c(all_seq)

  # Calculate codon usage with seqinr
  # This gives counts for each codon
  codon_counts <- uco(seq_vector, frame = 0, index = "rscu")

  # Get the RSCU values (already calculated by seqinr)
  rscu_values <- codon_counts

  return(rscu_values)
}

# ==============================================================================
# EFFECTIVE NUMBER OF CODONS (ENC)
# ==============================================================================

#' Calculate Effective Number of Codons (ENC/Nc)
#' Wright (1990) Gene 87:23-29
#'
#' @param sequence Single sequence (character or DNAString)
#' @return Numeric ENC value (2-61), or NA if sequence too short
#' @export
calculate_enc <- function(sequence) {
  # Delegates to seqinr-based implementation in optimized_analysis.R
  calculate_enc_fast(sequence)
}

# ==============================================================================
# GC CONTENT ANALYSIS
# ==============================================================================

#' Calculate GC content overall and by codon position
#'
#' @param sequence Character string or DNAString
#' @return Named vector with GC_overall, GC1, GC2, GC3
#' @export
calculate_gc_content <- function(sequence) {
  # Delegates to seqinr-based implementation in optimized_analysis.R
  calculate_gc_fast(sequence)
}

# ==============================================================================
# STATISTICAL COMPARISON
# ==============================================================================

#' Compare codon usage between two sets of sequences
#'
#' @param motif_sequences Character vector or DNAStringSet (motif sequences)
#' @param reference_sequences Character vector or DNAStringSet (reference)
#' @return Data frame with statistical comparison
#' @export
compare_codon_usage <- function(motif_sequences, reference_sequences) {
  # Delegates to vectorized seqinr-based implementation in optimized_analysis.R
  compare_codon_usage_fast(motif_sequences, reference_sequences)
}

# ==============================================================================
# WOBBLE POSITION ANALYSIS
# ==============================================================================

#' Check for wobble (3rd) position bias
#'
#' @param sequences Character vector or DNAStringSet
#' @return Data frame with chi-square test results
#' @export
check_wobble_bias <- function(sequences) {
  wobble_bases <- c()

  for (seq in as.character(sequences)) {
    seq <- toupper(gsub("U", "T", seq))
    for (i in seq(3, nchar(seq), by = 3)) {
      wobble_bases <- c(wobble_bases, substr(seq, i, i))
    }
  }

  counts <- table(wobble_bases)
  total <- sum(counts)
  expected <- total / 4

  results <- data.frame(
    base = c("A", "T", "G", "C"),
    observed = sapply(
      c("A", "T", "G", "C"),
      function(b) ifelse(b %in% names(counts), counts[[b]], 0)
    ),
    expected = expected,
    stringsAsFactors = FALSE
  )

  results$obs_freq <- results$observed / total
  results$deviation <- (results$observed - results$expected) / results$expected

  # Chi-square test
  if (total >= 20) {
    chi_result <- chisq.test(results$observed)
    message("\nWobble position bias test:\n")
    message(sprintf("  chi^2 = %.4f, p = %.4e\n", chi_result$statistic, chi_result$p.value))
  }

  return(results)
}

# ==============================================================================
# DINUCLEOTIDE FREQUENCY ANALYSIS
# ==============================================================================

#' Calculate dinucleotide frequencies (rho values)
#'
#' @param sequences Character vector or DNAStringSet
#' @return Data frame with observed/expected ratios
#' @export
calculate_dinucleotide_freq <- function(sequences) {
  all_seq <- paste(as.character(sequences), collapse = "")
  all_seq <- gsub("N", "", all_seq)
  all_seq <- toupper(gsub("U", "T", all_seq))

  # Mononucleotide frequencies
  mono_counts <- table(strsplit(all_seq, "")[[1]])
  total_mono <- sum(mono_counts)
  mono_freq <- mono_counts / total_mono

  # Dinucleotide counts
  dinucs <- substring(all_seq, 1:(nchar(all_seq) - 1), 2:nchar(all_seq))
  di_counts <- table(dinucs)
  total_di <- sum(di_counts)

  results <- data.frame(
    dinucleotide = character(),
    observed = numeric(),
    obs_freq = numeric(),
    exp_freq = numeric(),
    rho = numeric(),
    stringsAsFactors = FALSE
  )

  for (di in names(di_counts)) {
    base1 <- substr(di, 1, 1)
    base2 <- substr(di, 2, 2)

    obs_freq <- di_counts[[di]] / total_di
    exp_freq <- ifelse(base1 %in% names(mono_freq) && base2 %in% names(mono_freq),
      mono_freq[[base1]] * mono_freq[[base2]], 0
    )

    rho <- ifelse(exp_freq > 0, obs_freq / exp_freq, NA)

    results <- rbind(
      results,
      data.frame(
        dinucleotide = di,
        observed = di_counts[[di]],
        obs_freq = obs_freq,
        exp_freq = exp_freq,
        rho = rho,
        stringsAsFactors = FALSE
      )
    )
  }

  results <- results[order(results$rho, decreasing = TRUE), ]

  return(results)
}

# ==============================================================================
# MOTIF-SPECIFIC VARIANT ANALYSIS
# ==============================================================================

#' Analyze variation in 18nt motifs across isolates
#'
#' @param motif_sequences Character vector or DNAStringSet (must be 18nt)
#' @return Data frame with codon-by-codon variation
#' @export
analyze_motif_variants <- function(motif_sequences) {
  sequences <- as.character(motif_sequences)

  # Filter to 18nt sequences
  # sequences <- sequences[nchar(sequences) == 18]
  #
  # if (length(sequences) == 0) {
  #   stop("No valid 18nt sequences provided")
  # }
  #
  results <- data.frame(
    codon_position = integer(),
    nt_position = character(),
    n_codon_variants = integer(),
    n_aa_variants = integer(),
    codon_entropy = numeric(),
    aa_entropy = numeric(),
    synonymous_var_only = logical(),
    most_common_codon = character(),
    most_common_codon_freq = numeric(),
    most_common_aa = character(),
    stringsAsFactors = FALSE
  )

  for (codon_pos in 1:6) {
    start <- (codon_pos - 1) * 3 + 1
    end <- start + 2

    codons_at_pos <- substring(sequences, start, end)
    codons_at_pos <- toupper(codons_at_pos)

    amino_acids <- sapply(codons_at_pos, function(c) {
      ifelse(c %in% names(GENETIC_CODE), GENETIC_CODE[[c]], "X")
    })

    # Count variants
    codon_counts <- table(codons_at_pos)
    aa_counts <- table(amino_acids)

    # Calculate entropy
    codon_freqs <- as.numeric(codon_counts / sum(codon_counts))
    codon_entropy <- -sum(codon_freqs * log2(codon_freqs + 1e-10))

    aa_freqs <- as.numeric(aa_counts / sum(aa_counts))
    aa_entropy <- -sum(aa_freqs * log2(aa_freqs + 1e-10))

    n_unique_codons <- length(codon_counts)
    n_unique_aa <- length(aa_counts)

    synonymous_var_only <- n_unique_codons > n_unique_aa && n_unique_aa == 1

    most_common <- names(codon_counts)[which.max(codon_counts)]
    most_common_freq <- max(codon_counts) / length(codons_at_pos)
    most_common_aa <- names(aa_counts)[which.max(aa_counts)]

    results <- rbind(
      results,
      data.frame(
        codon_position = codon_pos,
        nt_position = sprintf("%d-%d", start, end),
        n_codon_variants = n_unique_codons,
        n_aa_variants = n_unique_aa,
        codon_entropy = codon_entropy,
        aa_entropy = aa_entropy,
        synonymous_var_only = synonymous_var_only,
        most_common_codon = most_common,
        most_common_codon_freq = most_common_freq,
        most_common_aa = most_common_aa,
        stringsAsFactors = FALSE
      )
    )
  }

  return(results)
}

# ==============================================================================
# RARE CODON DETECTION
# ==============================================================================

#' Identify rare codons in motif compared to reference
#'
#' @param motif_sequences Character vector or DNAStringSet
#' @param reference_freq Named vector of reference codon frequencies
#' @param threshold Frequency threshold for "rare" (default 0.1)
#' @return Data frame with rare codon analysis
#' @export
check_rare_codons <- function(motif_sequences, reference_freq, threshold = 0.1) {
  motif_codons <- get_all_codons(motif_sequences)
  motif_counts <- table(motif_codons)

  results <- data.frame(
    codon = character(),
    amino_acid = character(),
    motif_count = integer(),
    ref_frequency = numeric(),
    is_rare = logical(),
    is_suboptimal = logical(),
    optimal_codon = character(),
    optimal_freq = numeric(),
    stringsAsFactors = FALSE
  )

  for (codon in names(motif_counts)) {
    aa <- GENETIC_CODE[[codon]]
    ref_freq <- ifelse(codon %in% names(reference_freq),
      reference_freq[[codon]], 0
    )

    # Get synonymous codons and their frequencies
    synonymous <- AA_TO_CODONS[[aa]]
    syn_freqs <- sapply(synonymous, function(c) {
      ifelse(c %in% names(reference_freq), reference_freq[[c]], 0)
    })

    max_syn_freq <- max(syn_freqs)
    optimal_codon <- synonymous[which.max(syn_freqs)]

    is_rare <- ref_freq < threshold
    is_suboptimal <- ref_freq < max_syn_freq * 0.5

    results <- rbind(
      results,
      data.frame(
        codon = codon,
        amino_acid = aa,
        motif_count = motif_counts[[codon]],
        ref_frequency = ref_freq,
        is_rare = is_rare,
        is_suboptimal = is_suboptimal,
        optimal_codon = optimal_codon,
        optimal_freq = max_syn_freq,
        stringsAsFactors = FALSE
      )
    )
  }

  results <- results[order(results$ref_frequency), ]

  return(results)
}

# ==============================================================================
# MAIN ANALYSIS PIPELINE
# ==============================================================================

#' Complete analysis pipeline
#'
#' @param motif_sequences DNAStringSet or character vector (18nt motifs)
#' @param reference_sequences DNAStringSet or character vector (full ORFs)
#' @param host_codon_table Optional named vector of host codon frequencies
#' @param output_prefix Prefix for output files
#' @return List of analysis results
#' @export
run_complete_analysis <- function(motif_sequences,
                                  reference_sequences,
                                  host_codon_table = NULL,
                                  output_prefix = "viral_motif_analysis") {
  message("=" %s+% paste(rep("=", 80), collapse = "") %s+% "\n")
  message("VIRAL MOTIF CODON USAGE ANALYSIS\n")
  message("=" %s+% paste(rep("=", 80), collapse = "") %s+% "\n\n")

  # Convert to character if needed
  if (inherits(motif_sequences, "DNAStringSet")) {
    message("Converting DNAStringSet to character vectors...\n")
    motif_seqs <- as.character(motif_sequences)
  } else {
    motif_seqs <- motif_sequences
  }

  if (inherits(reference_sequences, "DNAStringSet")) {
    ref_seqs <- as.character(reference_sequences)
  } else {
    ref_seqs <- reference_sequences
  }

  message(sprintf("Dataset: %d isolates\n", length(motif_seqs)))
  message(sprintf(
    "Motif length: %d nt (%d codons)\n",
    nchar(motif_seqs[1]), nchar(motif_seqs[1]) / 3
  ))

  # 1. Motif variant analysis
  message("\n" %s+% paste(rep("=", 80), collapse = "") %s+% "\n")
  message("MOTIF VARIANT ANALYSIS\n")
  message(paste(rep("=", 80), collapse = "") %s+% "\n\n")

  motif_variants <- analyze_motif_variants(motif_seqs)
  print(motif_variants)
  write.csv(motif_variants,
    paste0(output_prefix, "_motif_variants.csv"),
    row.names = FALSE
  )

  # 2. RSCU analysis
  message("\n" %s+% paste(rep("=", 80), collapse = "") %s+% "\n")
  message("RELATIVE SYNONYMOUS CODON USAGE (RSCU)\n")
  message(paste(rep("=", 80), collapse = "") %s+% "\n\n")

  motif_rscu <- calculate_rscu(motif_seqs)
  ref_rscu <- calculate_rscu(ref_seqs)

  rscu_comparison <- merge(motif_rscu, ref_rscu,
    by = c("codon", "amino_acid"),
    suffixes = c("_motif", "_ref")
  )
  rscu_comparison$RSCU_diff <- rscu_comparison$RSCU_motif - rscu_comparison$RSCU_ref
  rscu_comparison <- rscu_comparison[order(abs(rscu_comparison$RSCU_diff),
    decreasing = TRUE
  ), ]

  message("Top 20 codons with largest RSCU differences:\n")
  print(head(rscu_comparison, 20))
  write.csv(rscu_comparison,
    paste0(output_prefix, "_rscu_comparison.csv"),
    row.names = FALSE
  )

  # 3. Codon usage comparison
  message("\n" %s+% paste(rep("=", 80), collapse = "") %s+% "\n")
  message("CODON USAGE COMPARISON (MOTIF vs REFERENCE)\n")
  message(paste(rep("=", 80), collapse = "") %s+% "\n\n")

  usage_comparison <- compare_codon_usage(motif_seqs, ref_seqs)
  significant <- usage_comparison[!is.na(usage_comparison$p_adjusted) &
    usage_comparison$p_adjusted < 0.05, ]

  if (nrow(significant) > 0) {
    message(sprintf(
      "%d codons show significant usage differences (FDR < 0.05):\n",
      nrow(significant)
    ))
    print(significant[, c("codon", "amino_acid", "fold_change", "p_adjusted")])
  } else {
    message("No significant codon usage differences detected\n")
  }

  write.csv(usage_comparison,
    paste0(output_prefix, "_codon_usage_comparison.csv"),
    row.names = FALSE
  )

  # 4. ENC calculation
  message("\n" %s+% paste(rep("=", 80), collapse = "") %s+% "\n")
  message("EFFECTIVE NUMBER OF CODONS (ENC)\n")
  message(paste(rep("=", 80), collapse = "") %s+% "\n\n")

  enc_motif <- sapply(motif_seqs, calculate_enc)
  enc_motif <- enc_motif[!is.na(enc_motif)]

  enc_ref <- sapply(ref_seqs, calculate_enc)
  enc_ref <- enc_ref[!is.na(enc_ref)]

  if (length(enc_motif) > 0) {
    message(sprintf(
      "Motif ENC: %.2f +/- %.2f\n",
      mean(enc_motif), sd(enc_motif)
    ))
    message(sprintf("  Range: %.2f - %.2f\n", min(enc_motif), max(enc_motif)))
  }

  if (length(enc_ref) > 0) {
    message(sprintf(
      "\nReference ENC: %.2f +/- %.2f\n",
      mean(enc_ref), sd(enc_ref)
    ))
    message(sprintf("  Range: %.2f - %.2f\n", min(enc_ref), max(enc_ref)))
  }

  if (length(enc_motif) > 0 && length(enc_ref) > 0) {
    t_result <- t.test(enc_motif, enc_ref)
    message(sprintf("\nT-test p-value: %.4e\n", t_result$p.value))
    if (t_result$p.value < 0.05) {
      message("  -> Significant difference in codon usage bias\n")
    }
  }

  message("\nInterpretation: ENC=20 (extreme bias) to ENC=61 (no bias)\n")

  # 5. GC content
  message("\n" %s+% paste(rep("=", 80), collapse = "") %s+% "\n")
  message("GC CONTENT ANALYSIS\n")
  message(paste(rep("=", 80), collapse = "") %s+% "\n\n")

  gc_motif <- t(sapply(motif_seqs, calculate_gc_content))
  gc_ref <- t(sapply(ref_seqs, calculate_gc_content))

  message("Motif GC content:\n")
  print(colMeans(gc_motif))

  message("\nReference GC content:\n")
  print(colMeans(gc_ref))

  # 6. Wobble position
  message("\n" %s+% paste(rep("=", 80), collapse = "") %s+% "\n")
  message("WOBBLE (3rd) POSITION BIAS\n")
  message(paste(rep("=", 80), collapse = "") %s+% "\n\n")

  message("Motif wobble position:\n")
  wobble_motif <- check_wobble_bias(motif_seqs)
  print(wobble_motif)

  message("\nReference wobble position:\n")
  wobble_ref <- check_wobble_bias(ref_seqs)
  print(wobble_ref)

  # 7. Dinucleotide analysis
  message("\n" %s+% paste(rep("=", 80), collapse = "") %s+% "\n")
  message("DINUCLEOTIDE FREQUENCY ANALYSIS\n")
  message(paste(rep("=", 80), collapse = "") %s+% "\n\n")

  dinuc_motif <- calculate_dinucleotide_freq(motif_seqs)
  message("Top 10 over-represented dinucleotides in motif (rho > 1):\n")
  print(head(dinuc_motif, 10))

  message("\nTop 10 under-represented dinucleotides in motif (rho < 1):\n")
  print(tail(dinuc_motif, 10))

  write.csv(dinuc_motif,
    paste0(output_prefix, "_dinucleotide_analysis.csv"),
    row.names = FALSE
  )

  # 8. Rare codon check
  message("\n" %s+% paste(rep("=", 80), collapse = "") %s+% "\n")
  message("RARE CODON ANALYSIS\n")
  message(paste(rep("=", 80), collapse = "") %s+% "\n\n")

  ref_freq <- calculate_codon_frequencies(ref_seqs)
  rare_codons <- check_rare_codons(motif_seqs, ref_freq, threshold = 0.1)

  rare_found <- rare_codons[rare_codons$is_rare, ]
  if (nrow(rare_found) > 0) {
    message("Rare codons detected in motif:\n")
    print(rare_found)
  } else {
    message("No rare codons detected\n")
  }

  suboptimal <- rare_codons[rare_codons$is_suboptimal, ]
  if (nrow(suboptimal) > 0) {
    message("\nSuboptimal codons:\n")
    print(suboptimal[, c(
      "codon", "amino_acid", "ref_frequency",
      "optimal_codon", "optimal_freq"
    )])
  }

  write.csv(rare_codons,
    paste0(output_prefix, "_rare_codons.csv"),
    row.names = FALSE
  )

  # 9. Host comparison (if provided)
  if (!is.null(host_codon_table)) {
    message("\n" %s+% paste(rep("=", 80), collapse = "") %s+% "\n")
    message("HOST CODON USAGE COMPARISON\n")
    message(paste(rep("=", 80), collapse = "") %s+% "\n\n")

    host_rare <- check_rare_codons(motif_seqs, host_codon_table, threshold = 0.1)
    host_rare_found <- host_rare[host_rare$is_rare, ]

    if (nrow(host_rare_found) > 0) {
      message("Codons rare in host:\n")
      print(host_rare_found)
    } else {
      message("No codons rare in host detected\n")
    }

    write.csv(host_rare,
      paste0(output_prefix, "_host_comparison.csv"),
      row.names = FALSE
    )
  }

  message("\n" %s+% paste(rep("=", 80), collapse = "") %s+% "\n")
  message("ANALYSIS COMPLETE\n")
  message(sprintf("Results saved with prefix: %s\n", output_prefix))
  message(paste(rep("=", 80), collapse = "") %s+% "\n")

  # Return results as list
  return(list(
    motif_variants = motif_variants,
    rscu_comparison = rscu_comparison,
    usage_comparison = usage_comparison,
    enc_motif = enc_motif,
    enc_ref = enc_ref,
    gc_motif = gc_motif,
    gc_ref = gc_ref,
    wobble_motif = wobble_motif,
    wobble_ref = wobble_ref,
    dinuc_motif = dinuc_motif,
    rare_codons = rare_codons
  ))
}
