# Extract Motifs with Variable Flanking Sequences
# For structural analysis and comparative studies

# ==============================================================================
# MOTIF + FLANKS EXTRACTION
# ==============================================================================

#' Extract motif with flanking sequences
#'
#' @param sequences DNAStringSet of genomes or aligned sequences
#' @param motif_start Start position of motif in first sequence
#' @param motif_length Length of core motif
#' @param flank_5prime Length of 5' flank
#' @param flank_3prime Length of 3' flank
#' @param extract_from "alignment" or "unaligned"
#' @return List with motif, 5' flank, 3' flank, and combined sequences
extract_motif_with_flanks <- function(sequences,
                                      motif_start,
                                      motif_length,
                                      flank_5prime = 0,
                                      flank_3prime = 0,
                                      extract_from = "alignment") {
  if (extract_from == "alignment") {
    # Extract from aligned sequences (positions are alignment coordinates)

    # Calculate positions
    total_start <- motif_start - flank_5prime
    total_end <- motif_start + motif_length - 1 + flank_3prime

    # Ensure valid range
    if (total_start < 1) {
      warning("5' flank extends before sequence start, truncating")
      total_start <- 1
    }

    if (total_end > width(sequences)[1]) {
      warning("3' flank extends beyond sequence end, truncating")
      total_end <- width(sequences)[1]
    }

    # Extract combined sequence
    combined <- subseq(sequences, start = total_start, end = total_end)

    # Extract individual parts
    flank_5 <- if (flank_5prime > 0) {
      subseq(sequences,
        start = total_start,
        end = motif_start - 1
      )
    } else {
      DNAStringSet(rep("", length(sequences)))
    }

    motif <- subseq(sequences,
      start = motif_start,
      end = motif_start + motif_length - 1
    )

    flank_3 <- if (flank_3prime > 0) {
      subseq(sequences,
        start = motif_start + motif_length,
        end = total_end
      )
    } else {
      DNAStringSet(rep("", length(sequences)))
    }
  } else {
    # Extract from unaligned sequences
    # Assumes motif position is known for each sequence

    combined <- DNAStringSet()
    flank_5 <- DNAStringSet()
    motif <- DNAStringSet()
    flank_3 <- DNAStringSet()

    for (i in seq_along(sequences)) {
      seq <- sequences[i]

      # For unaligned, motif_start should be a vector
      if (length(motif_start) > 1) {
        start_pos <- motif_start[i]
      } else {
        start_pos <- motif_start
      }

      total_start <- max(1, start_pos - flank_5prime)
      total_end <- min(width(seq), start_pos + motif_length - 1 + flank_3prime)

      combined <- c(combined, subseq(seq, start = total_start, end = total_end))

      if (flank_5prime > 0) {
        flank_5 <- c(flank_5, subseq(seq,
          start = total_start,
          end = start_pos - 1
        ))
      }

      motif <- c(motif, subseq(seq,
        start = start_pos,
        end = start_pos + motif_length - 1
      ))

      if (flank_3prime > 0) {
        flank_3 <- c(flank_3, subseq(seq,
          start = start_pos + motif_length,
          end = total_end
        ))
      }
    }

    names(combined) <- names(sequences)
    names(flank_5) <- names(sequences)
    names(motif) <- names(sequences)
    names(flank_3) <- names(sequences)
  }

  cat(sprintf("Extracted motif with flanks:\n"))
  cat(sprintf("  Core motif: %d nt\n", motif_length))
  cat(sprintf("  5' flank: %d nt\n", flank_5prime))
  cat(sprintf("  3' flank: %d nt\n", flank_3prime))
  cat(sprintf(
    "  Total length: %d nt\n",
    motif_length + flank_5prime + flank_3prime
  ))
  cat(sprintf("  Number of sequences: %d\n", length(sequences)))

  return(list(
    combined = combined,
    flank_5 = flank_5,
    motif = motif,
    flank_3 = flank_3,
    config = list(
      motif_start = motif_start,
      motif_length = motif_length,
      flank_5prime = flank_5prime,
      flank_3prime = flank_3prime
    )
  ))
}

#' Generate multiple motif+flank combinations
#'
#' @param sequences DNAStringSet
#' @param motif_start Motif start position
#' @param motif_length Motif length
#' @param flank_sizes Vector of flank sizes to test
#' @return Named list of extraction results
generate_flank_variants <- function(sequences,
                                    motif_start,
                                    motif_length,
                                    flank_sizes = c(0, 10, 20, 30, 50)) {
  variants <- list()

  for (size in flank_sizes) {
    variant_name <- sprintf("flank_%d", size)

    cat(sprintf("\nGenerating %s variant...\n", variant_name))

    variants[[variant_name]] <- extract_motif_with_flanks(
      sequences,
      motif_start = motif_start,
      motif_length = motif_length,
      flank_5prime = size,
      flank_3prime = size
    )
  }

  cat(sprintf("\nGenerated %d flank variants\n", length(variants)))

  return(variants)
}

#' Save motif+flank sequences to FASTA files
#'
#' @param extraction_result Result from extract_motif_with_flanks
#' @param output_prefix Prefix for output files
save_motif_flanks <- function(extraction_result, output_prefix = "motif") {
  # Save combined
  writeXStringSet(extraction_result$combined,
    filepath = paste0(output_prefix, "_combined.fasta")
  )

  # Save motif only
  writeXStringSet(extraction_result$motif,
    filepath = paste0(output_prefix, "_core.fasta")
  )

  # Save flanks if they exist
  if (extraction_result$config$flank_5prime > 0) {
    writeXStringSet(extraction_result$flank_5,
      filepath = paste0(output_prefix, "_5prime_flank.fasta")
    )
  }

  if (extraction_result$config$flank_3prime > 0) {
    writeXStringSet(extraction_result$flank_3,
      filepath = paste0(output_prefix, "_3prime_flank.fasta")
    )
  }

  cat(sprintf("Saved sequences with prefix: %s\n", output_prefix))
}

# ==============================================================================
# HOMOLOG DETECTION AND EXTRACTION
# ==============================================================================

#' Find motif homologs in sequences using BLAST-like approach
#'
#' @param query_motif DNAString of query motif
#' @param target_sequences DNAStringSet to search
#' @param max_mismatches Maximum allowed mismatches
#' @return Data frame with positions of matches
find_motif_homologs <- function(query_motif,
                                target_sequences,
                                max_mismatches = 2) {
  if (inherits(query_motif, "DNAStringSet")) {
    query_motif <- query_motif[[1]]
  }
  if (is.character(query_motif)) {
    query_motif <- Biostrings::DNAString(query_motif)
  }
  if (is.character(target_sequences)) {
    target_sequences <- Biostrings::DNAStringSet(target_sequences)
  }

  query_char <- as.character(query_motif)
  query_length <- length(query_motif)

  # Biostrings does the fixed-mismatch sliding-window search in C - far faster
  # than an R-level per-position substr()/strsplit() loop.
  hits <- Biostrings::vmatchPattern(query_motif, target_sequences,
    max.mismatch = max_mismatches
  )

  seq_names <- names(target_sequences)
  if (is.null(seq_names)) seq_names <- as.character(seq_along(target_sequences))

  # Assemble one data frame per target, then bind once (no rbind-in-loop).
  per_seq <- lapply(seq_along(target_sequences), function(i) {
    mi <- hits[[i]]
    if (length(mi) == 0) {
      return(NULL)
    }
    starts <- BiocGenerics::start(mi)
    windows <- as.character(Biostrings::extractAt(target_sequences[[i]], mi))
    # Exact mismatch count per window (vmatchPattern already bounded it)
    q_chars <- strsplit(query_char, "", fixed = TRUE)[[1]]
    mismatches <- vapply(windows, function(w) {
      sum(q_chars != strsplit(w, "", fixed = TRUE)[[1]])
    }, integer(1))
    data.frame(
      sequence = seq_names[i],
      start = starts,
      end = starts + query_length - 1L,
      mismatches = mismatches,
      identity = 1 - (mismatches / query_length),
      match_seq = windows,
      stringsAsFactors = FALSE
    )
  })

  matches <- do.call(rbind, per_seq)
  if (is.null(matches)) matches <- data.frame()

  if (nrow(matches) > 0) {
    message(sprintf("Found %d homologous regions", nrow(matches)))
  } else {
    message("No homologs found")
  }

  matches
}

#' Extract homologous motifs with flanks
#'
#' @param sequences DNAStringSet
#' @param homolog_positions Data frame from find_motif_homologs
#' @param flank_5prime 5' flank length
#' @param flank_3prime 3' flank length
#' @return DNAStringSet of extracted sequences
extract_homologs_with_flanks <- function(sequences,
                                         homolog_positions,
                                         flank_5prime = 0,
                                         flank_3prime = 0) {
  extracted <- DNAStringSet()

  for (i in 1:nrow(homolog_positions)) {
    hit <- homolog_positions[i, ]

    seq_idx <- which(names(sequences) == hit$sequence)
    if (length(seq_idx) == 0) next

    seq <- sequences[seq_idx]

    total_start <- max(1, hit$start - flank_5prime)
    total_end <- min(width(seq), hit$end + flank_3prime)

    extracted_seq <- subseq(seq, start = total_start, end = total_end)
    names(extracted_seq) <- sprintf(
      "%s_pos%d_id%.2f",
      hit$sequence, hit$start, hit$identity
    )

    extracted <- c(extracted, extracted_seq)
  }

  cat(sprintf("Extracted %d homologous sequences\n", length(extracted)))

  return(extracted)
}

# ==============================================================================
# COMPLETE WORKFLOW
# ==============================================================================

#' Complete motif extraction workflow
#'
#' @param sequences DNAStringSet
#' @param motif_start Start position
#' @param motif_length Length of motif
#' @param flank_variants Vector of flank sizes to generate
#' @param find_homologs Whether to find homologs
#' @param output_dir Output directory
#' @return List of results
workflow_motif_extraction <- function(sequences,
                                      motif_start,
                                      motif_length,
                                      flank_variants = c(0, 10, 20, 30),
                                      find_homologs = TRUE,
                                      output_dir = "motif_extraction") {
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("MOTIF + FLANKS EXTRACTION WORKFLOW\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Generate variants
  cat("Generating flank variants...\n")
  variants <- generate_flank_variants(
    sequences,
    motif_start,
    motif_length,
    flank_sizes = flank_variants
  )

  # Save each variant
  for (variant_name in names(variants)) {
    cat(sprintf("\nSaving %s...\n", variant_name))
    save_motif_flanks(
      variants[[variant_name]],
      output_prefix = file.path(output_dir, variant_name)
    )
  }

  # Find homologs if requested
  homolog_results <- NULL
  if (find_homologs) {
    cat("\nSearching for motif homologs...\n")

    query_motif <- variants[[1]]$motif[[1]]

    homologs <- find_motif_homologs(
      query_motif,
      sequences,
      max_mismatches = 2
    )

    if (nrow(homologs) > 0) {
      # Extract homologs with flanks
      for (variant_name in names(variants)) {
        flank_size <- variants[[variant_name]]$config$flank_5prime

        homolog_seqs <- extract_homologs_with_flanks(
          sequences,
          homologs,
          flank_5prime = flank_size,
          flank_3prime = flank_size
        )

        writeXStringSet(
          homolog_seqs,
          filepath = file.path(
            output_dir,
            paste0("homologs_", variant_name, ".fasta")
          )
        )
      }

      homolog_results <- homologs
    }
  }

  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("COMPLETE\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  cat(sprintf("Output saved to: %s/\n", output_dir))

  return(list(
    variants = variants,
    homologs = homolog_results
  ))
}
