# Optimized codon usage analysis functions (seqinr-based)
# Ported from 2026_wsl/R/optimized_functions.R

#' Calculate RSCU using seqinr
#'
#' @param sequences DNAStringSet or character vector
#' @return Data frame with codon, amino_acid, count, frequency, RSCU
#' @export
calculate_rscu_correct <- function(sequences) {
  if (inherits(sequences, "DNAStringSet")) sequences <- as.character(sequences)
  all_seq <- tolower(paste(sequences, collapse = ""))
  seq_vector <- seqinr::s2c(all_seq)
  rscu_values <- seqinr::uco(seq_vector, frame = 0, index = "rscu")
  codon_counts_raw <- seqinr::uco(seq_vector, frame = 0, index = "eff")
  results <- data.frame(
    codon = toupper(names(rscu_values)),
    RSCU = as.numeric(rscu_values),
    stringsAsFactors = FALSE
  )
  results$amino_acid <- sapply(results$codon, function(c) {
    aa <- GENETIC_CODE[[c]]
    if (is.null(aa)) "X" else aa
  })
  results$count <- codon_counts_raw[tolower(results$codon)]
  total_codons <- sum(results$count, na.rm = TRUE)
  results$frequency <- results$count / total_codons
  results <- results[, c("codon", "amino_acid", "count", "frequency", "RSCU")]
  results <- results[results$amino_acid != "*", ]
  results[order(results$amino_acid, results$codon), ]
}

#' Calculate RSCU matrix across isolates (for PCA/clustering)
#'
#' @param sequences DNAStringSet or character vector
#' @return Matrix: sequences x codons with RSCU values
#' @export
calculate_rscu_matrix_fast <- function(sequences) {
  if (inherits(sequences, "DNAStringSet")) {
    seq_names <- names(sequences)
    sequences <- as.character(sequences)
  } else {
    seq_names <- names(sequences)
    if (is.null(seq_names)) seq_names <- paste0("seq_", seq_along(sequences))
  }
  all_codons <- names(GENETIC_CODE)[GENETIC_CODE != "*"]
  rscu_matrix <- matrix(0,
    nrow = length(sequences), ncol = length(all_codons),
    dimnames = list(seq_names, all_codons)
  )
  for (i in seq_along(sequences)) {
    seq_vector <- seqinr::s2c(tolower(sequences[i]))
    tryCatch(
      {
        rscu_vals <- seqinr::uco(seq_vector, frame = 0, index = "rscu")
        codon_upper <- toupper(names(rscu_vals))
        idx <- codon_upper %in% all_codons
        rscu_matrix[i, codon_upper[idx]] <- rscu_vals[idx]
      },
      error = function(e) {
        warning(sprintf("Could not calculate RSCU for sequence %d: %s", i, e$message))
      }
    )
  }
  rscu_matrix[is.na(rscu_matrix)] <- 0
  rscu_matrix
}

#' PCA on per-isolate RSCU values
#'
#' @param sequences DNAStringSet or character vector
#' @param host_classification Optional data frame with columns isolate, host_type
#' @return List with pca, data, variance, rscu_matrix
#' @export
perform_rscu_pca_improved <- function(sequences, host_classification = NULL) {
  rscu_matrix_unfiltered <- calculate_rscu_matrix_fast(sequences)
  vars <- apply(rscu_matrix_unfiltered, 2, var)
  rscu_matrix <- rscu_matrix_unfiltered[, vars > 0, drop = FALSE]
  pca_result <- stats::prcomp(rscu_matrix, scale. = TRUE)
  pca_data <- data.frame(
    isolate = rownames(rscu_matrix),
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    PC3 = pca_result$x[, 3],
    stringsAsFactors = FALSE
  )
  if (!is.null(host_classification)) {
    pca_data <- merge(pca_data, host_classification[, c("isolate", "host_type")],
      by = "isolate", all.x = TRUE
    )
  }
  var_explained <- summary(pca_result)$importance[2, ] * 100
  list(
    pca = pca_result, data = pca_data, variance = var_explained,
    rscu_matrix = rscu_matrix, rscu_matrix_unfiltered = rscu_matrix_unfiltered
  )
}

#' Fast GC content by codon position
#'
#' @param sequence Character string or DNAString
#' @return Named vector: GC_overall, GC1, GC2, GC3
#' @export
calculate_gc_fast <- function(sequence) {
  if (inherits(sequence, "DNAString")) sequence <- as.character(sequence)
  sv <- seqinr::s2c(tolower(sequence))
  len <- length(sv)
  gc_pos <- function(pos) if (length(pos) > 0) seqinr::GC(pos) else 0
  c(
    GC_overall = seqinr::GC(sv),
    GC1 = gc_pos(sv[seq(1, len, by = 3)]),
    GC2 = gc_pos(sv[seq(2, len, by = 3)]),
    GC3 = gc_pos(sv[seq(3, len, by = 3)])
  )
}

#' GC content matrix for multiple sequences
#'
#' @param sequences DNAStringSet or character vector
#' @return Matrix: sequences x c(GC_overall, GC1, GC2, GC3)
#' @export
calculate_gc_matrix_fast <- function(sequences) {
  if (!inherits(sequences, "DNAStringSet")) {
    sequences <- Biostrings::DNAStringSet(as.character(sequences))
  }
  seq_names <- names(sequences)
  n <- length(sequences)

  # Overall GC in one vectorized C-level call.
  af <- Biostrings::alphabetFrequency(sequences, baseOnly = TRUE)
  acgt <- rowSums(af[, c("A", "C", "G", "T"), drop = FALSE])
  gc_overall <- ifelse(acgt > 0, (af[, "G"] + af[, "C"]) / acgt, 0)

  # Per-codon-position GC (GC1/GC2/GC3).
  w <- Biostrings::width(sequences)
  pos_gc <- if (n > 0 && length(unique(w)) == 1L && w[1] > 0L) {
    # Equal-length input (e.g. an alignment): fully vectorized matrix path.
    m <- as.matrix(sequences) # n x L character matrix, built in C
    is_gc <- m == "G" | m == "C"
    is_acgt <- is_gc | m == "A" | m == "T"
    len <- ncol(m)
    matrix(vapply(1:3, function(p) {
      cols <- seq(p, len, by = 3)
      num <- rowSums(is_gc[, cols, drop = FALSE])
      den <- rowSums(is_acgt[, cols, drop = FALSE])
      ifelse(den > 0, num / den, 0)
    }, numeric(n)), nrow = n)
  } else if (n > 0) {
    # Variable-length fallback: lightweight base-R per sequence (no seqinr).
    t(vapply(as.character(sequences), function(s) {
      sv <- strsplit(toupper(s), "", fixed = TRUE)[[1]]
      is_gc <- sv == "G" | sv == "C"
      is_acgt <- is_gc | sv == "A" | sv == "T"
      vapply(1:3, function(p) {
        idx <- seq(p, length(sv), by = 3)
        den <- sum(is_acgt[idx])
        if (den > 0) sum(is_gc[idx]) / den else 0
      }, numeric(1))
    }, numeric(3), USE.NAMES = FALSE))
  } else {
    matrix(numeric(0), nrow = 0, ncol = 3)
  }

  gc_matrix <- cbind(GC_overall = gc_overall, GC1 = pos_gc[, 1],
                     GC2 = pos_gc[, 2], GC3 = pos_gc[, 3])
  rownames(gc_matrix) <- seq_names
  gc_matrix
}

#' Fast ENC calculation via seqinr codon counts
#'
#' @param sequence Character string or DNAString
#' @return Numeric ENC value (NA if sequence too short)
#' @export
calculate_enc_fast <- function(sequence) {
  if (inherits(sequence, "DNAString")) sequence <- as.character(sequence)
  sv <- seqinr::s2c(tolower(sequence))
  if (length(sv) < 15) {
    return(NA_real_)
  }
  tryCatch(
    {
      codon_counts <- seqinr::uco(sv, frame = 0, index = "eff")
      aa_codon_counts <- list()
      for (codon in names(codon_counts)) {
        aa <- GENETIC_CODE[[toupper(codon)]]
        if (!is.null(aa) && aa != "*") {
          aa_codon_counts[[aa]][[toupper(codon)]] <- codon_counts[codon]
        }
      }
      families <- list("2" = c(), "3" = c(), "4" = c(), "6" = c())
      for (aa in names(aa_codon_counts)) {
        n_syn <- length(AA_TO_CODONS[[aa]])
        if (n_syn %in% c(2, 3, 4, 6)) {
          cts <- unlist(aa_codon_counts[[aa]])
          tot <- sum(cts)
          if (tot > 1) {
            families[[as.character(n_syn)]] <- c(
              families[[as.character(n_syn)]],
              sum((cts / tot)^2)
            )
          }
        }
      }
      # Degenerate input (very short, or dominated by single-codon amino
      # acids such as Met/Trp): no synonymous family has usable data, so
      # Wright's Nc is undefined. Return NA rather than the misleading base
      # value of 2 (the [20, 61] range assumes a realistic amino-acid spread).
      if (!any(lengths(families) > 0)) {
        return(NA_real_)
      }
      enc <- 2
      for (fs in c("2", "3", "4", "6")) {
        F <- if (length(families[[fs]]) > 0) mean(families[[fs]]) else 0
        if (F > 0) {
          n_aa <- sum(sapply(
            names(aa_codon_counts),
            function(aa) length(AA_TO_CODONS[[aa]]) == as.numeric(fs)
          ))
          enc <- enc + n_aa / F
        }
      }
      enc
    },
    error = function(e) NA_real_
  )
}

#' Fast codon usage comparison (motif vs reference)
#'
#' @param motif_sequences DNAStringSet or character vector
#' @param reference_sequences DNAStringSet or character vector
#' @return Data frame with per-codon comparison statistics
#' @export
compare_codon_usage_fast <- function(motif_sequences, reference_sequences) {
  motif_rscu <- calculate_rscu_correct(motif_sequences)
  ref_rscu <- calculate_rscu_correct(reference_sequences)
  comparison <- merge(
    motif_rscu[, c("codon", "amino_acid", "count", "frequency", "RSCU")],
    ref_rscu[, c("codon", "count", "frequency", "RSCU")],
    by = "codon", suffixes = c("_motif", "_ref"), all = TRUE
  )
  comparison[is.na(comparison)] <- 0
  comparison$fold_change <- ifelse(
    comparison$frequency_ref > 0,
    comparison$frequency_motif / comparison$frequency_ref,
    ifelse(comparison$frequency_motif > 0, Inf, 1.0)
  )
  comparison$p_value <- NA
  motif_total <- sum(comparison$count_motif)
  ref_total <- sum(comparison$count_ref)
  for (i in seq_len(nrow(comparison))) {
    if ((comparison$count_motif[i] + comparison$count_ref[i]) >= 5) {
      ct <- matrix(
        c(
          comparison$count_motif[i], motif_total - comparison$count_motif[i],
          comparison$count_ref[i], ref_total - comparison$count_ref[i]
        ),
        nrow = 2, byrow = TRUE
      )
      tryCatch(
        {
          comparison$p_value[i] <- stats::fisher.test(ct)$p.value
        },
        error = function(e) {}
      )
    }
  }
  comparison$p_adjusted <- stats::p.adjust(comparison$p_value, method = "BH")
  comparison[order(-comparison$fold_change), ]
}

#' Complete optimized codon usage analysis
#'
#' @param motif_sequences DNAStringSet - the motif/query sequences
#' @param reference_sequences DNAStringSet - background reference sequences
#' @param host_codon_table Optional host codon table (unused, kept for compatibility)
#' @param output_prefix Prefix for CSV output files
#' @return Named list: rscu_comparison, usage_comparison, motif_rscu, ref_rscu,
#'   enc_motif, enc_ref, gc_motif, gc_ref
#' @export
run_complete_analysis_optimized <- function(motif_sequences,
                                            reference_sequences,
                                            host_codon_table = NULL,
                                            output_prefix = "analysis") {
  motif_char <- as.character(motif_sequences)
  ref_char <- as.character(reference_sequences)

  message(sprintf(
    "Dataset: %d isolates; motif length: %d nt",
    length(motif_char), nchar(motif_char[1])
  ))

  motif_rscu <- calculate_rscu_correct(motif_char)
  ref_rscu <- calculate_rscu_correct(ref_char)

  rscu_comparison <- merge(
    motif_rscu[, c("codon", "amino_acid", "RSCU")],
    ref_rscu[, c("codon", "RSCU")],
    by = "codon", suffixes = c("_motif", "_ref")
  )
  rscu_comparison$RSCU_diff <- rscu_comparison$RSCU_motif - rscu_comparison$RSCU_ref
  utils::write.csv(rscu_comparison, paste0(output_prefix, "_rscu_comparison.csv"),
    row.names = FALSE
  )

  usage_comparison <- compare_codon_usage_fast(motif_char, ref_char)
  utils::write.csv(usage_comparison, paste0(output_prefix, "_usage_comparison.csv"),
    row.names = FALSE
  )

  enc_motif <- stats::na.omit(sapply(motif_char, calculate_enc_fast))
  enc_ref <- stats::na.omit(sapply(ref_char, calculate_enc_fast))

  gc_motif <- calculate_gc_matrix_fast(motif_char)
  gc_ref <- calculate_gc_matrix_fast(ref_char)

  list(
    rscu_comparison  = rscu_comparison,
    usage_comparison = usage_comparison,
    motif_rscu       = motif_rscu,
    ref_rscu         = ref_rscu,
    enc_motif        = enc_motif,
    enc_ref          = enc_ref,
    gc_motif         = gc_motif,
    gc_ref           = gc_ref
  )
}
