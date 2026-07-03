# Sliding-window composition statistics over a nucleotide alignment
# Ported from 2026_cleanPoty/enhanced_dashboardwphylo_extended.qmd

# Internal metric helpers ------------------------------------------------------

.gc_content <- function(x) mean(x %in% c("G", "C"), na.rm = TRUE)
.at_content <- function(x) mean(x %in% c("A", "T", "U"), na.rm = TRUE)

.shannon_entropy <- function(x) {
  x <- x[x %in% c("A", "C", "G", "T", "U")]
  if (length(x) == 0) {
    return(NA_real_)
  }
  p <- table(x) / length(x)
  -sum(p * log2(p))
}

.cpg_oe <- function(x) {
  x <- x[x %in% c("A", "C", "G", "T", "U")]
  if (length(x) < 2) {
    return(NA_real_)
  }
  obs <- sum(x[-length(x)] == "C" & x[-1] == "G")
  expct <- sum(x == "C") * sum(x == "G") / length(x)
  if (expct == 0) {
    return(NA_real_)
  }
  obs / expct
}

#' Sliding-window composition statistics across an alignment
#'
#' Slides a window along the columns of a nucleotide alignment and, pooling all
#' sequences within each window, computes GC content, AT content, GC at the
#' third codon position (GC3), the CpG observed/expected ratio, Shannon entropy,
#' gap fraction, and per-codon-position GC. Useful for mapping regional
#' compositional bias (e.g. CpG depletion, GC3 wobble adaptation) along a genome
#' or gene alignment.
#'
#' @param aln Character matrix (rows = sequences, columns = sites), e.g.
#'   `as.matrix(Biostrings::readDNAStringSet(...))`.
#' @param window Window width in columns (default 300).
#' @param step Step between window starts, in columns (default 30).
#' @param codon_aligned If `TRUE` (default), `window` and `step` must be
#'   multiples of 3 so codon positions line up.
#' @param remove_gaps If `TRUE`, drop gap characters before computing metrics.
#' @return Data frame, one row per window, with columns `window`, `start`,
#'   `end`, `GC`, `AT`, `GC3`, `CpG_OE`, `gap_fraction`, `entropy`,
#'   `GC_pos1`, `GC_pos2`, `GC_pos3`.
#' @export
rolling_alignment_stats <- function(aln, window = 300, step = 30,
                                    codon_aligned = TRUE, remove_gaps = FALSE) {
  if (!is.matrix(aln)) {
    stop("Alignment must be a character matrix (rows = seqs, cols = sites)")
  }
  if (codon_aligned && (window %% 3 != 0 || step %% 3 != 0)) {
    stop("For codon-aligned analysis, window and step must be multiples of 3")
  }

  n_sites <- ncol(aln)
  results <- list()
  win_id <- 1L

  for (start in seq(1, n_sites - window + 1, by = step)) {
    end <- start + window - 1
    vec <- as.vector(aln[, start:end, drop = FALSE])
    gap_prop <- mean(vec == "-")
    if (remove_gaps) vec <- vec[vec != "-"]

    codon_pos <- rep(1:3, length.out = window)
    pos1 <- vec[codon_pos == 1]
    pos2 <- vec[codon_pos == 2]
    pos3 <- vec[codon_pos == 3]

    results[[win_id]] <- data.frame(
      window = win_id, start = start, end = end,
      GC = .gc_content(vec), AT = .at_content(vec), GC3 = .gc_content(pos3),
      CpG_OE = .cpg_oe(vec), gap_fraction = gap_prop,
      entropy = .shannon_entropy(vec),
      GC_pos1 = .gc_content(pos1), GC_pos2 = .gc_content(pos2),
      GC_pos3 = .gc_content(pos3)
    )
    win_id <- win_id + 1L
  }

  do.call(rbind, results)
}

#' Per-sequence sliding-window composition statistics
#'
#' As [rolling_alignment_stats()] but computes the metrics separately for each
#' sequence (row) of the alignment, so isolate-to-isolate variation across the
#' alignment can be plotted as a multi-line trace.
#'
#' @inheritParams rolling_alignment_stats
#' @return Data frame, one row per sequence x window, with a leading
#'   `sequence` column plus `window`, `start`, `end`, `GC`, `GC3`, `CpG_OE`,
#'   `entropy`, `gap_fraction`, `GC_pos1`, `GC_pos2`, `GC_pos3`.
#' @export
rolling_alignment_stats_per_sequence <- function(aln, window = 300, step = 30,
                                                 codon_aligned = TRUE,
                                                 remove_gaps = FALSE) {
  if (!is.matrix(aln)) stop("Alignment must be a character matrix")
  if (codon_aligned && (window %% 3 != 0 || step %% 3 != 0)) {
    stop("window and step must be multiples of 3 for codon-aligned CDS")
  }

  n_sites <- ncol(aln)
  results <- vector("list", length = nrow(aln))

  for (i in seq_len(nrow(aln))) {
    seq_name <- rownames(aln)[i]
    seq_i <- aln[i, ]
    win_id <- 1L
    seq_res <- list()

    for (start in seq(1, n_sites - window + 1, by = step)) {
      end <- start + window - 1
      window_seq <- seq_i[start:end]
      gap_fraction <- mean(window_seq == "-")
      if (remove_gaps) window_seq <- window_seq[window_seq != "-"]

      codon_pos <- rep(1:3, length.out = window)
      pos1 <- window_seq[codon_pos == 1]
      pos2 <- window_seq[codon_pos == 2]
      pos3 <- window_seq[codon_pos == 3]

      seq_res[[win_id]] <- data.frame(
        sequence = seq_name, window = win_id, start = start, end = end,
        GC = .gc_content(window_seq), GC3 = .gc_content(pos3),
        CpG_OE = .cpg_oe(window_seq), entropy = .shannon_entropy(window_seq),
        gap_fraction = gap_fraction,
        GC_pos1 = .gc_content(pos1), GC_pos2 = .gc_content(pos2),
        GC_pos3 = .gc_content(pos3),
        stringsAsFactors = FALSE
      )
      win_id <- win_id + 1L
    }
    results[[i]] <- do.call(rbind, seq_res)
  }

  do.call(rbind, results)
}
