# Sliding-window composition statistics over a nucleotide alignment.
#
# Vectorised cumulative-sum implementation: each metric is derived from prefix
# sums of per-site indicator vectors, so a whole sequence is processed in O(L)
# regardless of window/step. This replaces an earlier per-window R loop that
# built one data.frame per window (O(n_seq * n_windows) allocations) and was
# orders of magnitude slower on genome-length alignments.

# Prefix sum with a leading zero: cs(x)[k + 1] == sum(x[1:k]); window sum over
# columns [a, b] is cs[b + 1] - cs[a].
.prefix <- function(x) c(0, cumsum(as.numeric(x)))
.winsum <- function(cs, a, b) cs[b + 1L] - cs[a]

# Codon-position masks (1-based, global frame) for a length-n sequence.
.codon_mask <- function(n, k) (((seq_len(n) - 1L) %% 3L) == (k - 1L))

# Core: given per-site indicator prefix sums for one sequence (or pooled across
# sequences) plus window bounds, assemble the per-window statistics data.frame.
.roll_from_prefix <- function(pref, starts, ends, window, remove_gaps) {
  ws <- function(cs) .winsum(cs, starts, ends)

  n_gc <- ws(pref$gc)
  n_at <- ws(pref$at)
  n_c <- ws(pref$c)
  n_g <- ws(pref$g)
  n_gap <- ws(pref$gap)
  n_valid <- ws(pref$valid)
  # CpG dinucleotides fully inside a window occupy query positions start..end-1.
  n_cpg <- pref$cpg[ends] - pref$cpg[starts]

  # Base denominators: full window, or non-gap count when remove_gaps.
  base_n <- if (remove_gaps) pmax(n_valid, 1) else window
  gc <- n_gc / base_n
  at <- n_at / base_n

  # Per-codon-position GC (third position = GC3).
  gcp <- lapply(1:3, function(k) {
    num <- ws(pref$gc_pos[[k]])
    den <- if (remove_gaps) pmax(ws(pref$n_pos[[k]]), 1) else (window / 3)
    num / den
  })

  # Shannon entropy over A/C/G/T within the window (valid bases only).
  n_a <- ws(pref$a)
  n_t <- ws(pref$t)
  ent <- vapply(seq_along(starts), function(i) {
    counts <- c(n_a[i], n_c[i], n_g[i], n_t[i])
    tot <- sum(counts)
    if (tot == 0) {
      return(NA_real_)
    }
    p <- counts[counts > 0] / tot
    -sum(p * log2(p))
  }, numeric(1))

  # CpG observed/expected: obs / (C * G / n_valid).
  expct <- (n_c * n_g) / pmax(n_valid, 1)
  cpg_oe <- ifelse(expct > 0, n_cpg / expct, NA_real_)

  list(
    gc = gc, at = at, gc3 = gcp[[3]], cpg_oe = cpg_oe,
    gap_fraction = n_gap / window, entropy = ent,
    gc_pos1 = gcp[[1]], gc_pos2 = gcp[[2]], gc_pos3 = gcp[[3]]
  )
}

# Prefix sums of the site indicators for a single character vector (one row).
.seq_prefix <- function(v) {
  v[v == "U"] <- "T"
  n <- length(v)
  is_g <- v == "G"
  is_c <- v == "C"
  is_a <- v == "A"
  is_t <- v == "T"
  is_gc <- is_g | is_c
  valid <- is_a | is_c | is_g | is_t
  masks <- lapply(1:3, .codon_mask, n = n)
  cpg <- if (n >= 2) c(is_c[-n] & is_g[-1], FALSE) else rep(FALSE, n)
  list(
    gc = .prefix(is_gc), at = .prefix(is_a | is_t),
    a = .prefix(is_a), c = .prefix(is_c), g = .prefix(is_g), t = .prefix(is_t),
    gap = .prefix(v == "-"), valid = .prefix(valid),
    cpg = .prefix(cpg),
    gc_pos = lapply(masks, function(m) .prefix(is_gc & m)),
    n_pos = lapply(masks, function(m) .prefix(valid & m))
  )
}

.window_bounds <- function(n_sites, window, step) {
  starts <- seq(1L, n_sites - window + 1L, by = step)
  list(starts = starts, ends = starts + window - 1L)
}

.check_args <- function(aln, window, step, codon_aligned) {
  if (!is.matrix(aln)) {
    stop("Alignment must be a character matrix (rows = seqs, cols = sites)")
  }
  if (codon_aligned && (window %% 3 != 0 || step %% 3 != 0)) {
    stop("For codon-aligned analysis, window and step must be multiples of 3")
  }
}

#' Sliding-window composition statistics across an alignment
#'
#' Slides a window along the columns of a nucleotide alignment and, pooling all
#' sequences within each window, computes GC content, AT content, GC at the
#' third codon position (GC3), the CpG observed/expected ratio, Shannon entropy,
#' gap fraction, and per-codon-position GC. Useful for mapping regional
#' compositional bias (e.g. CpG depletion, GC3 wobble adaptation) along a genome
#' or gene alignment. Implemented with vectorised prefix sums (fast on
#' genome-length alignments).
#'
#' @param aln Character matrix (rows = sequences, columns = sites), e.g.
#'   `as.matrix(Biostrings::readDNAStringSet(...))`.
#' @param window Window width in columns (default 300).
#' @param step Step between window starts, in columns (default 30).
#' @param codon_aligned If `TRUE` (default), `window` and `step` must be
#'   multiples of 3 so codon positions line up.
#' @param remove_gaps If `TRUE`, gap characters are excluded from the metric
#'   denominators (composition computed over aligned bases only).
#' @return Data frame, one row per window, with columns `window`, `start`,
#'   `end`, `GC`, `AT`, `GC3`, `CpG_OE`, `gap_fraction`, `entropy`,
#'   `GC_pos1`, `GC_pos2`, `GC_pos3`.
#' @export
rolling_alignment_stats <- function(aln, window = 300, step = 30,
                                    codon_aligned = TRUE, remove_gaps = FALSE) {
  .check_args(aln, window, step, codon_aligned)
  n_sites <- ncol(aln)
  b <- .window_bounds(n_sites, window, step)

  # Pooled per-column counts across all sequences, then prefix sums.
  mat <- aln
  mat[mat == "U"] <- "T"
  col_g <- colSums(mat == "G")
  col_c <- colSums(mat == "C")
  col_a <- colSums(mat == "A")
  col_t <- colSums(mat == "T")
  col_gc <- col_g + col_c
  col_valid <- col_a + col_c + col_g + col_t
  col_gap <- colSums(mat == "-")
  masks <- lapply(1:3, .codon_mask, n = n_sites)
  # Pooled CpG per query position i: C at column i and G at column i+1, any seq.
  col_cpg <- numeric(n_sites)
  if (n_sites >= 2) {
    is_cg <- (mat[, -n_sites, drop = FALSE] == "C") &
      (mat[, -1, drop = FALSE] == "G")
    col_cpg[-n_sites] <- colSums(is_cg)
  }
  pref <- list(
    gc = .prefix(col_gc), at = .prefix(col_a + col_t),
    a = .prefix(col_a), c = .prefix(col_c), g = .prefix(col_g),
    t = .prefix(col_t),
    gap = .prefix(col_gap), valid = .prefix(col_valid), cpg = .prefix(col_cpg),
    gc_pos = lapply(masks, function(m) .prefix(col_gc * m)),
    n_pos = lapply(masks, function(m) .prefix(col_valid * m))
  )

  # For the pooled case the "window length" is n_seq * window columns.
  n_seq <- nrow(aln)
  s <- .roll_from_prefix(pref, b$starts, b$ends, window * n_seq, remove_gaps)

  data.frame(
    window = seq_along(b$starts), start = b$starts, end = b$ends,
    GC = s$gc, AT = s$at, GC3 = s$gc3, CpG_OE = s$cpg_oe,
    gap_fraction = s$gap_fraction, entropy = s$entropy,
    GC_pos1 = s$gc_pos1, GC_pos2 = s$gc_pos2, GC_pos3 = s$gc_pos3
  )
}

#' Per-sequence sliding-window composition statistics
#'
#' As [rolling_alignment_stats()] but computes the metrics separately for each
#' sequence (row) of the alignment, so isolate-to-isolate variation across the
#' alignment can be plotted as a multi-line trace. Also vectorised per sequence.
#'
#' @inheritParams rolling_alignment_stats
#' @return Data frame, one row per sequence x window, with a leading
#'   `sequence` column plus `window`, `start`, `end`, `GC`, `GC3`, `CpG_OE`,
#'   `entropy`, `gap_fraction`, `GC_pos1`, `GC_pos2`, `GC_pos3`.
#' @export
rolling_alignment_stats_per_sequence <- function(aln, window = 300, step = 30,
                                                 codon_aligned = TRUE,
                                                 remove_gaps = FALSE) {
  .check_args(aln, window, step, codon_aligned)
  n_sites <- ncol(aln)
  b <- .window_bounds(n_sites, window, step)
  seq_names <- rownames(aln)
  if (is.null(seq_names)) seq_names <- as.character(seq_len(nrow(aln)))

  parts <- lapply(seq_len(nrow(aln)), function(i) {
    s <- .roll_from_prefix(
      .seq_prefix(aln[i, ]), b$starts, b$ends,
      window, remove_gaps
    )
    data.frame(
      sequence = seq_names[i], window = seq_along(b$starts),
      start = b$starts, end = b$ends,
      GC = s$gc, GC3 = s$gc3, CpG_OE = s$cpg_oe, entropy = s$entropy,
      gap_fraction = s$gap_fraction,
      GC_pos1 = s$gc_pos1, GC_pos2 = s$gc_pos2, GC_pos3 = s$gc_pos3,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parts)
}
