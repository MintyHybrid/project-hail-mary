# In-silico assessment of 10-23 DNAzyme catalytic activity for DNAzyme-like
# homologs, based on conservation of the catalytic core.
#
# The 10-23 RNA-cleaving deoxyribozyme has a 15-nt catalytic core,
# 5'-GGCTAGCTACAACGA-3' (Santoro & Joyce 1997, PNAS). The core is highly
# intolerant to substitution (Zaborowska et al. 2002, JBC); position G14 acts as
# the catalytic general base that activates the 2'-OH nucleophile (Borggraefe et
# al. 2022, Nature; Cramer et al. 2023), and the two core thymines (T4, T8)
# contribute functional groups that modulate catalysis. We encode this as a
# per-position criticality weight and score each homolog's core against the
# canonical sequence. This is a sequence-based heuristic, NOT a measurement of
# catalytic rate.

#' Canonical 10-23 DNAzyme catalytic core (15 nt, 5'->3').
#' @export
CANONICAL_1023_CORE <- "GGCTAGCTACAACGA"

# Per-position criticality weight (position 1..15), 0..1; higher = substitution
# more strongly reduces activity. G14 (general base) = 1.0; the two thymines
# (T4, T8) modulatory = 0.80; all other core positions highly conserved = 0.90.
.DNAZYME_CORE_WEIGHTS <- c(0.90, 0.90, 0.90, 0.80, 0.90, 0.90, 0.90, 0.80,
                          0.90, 0.90, 0.90, 0.90, 0.90, 1.00, 0.90)

#' Annotated 10-23 catalytic-core criticality table
#'
#' @return Data frame with one row per core position: position, base,
#'   weight (criticality), tier, and a functional note.
#' @export
dnazyme_core_table <- function() {
  base <- strsplit(CANONICAL_1023_CORE, "")[[1]]
  note <- rep("Conserved core position (substitutions reduce activity)",
              length(base))
  note[c(4, 8)] <- "Core thymine; 5-substituents modulate catalysis"
  note[14]      <- "Catalytic general base (activates 2'-OH nucleophile)"
  data.frame(
    position = seq_along(base),
    base     = base,
    weight   = .DNAZYME_CORE_WEIGHTS,
    tier     = ifelse(.DNAZYME_CORE_WEIGHTS >= 1.0, "essential",
                      ifelse(.DNAZYME_CORE_WEIGHTS >= 0.9, "critical", "modulatory")),
    note     = note,
    stringsAsFactors = FALSE
  )
}

# Best-matching core window (fewest mismatches) within one sequence.
.best_core_window <- function(seq_chr, core) {
  L <- nchar(core); n <- nchar(seq_chr)
  core_v <- strsplit(core, "")[[1]]
  if (n < L) return(list(window = NA_character_, offset = NA_integer_,
                         mism = NA_integer_))
  starts <- seq_len(n - L + 1)
  wins   <- substring(seq_chr, starts, starts + L - 1)
  mism   <- vapply(wins, function(w) sum(strsplit(w, "")[[1]] != core_v), integer(1))
  i <- which.min(mism)
  list(window = wins[i], offset = starts[i], mism = mism[i])
}

#' Score predicted 10-23 DNAzyme activity from core conservation
#'
#' For each sequence, the best-matching 15-nt window is compared to the canonical
#' catalytic core position-by-position. Two summaries are returned:
#' \itemize{
#'   \item \code{weighted_identity}: criticality-weighted fraction of core
#'     positions matching the canonical core (0..1).
#'   \item \code{activity_score}: multiplicative model of retained activity —
#'     the product over substituted positions of \code{(1 - weight)}, so a
#'     substitution at the essential general base (G14, weight 1) predicts an
#'     inactive enzyme (score 0).
#' }
#'
#' @param sequences DNAStringSet or character vector of DNAzyme-like homologs.
#' @param core Canonical core (default \code{CANONICAL_1023_CORE}).
#' @param weights Per-position criticality weights (length 15).
#' @return Data frame: id, core_window, n_substitutions, sub_positions,
#'   critical_substitutions, weighted_identity, activity_score, predicted_active.
#' @export
score_dnazyme_activity <- function(sequences,
                                   core = CANONICAL_1023_CORE,
                                   weights = .DNAZYME_CORE_WEIGHTS) {
  if (inherits(sequences, "DNAStringSet")) {
    ids  <- names(sequences)
    seqs <- as.character(sequences)
  } else {
    seqs <- as.character(sequences)
    ids  <- names(sequences)
  }
  if (is.null(ids)) ids <- paste0("seq_", seq_along(seqs))
  seqs   <- toupper(gsub("U", "T", seqs))
  core   <- toupper(core)
  core_v <- strsplit(core, "")[[1]]
  stopifnot(length(weights) == nchar(core))

  rows <- lapply(seq_along(seqs), function(k) {
    bw <- .best_core_window(seqs[k], core)
    if (is.na(bw$window)) {
      return(data.frame(id = ids[k], core_window = NA_character_,
                        n_substitutions = NA_integer_, sub_positions = NA_character_,
                        critical_substitutions = NA_character_,
                        weighted_identity = NA_real_, activity_score = NA_real_,
                        predicted_active = NA, stringsAsFactors = FALSE))
    }
    win_v   <- strsplit(bw$window, "")[[1]]
    match_i <- win_v == core_v
    sub_pos <- which(!match_i)

    weighted_identity <- sum(weights[match_i]) / sum(weights)
    # retained activity = product of tolerance at substituted positions
    activity_score <- if (length(sub_pos) == 0) 1 else
      prod(1 - weights[sub_pos])
    crit_pos <- sub_pos[weights[sub_pos] >= 0.9]

    data.frame(
      id                    = ids[k],
      core_window           = bw$window,
      n_substitutions       = length(sub_pos),
      sub_positions         = if (length(sub_pos)) paste(sub_pos, collapse = ";") else "",
      critical_substitutions = if (length(crit_pos)) paste(crit_pos, collapse = ";") else "",
      weighted_identity     = round(weighted_identity, 4),
      activity_score        = round(activity_score, 4),
      predicted_active      = activity_score > 0.05 & !any(weights[sub_pos] >= 1.0),
      stringsAsFactors      = FALSE
    )
  })
  do.call(rbind, rows)
}
