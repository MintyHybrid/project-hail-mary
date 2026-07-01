# Number of junction-masking residues on each side of a boundary.
# Set to the length of the longest ELM regex match (conservatively 10 aa,
# but 15 covers essentially all ELM motifs).
JUNCTION_MASK_AA <- 15L

# Linker inserted between fragments. Must be:
#   (a) long enough to clearly demarcate boundaries
#   (b) not itself matching common ELM motifs
#   (c) not a valid amino acid run (X = unknown, tolerated by ELM but flagged)
# We use a run of X that is 2 × JUNCTION_MASK_AA so masking from both sides
# never overlaps the real sequence.
LINKER <- strrep("X", 2L * JUNCTION_MASK_AA)
LINKER_LEN <- nchar(LINKER)
# ── 1. Build concatenated sequence + junction map ─────────────────────────────

#' Build a concatenated sequence from a named character vector of fragments.
#'
#' @param fragments Named character vector. Names = isolate / genus IDs.
#'                  Values = amino acid sequences (single-letter, no gaps).
#' @return A list with:
#'   $concat    : single concatenated string (fragments joined by LINKER)
#'   $map       : data frame with one row per fragment:
#'                  isolate, start, stop, frag_len
#'                  (start/stop are 1-based positions in $concat)
build_concat <- function(fragments) {
  stopifnot(is.character(fragments), !is.null(names(fragments)))
  
  # Remove alignment gaps if sequences came straight from a MSA
  fragments <- gsub("[-.]", "", fragments)
  
  # Remove any existing X runs to avoid confusion with our linker
  has_x <- any(grepl("X", fragments, ignore.case = TRUE))
  if (has_x) {
    warning(
      "Some sequences contain 'X' residues. These are retained but may ",
      "interfere with junction detection. Consider replacing with 'A' ",
      "before running if this is a concern."
    )
  }
  
  map_rows <- vector("list", length(fragments))
  cursor <- 1L
  
  for (i in seq_along(fragments)) {
    flen <- nchar(fragments[[i]])
    map_rows[[i]] <- data.frame(
      isolate = names(fragments)[i],
      start = cursor,
      stop = cursor + flen - 1L,
      frag_len = flen,
      stringsAsFactors = FALSE
    )
    cursor <- cursor + flen + LINKER_LEN
  }
  
  concat <- paste(fragments, collapse = LINKER)
  
  list(
    concat = concat,
    map    = bind_rows(map_rows)
  )
}
