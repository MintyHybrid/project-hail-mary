library(httr2)

# ── Constants ────────────────────────────────────────────────────────────────
ELM_BASE_URL  <- "http://elm.eu.org"
ELM_ENDPOINT  <- "start_search"
MAX_AA_LENGTH <- 2000L
RATE_LIMIT_S  <- 61L        # 1 request/min + 1s buffer
LINKER        <- strrep("X", 15)

# ── TSV column schema from ELM API docs ──────────────────────────────────────
ELM_TSV_COLS <- c(
  "elm_identifier",
  "start",
  "stop",
  "is_annotated",
  "is_phiblastmatch",
  "is_filtered",
  "structure",
  "topodomfilter",
  "taxonfilter"
)

# =============================================================================
# PART 1 — SPLITTING: decompose concatenated sequence at linker boundaries
# =============================================================================

#' Split a concatenated sequence into chunks at XXXXXXXXXXXXXXX linker sites,
#' keeping chunks within the 2000 aa ELM limit.
#'
#' Splitting strategy:
#'   - Find all linker positions in the concatenated string
#'   - Greedily accumulate individual sequences into chunks, starting a new
#'     chunk whenever the next sequence would push the chunk over MAX_AA_LENGTH
#'   - Each chunk therefore contains one or more complete original sequences
#'     separated by intact linkers — no sequence is ever split mid-residue
#'
#' @param concat_seq   Character string. The full concatenated AA sequence.
#' @param sample_df    Data frame with columns: isolate, start, stop, frag_len
#'                     (1-based positions in concat_seq).
#'
#' @return A list of chunk descriptors, each containing:
#'   $chunk_id      Character label
#'   $sequence      The chunk's AA sequence (substring of concat_seq)
#'   $offset        Integer. concat_seq position of chunk's first residue (1-based)
#'   $samples       Subset of sample_df for sequences in this chunk
elm_split_at_linkers <- function(concat_seq, sample_df) {

  stopifnot(
    is.character(concat_seq), length(concat_seq) == 1,
    is.data.frame(sample_df),
    all(c("isolate", "start", "end", "frag_len") %in% colnames(sample_df))
  )

  n_samples <- nrow(sample_df)

  # ── Greedily bin samples into chunks ─────────────────────────────────────
  chunks    <- list()
  chunk_idx <- 1L
  bin_start <- 1L   # index into sample_df rows

  while (bin_start <= n_samples) {

    bin_stop <- bin_start

    # Extstop bin as far as possible without exceeding MAX_AA_LENGTH
    # Chunk spans from start of first sample to stop of last sample in bin,
    # inclusive of any linkers between them.
    while (bin_stop < n_samples) {
      proposed_stop <- sample_df$end[bin_stop + 1L]
      # Also include the linker between current last and next sample
      # (linker sits between stop[i] and start[i+1])
      chunk_len <- proposed_stop - sample_df$start[bin_start] + 1L

      if (chunk_len > MAX_AA_LENGTH) break
      bin_stop <- bin_stop + 1L
    }

    # Safety: single sequence itself exceeds limit — warn and include anyway
    # (ELM will truncate, but that's a data issue the caller must handle upstream)
    chunk_start_pos <- sample_df$start[bin_start]
    chunk_stop_pos  <- sample_df$end[bin_stop]
    chunk_seq       <- substr(concat_seq, chunk_start_pos, chunk_stop_pos)
    chunk_len       <- nchar(chunk_seq)

    if (chunk_len > MAX_AA_LENGTH) {
      warning(sprintf(
        "Chunk %d contains a single sequence (%s) that is %d aa > %d aa limit. ",
        chunk_idx, sample_df$sample_name[bin_start], chunk_len, MAX_AA_LENGTH,
        "ELM will truncate — consider filtering this sequence before querying."
      ))
    }

    chunks[[chunk_idx]] <- list(
      chunk_id  = sprintf("chunk_%03d", chunk_idx),
      sequence  = chunk_seq,
      offset    = chunk_start_pos,   # 1-based position in original concat_seq
      samples   = sample_df[bin_start:bin_stop, , drop = FALSE]
    )

    chunk_idx <- chunk_idx + 1L
    bin_start <- bin_stop + 1L
  }

  message(sprintf(
    "Split %d sequences into %d chunk(s) (max %d aa each).",
    n_samples, length(chunks), MAX_AA_LENGTH
  ))

  chunks
}


# =============================================================================
# PART 2 — QUERYING: GET request per chunk, timed, rate-limited
# =============================================================================

#' Query ELM for a single sequence chunk
#'
#' @param chunk   A chunk descriptor from elm_split_at_linkers()
#' @return A list: $chunk_id, $offset, $samples, $elapsed_s, $status, $hits (data frame or NULL)
elm_query_chunk <- function(chunk) {

  req <- request(ELM_BASE_URL) |>
    req_url_path_append(ELM_ENDPOINT) |>
    req_url_path_append(chunk$sequence) |>
    req_headers("Accept" = "text/plain") |>
    req_timeout(120) |>
    req_retry(
      max_tries = 3,
      retry_on_failure = TRUE,
      backoff = \(attempt) 30 * attempt   # 30s, 60s
    )

  t_start  <- proc.time()[["elapsed"]]
  response <- tryCatch(
    req_perform(req),
    error = function(e) {
      message(sprintf("[%s] Request failed: %s", chunk$chunk_id, conditionMessage(e)))
      NULL
    }
  )
  elapsed <- round(proc.time()[["elapsed"]] - t_start, 2)

  if (is.null(response)) {
    return(list(
      chunk_id = chunk$chunk_id, offset = chunk$offset,
      samples = chunk$samples, elapsed_s = elapsed,
      status = NA_integer_, hits = NULL
    ))
  }

  status <- resp_status(response)
  hits   <- tryCatch(
    elm_parse_tsv(resp_body_string(response)),
    error = function(e) {
      message(sprintf("[%s] TSV parse error: %s", chunk$chunk_id, conditionMessage(e)))
      NULL
    }
  )

  message(sprintf("[%s] HTTP %d | %.1fs | %d aa | %d hits",
    chunk$chunk_id, status, elapsed,
    nchar(chunk$sequence), if (is.null(hits)) 0L else nrow(hits)
  ))

  list(
    chunk_id  = chunk$chunk_id,
    offset    = chunk$offset,
    samples   = chunk$samples,
    elapsed_s = elapsed,
    status    = status,
    hits      = hits
  )
}


#' Parse ELM TSV response body into a data frame
#'
#' @param tsv_text  Raw character string from resp_body_string()
#' @return Data frame with ELM_TSV_COLS columns, or empty data frame if no hits
elm_parse_tsv <- function(tsv_text) {

  writeLines(head(strsplit(tsv_text, "\n", fixed = TRUE)[[1]], 20), con = stderr())
  lines <- strsplit(tsv_text, "\n", fixed = TRUE)[[1]]

  # Drop blank lines and lines starting with #
  lines <- lines[nzchar(trimws(lines)) & !startsWith(lines, "#")]

  if (length(lines) == 0L) {
    return(elm_empty_hits())
  }

  # ── Detect header row ─────────────────────────────────────────────────────
  # ELM may return a header line containing "elm_identifier" or starting with
  # a known column name. Use it if present; fall back to hardcoded schema.
  first_fields <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]

  if (tolower(first_fields[1]) %in% c("elm_identifier", "elm identifier", "#elm_identifier")) {
    # Header present — parse col names from it, drop the row from data
    col_names <- tolower(gsub("[^a-zA-Z0-9_]", "_", trimws(first_fields)))
    data_lines <- lines[-1]
  } else {
    col_names  <- ELM_TSV_COLS
    data_lines <- lines
  }

  if (length(data_lines) == 0L) return(elm_empty_hits())

  # ── Count actual columns in data to catch schema mismatches early ─────────
  n_cols_data   <- lengths(strsplit(data_lines[1], "\t", fixed = TRUE))
  n_cols_schema <- length(col_names)

  if (n_cols_data != n_cols_schema) {
    warning(sprintf(
      "ELM TSV column count mismatch: got %d, expected %d. ",
      n_cols_data, n_cols_schema,
      "Using positional assignment; check ELM_TSV_COLS against current API schema."
    ))
    # Truncate or pad col_names to match actual data width
    if (n_cols_data > n_cols_schema) {
      col_names <- c(col_names, paste0("extra_col_", seq_len(n_cols_data - n_cols_schema)))
    } else {
      col_names <- col_names[seq_len(n_cols_data)]
    }
  }

  df <- read.table(
    text             = paste(data_lines, collapse = "\n"),
    sep              = "\t",
    header           = FALSE,
    col.names        = col_names,
    quote            = "",
    comment.char     = "",
    stringsAsFactors = FALSE,
    fill             = TRUE    # tolerate ragged rows rather than erroring
  )

  df$start <- as.integer(df$start)
  df$stop  <- as.integer(df$stop)

  df
}

elm_empty_hits <- function() {
  df        <- as.data.frame(setNames(lapply(ELM_TSV_COLS, \(.) character(0)), ELM_TSV_COLS))
  df$start  <- integer(0)
  df$stop   <- integer(0)
  df
}

# =============================================================================
# PART 3 — RECONSTITUTION: lift chunk-local coordinates back to concat space,
#           flag linker-overlapping hits, join sample names
# =============================================================================

LINKER_LEN <- nchar(LINKER)  # 15

#' Reconstitute results from all queried chunks into a single annotated data frame
#'
#' @param query_results  List of results from elm_query_chunk()
#' @param sample_df      Original sample data frame (sample_name, start, stop)
#' @param concat_seq     The original concatenated sequence string
#'
#' @return Data frame with all hits, absolute coordinates, sample assignments,
#'         and linker-overlap flags
elm_reconstitute <- function(query_results, sample_df, concat_seq) {
  # ── Precompute linker intervals in concat_seq ─────────────────────────────
  linker_intervals <- elm_find_linkers(concat_seq)

  all_hits <- lapply(query_results, function(res) {

    if (is.null(res$hits) || nrow(res$hits) == 0L) return(NULL)

    hits <- res$hits

    # Lift coordinates: chunk-local (1-based) → concat_seq absolute (1-based)
    # chunk$offset is the 1-based start of this chunk in concat_seq
    hits$abs_start <- hits$start + res$offset - 1L
    hits$abs_stop  <- hits$stop  + res$offset - 1L
    hits$chunk_id  <- res$chunk_id

    # ── Flag hits overlapping any linker ────────────────────────────────────
    hits$in_linker <- mapply(
      function(s, e) elm_overlaps_linker(s, e, linker_intervals),
      hits$abs_start, hits$abs_stop
    )

    hits
  })

  all_hits <- do.call(rbind, Filter(Negate(is.null), all_hits))

  if (is.null(all_hits) || nrow(all_hits) == 0L) {
    message("No hits returned across all chunks.")
    return(invisible(NULL))
  }

  # ── Assign sample names by overlap with sample_df intervals ──────────────
  # A hit belongs to the sample whose [start, stop] interval contains abs_start
  # (using start as the anchor; hits spanning a linker are already flagged)
  all_hits$sample_name <- vapply(all_hits$abs_start, function(pos) {
    idx <- which(sample_df$start <= pos & sample_df$end >= pos)
    if (length(idx) == 0L) NA_character_
    else sample_df$sample_name[idx[1L]]
  }, character(1L))

  # ── Compute within-sample local coordinates ───────────────────────────────
  all_hits$local_start <- all_hits$abs_start - sample_df$start[
    match(all_hits$sample_name, sample_df$sample_name)
  ] + 1L
  all_hits$local_stop  <- all_hits$abs_stop - sample_df$end[
    match(all_hits$sample_name, sample_df$sample_name)
  ] + 1L

  # ── Reorder columns for readability ──────────────────────────────────────
  first_cols <- c("sample_name", "elm_identifier",
    "local_start", "local_stop",
    "abs_start",   "abs_stop",
    "in_linker",   "chunk_id",
    "is_annotated", "is_phiblastmatch", "is_filtered",
    "structure",    "topodomfilter",    "taxonfilter")
  remaining  <- setdiff(names(all_hits), first_cols)
  all_hits[, c(first_cols, remaining)]
}


#' Find all linker intervals in a concatenated sequence
elm_find_linkers <- function(concat_seq) {
  starts <- gregexpr(LINKER, concat_seq, fixed = TRUE)[[1]]
  if (identical(as.integer(starts), -1L)) return(data.frame(start = integer(0), stop = integer(0)))
  data.frame(
    start = as.integer(starts),
    stop  = as.integer(starts) + LINKER_LEN - 1L
  )
}

#' Test whether a hit [s, e] overlaps any linker interval
elm_overlaps_linker <- function(s, e, linker_intervals) {
  if (nrow(linker_intervals) == 0L) return(FALSE)
  any(s <= linker_intervals$stop & e >= linker_intervals$start)
}


# =============================================================================
# PART 4 — ORCHESTRATOR: the single entry point
# =============================================================================

#' Query ELM for a concatenated sequence, splitting at linkers, with rate limiting
#'
#' @param elm_input   A list with two elements:
#'                    [[1]] character: the concatenated AA sequence
#'                    [[2]] data.frame: sample metadata with columns
#'                          sample_name, start, stop (1-based, in concat seq)
#'
#' @return Data frame of reconstituted, annotated ELM hits
elm_query_concatenated <- function(elm_input) {

  stopifnot(
    is.list(elm_input), length(elm_input) >= 2L,
    is.character(elm_input[[1]]),
    is.data.frame(elm_input[[2]])
  )

  concat_seq <- elm_input[[1]]
  sample_df  <- elm_input[[2]]

  # ── 1. Split ──────────────────────────────────────────────────────────────
  chunks <- elm_split_at_linkers(concat_seq, sample_df)
  n      <- length(chunks)

  # ── 2. Query with rate limiting ───────────────────────────────────────────
  results <- vector("list", n)

  for (i in seq_len(n)) {
    message(sprintf("\n── Chunk %d / %d: %s (%d aa, %d sequences) ──",
      i, n, chunks[[i]]$chunk_id,
      nchar(chunks[[i]]$sequence),
      nrow(chunks[[i]]$samples)
    ))

    results[[i]] <- elm_query_chunk(chunks[[i]])

    if (i < n) {
      message(sprintf("  Waiting %ds (rate limit)...", RATE_LIMIT_S))
      Sys.sleep(RATE_LIMIT_S)
    }
  }

  # ── 3. Reconstitute ───────────────────────────────────────────────────────
  elm_reconstitute(results, sample_df, concat_seq)
}
