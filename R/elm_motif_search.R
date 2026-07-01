# elm_motif_search.R
# ─────────────────────────────────────────────────────────────────────────────
# Batch ELM motif search across concatenated protein sequence fragments.
#
# Designed for: unbiased ELM motif conservation analysis across Potyviridae
# genera, using CI motif-flanking regions (± 20 aa) extracted from the ICTV
# alignment. Handles junction artefact filtering and maps hits back to source
# isolates.
#
# Dependencies: httr, jsonlite, dplyr, stringr, purrr
# ─────────────────────────────────────────────────────────────────────────────

library(httr)
library(jsonlite)
library(dplyr)
library(stringr)
library(purrr)

# ── 0. Constants ──────────────────────────────────────────────────────────────

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

ELM_BASE_URL <- "https://elm.eu.org"

# ── 1. Build concatenated sequence + junction map ─────────────────────────────

#' Build a concatenated sequence from a named character vector of fragments.
#'
#' @param fragments Named character vector. Names = isolate / genus IDs.
#'                  Values = amino acid sequences (single-letter, no gaps).
#' @return A list with:
#'   $concat    : single concatenated string (fragments joined by LINKER)
#'   $map       : data frame with one row per fragment:
#'                  isolate, start, end, frag_len
#'                  (start/end are 1-based positions in $concat)
build_concat <- function(fragments) {
  stopifnot(is.character(fragments), !is.null(names(fragments)))

  # Remove alignment gaps if sequences came straight from a MSA
  fragments <- gsub("[-]", "", fragments)

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
      end = cursor + flen - 1L,
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

#' # ── 2. ELM connectivity check ────────────────────────────────────────────────
#'
#' #' Ping the ELM REST API and verify that key endpoints are reachable.
#' #' !!!!!!!!! NEED TO REPLACE WITH WRAPPER FROM GGET
#' #' Checks three things:
#' #'   1. General reachability of elm.eu.org
#' #'   2. That the /elms.json endpoint returns a valid motif list
#' #'   3. That the /search/sequence/ endpoint accepts a POST
#' #'      (using a minimal 10-aa test sequence)
#' #'
#' #' @param verbose  If TRUE (default), prints a formatted status report.
#' #' @param timeout  Seconds before a check is considered failed.
#' #' @return Invisibly returns a named logical vector:
#' #'   c(reachable, elms_endpoint, search_endpoint)
#' #'   All TRUE = safe to proceed. Any FALSE = do not submit real query.
#' check_elm <- function(verbose = TRUE, timeout = 15L) {
#'   .ok <- "\u2713" # ✓
#'   .fail <- "\u2717" # ✗
#'   .warn <- "\u26a0" # ⚠
#'
#'   status <- c(
#'     reachable = FALSE,
#'     elms_endpoint = FALSE,
#'     search_endpoint = FALSE
#'   )
#'
#'   if (verbose) {
#'     message("── ELM API connectivity check ────────────────────────────────────")
#'     message("  Base URL : ", ELM_BASE_URL)
#'   }
#'
#'   # ── Check 1: General reachability ──────────────────────────────────────────
#'   reach <- tryCatch(
#'     GET(ELM_BASE_URL, timeout(timeout)),
#'     error = function(e) NULL
#'   )
#'
#'   if (!is.null(reach) && status_code(reach) < 500L) {
#'     status["reachable"] <- TRUE
#'     if (verbose) message("  ", .ok, " Host reachable (HTTP ", status_code(reach), ")")
#'   } else {
#'     if (verbose) message("  ", .fail, " Host unreachable — check your internet connection")
#'     if (verbose) message("  Aborting remaining checks.")
#'     return(invisible(status))
#'   }
#'
#'   # ── Check 2: /elms.json ────────────────────────────────────────────────────
#'   elms_url <- paste0(ELM_BASE_URL, "/elms.json")
#'   elms_resp <- tryCatch(
#'     GET(elms_url, timeout(timeout)),
#'     error = function(e) NULL
#'   )
#'
#'   if (!is.null(elms_resp) && status_code(elms_resp) == 200L) {
#'     elms_data <- tryCatch(
#'       fromJSON(content(elms_resp, as = "text", encoding = "UTF-8")),
#'       error = function(e) NULL
#'     )
#'     n_classes <- if (!is.null(elms_data)) {
#'       if (is.data.frame(elms_data)) {
#'         nrow(elms_data)
#'       } else if (is.list(elms_data)) {
#'         length(elms_data)
#'       } else {
#'         NA_integer_
#'       }
#'     } else {
#'       NA_integer_
#'     }
#'
#'     status["elms_endpoint"] <- TRUE
#'     if (verbose) {
#'       msg <- paste0("  ", .ok, " /elms.json OK")
#'       if (!is.na(n_classes)) msg <- paste0(msg, " (", n_classes, " ELM classes loaded)")
#'       message(msg)
#'     }
#'   } else {
#'     http_code <- if (!is.null(elms_resp)) status_code(elms_resp) else "no response"
#'     if (verbose) message("  ", .fail, " /elms.json failed (HTTP ", http_code, ")")
#'   }
#'
#'   # ── Check 3: /search/sequence/ POST ───────────────────────────────────────
#'   # Use a minimal valid peptide — 10 aa, no known pathological motifs,
#'   # short enough that ELM returns quickly
#'   test_seq <- "ACDEFGHIKL"
#'   test_url <- paste0(ELM_BASE_URL, "/search/sequence/")
#'   test_resp <- tryCatch(
#'     POST(test_url,
#'       body   = list(sequence = test_seq, taxon = "Viruses"),
#'       encode = "form",
#'       timeout(timeout)
#'     ),
#'     error = function(e) NULL
#'   )
#'
#'   if (!is.null(test_resp)) {
#'     sc <- status_code(test_resp)
#'     if (sc == 200L) {
#'       # Verify the response actually parses as JSON with expected structure
#'       parsed <- tryCatch(
#'         fromJSON(content(test_resp, as = "text", encoding = "UTF-8")),
#'         error = function(e) NULL
#'       )
#'       if (!is.null(parsed)) {
#'         status["search_endpoint"] <- TRUE
#'         if (verbose) message("  ", .ok, " /search/sequence/ POST OK (response parses as JSON)")
#'       } else {
#'         if (verbose) {
#'           message(
#'             "  ", .warn, " /search/sequence/ returned HTTP 200 but response ",
#'             "is not valid JSON — API may have changed format"
#'           )
#'         }
#'       }
#'     } else if (sc == 405L) {
#'       # 405 = Method Not Allowed — endpoint exists but POST may need different params
#'       if (verbose) {
#'         message(
#'           "  ", .warn, " /search/sequence/ returned HTTP 405 ",
#'           "(Method Not Allowed) — check API docs at ", ELM_BASE_URL, "/api/"
#'         )
#'       }
#'     } else {
#'       if (verbose) message("  ", .fail, " /search/sequence/ returned HTTP ", sc)
#'     }
#'   } else {
#'     if (verbose) message("  ", .fail, " /search/sequence/ POST timed out or connection refused")
#'   }
#'
#'   # ── Summary ────────────────────────────────────────────────────────────────
#'   if (verbose) {
#'     message("─────────────────────────────────────────────────────────────────")
#'     if (all(status)) {
#'       message("  ", .ok, " All checks passed — safe to submit query.")
#'     } else if (status["reachable"] && !status["search_endpoint"]) {
#'       message("  ", .warn, " ELM is reachable but the search endpoint may have moved.")
#'       message("  Check: ", ELM_BASE_URL, "/api/ or ", ELM_BASE_URL, "/swagger/")
#'     } else {
#'       message("  ", .fail, " One or more checks failed — do not submit real query yet.")
#'     }
#'     message("─────────────────────────────────────────────────────────────────")
#'   }
#'
#'   invisible(status)
#' }
#'
#' # ── 3. ELM API query ──────────────────────────────────────────────────────────
#'
#' #' Submit a sequence to ELM and return raw parsed JSON.
#' #'
#' #' @param sequence  Character string. Amino acid sequence.
#' #' @param taxon     Taxonomic context string passed to ELM.
#' #'                  Use "Viruses" for viral sequences.
#' #'                  Use "Viridiplantae" for plant proteins.
#' #' @param timeout   Seconds before giving up (ELM can be slow).
#' #' @return Parsed list from ELM JSON response, or NULL on failure.
#' query_elm <- function(sequence,
#'                       taxon = "Viruses",
#'                       timeout = 120L) {
#'   url <- paste0(ELM_BASE_URL, "/search/sequence/")
#'
#'   resp <- tryCatch(
#'     POST(url,
#'       body   = list(sequence = sequence, taxon = taxon),
#'       encode = "form",
#'       timeout(timeout)
#'     ),
#'     error = function(e) {
#'       message("ELM query failed: ", conditionMessage(e))
#'       NULL
#'     }
#'   )
#'
#'   if (is.null(resp)) {
#'     return(NULL)
#'   }
#'
#'   status <- status_code(resp)
#'   if (status != 200L) {
#'     message("ELM returned HTTP ", status, " — check sequence and taxon.")
#'     return(NULL)
#'   }
#'
#'   raw <- content(resp, as = "text", encoding = "UTF-8")
#'
#'   tryCatch(
#'     fromJSON(raw, simplifyDataFrame = TRUE),
#'     error = function(e) {
#'       message("JSON parse error: ", conditionMessage(e))
#'       NULL
#'     }
#'   )
#' }

# ── 4. Parse ELM hits into a tidy data frame ──────────────────────────────────

#' Extract hit table from ELM response object.
#'
#' @param elm_result  List returned by query_elm().
#' @return Data frame with columns:
#'   elm_identifier, elm_name, start, end, sequence_match,
#'   probability, logic, phiblast_match (where available)
parse_elm_hits <- function(elm_result) {
  if (is.null(elm_result)) {
    return(tibble())
  }

  # ELM response structure varies slightly by version — handle both
  hits <- elm_result$instances %||% elm_result$results %||% elm_result$Instances %||% NULL

  if (is.null(hits) || length(hits) == 0) {
    message("No hits returned by ELM.")
    return(tibble())
  }

  # Coerce to data frame if it came back as a list of lists
  if (is.data.frame(hits)) {
    df <- as_tibble(hits)
  } else {
    df <- map_dfr(hits, as_tibble)
  }

  # Normalise column names to snake_case
  names(df) <- tolower(gsub("([A-Z])", "_\\1", names(df)))
  names(df) <- gsub("^_", "", names(df))
  names(df) <- gsub("__", "_", names(df))

  # Ensure start/end are integer
  if ("start" %in% colnames(df)) df$start <- as.integer(df$start)
  if ("end" %in% colnames(df)) df$end <- as.integer(df$end)

  df
}

# ── 5. Filter junction-spanning hits ─────────────────────────────────────────

#' Remove hits that overlap junction (linker) regions.
#'
#' A hit is flagged if it overlaps any linker region by at least 1 residue,
#' OR if either endpoint falls within JUNCTION_MASK_AA residues of a boundary.
#'
#' @param hits   Data frame from parse_elm_hits().
#' @param jmap   $map data frame from build_concat().
#' @return hits with an added logical column `junction_artefact`;
#'         rows where junction_artefact == TRUE are kept but flagged
#'         so the caller can inspect or drop them.
flag_junction_hits <- function(hits, jmap) {
  if (nrow(hits) == 0L) {
    return(hits)
  }

  # Build a list of linker intervals [linker_start, linker_end]
  # Linker follows each fragment except the last
  n_frags <- nrow(jmap)
  linker_intervals <- if (n_frags > 1L) {
    map_dfr(seq_len(n_frags - 1L), function(i) {
      tibble(
        lstart = jmap$end[i] + 1L,
        lend   = jmap$start[i + 1L] - 1L
      )
    })
  } else {
    tibble(lstart = integer(0), lend = integer(0))
  }

  # Masking zone: JUNCTION_MASK_AA residues into real sequence on each side
  mask_intervals <- map_dfr(seq_len(n_frags - 1L), function(i) {
    bind_rows(
      tibble(
        lstart = jmap$end[i] - JUNCTION_MASK_AA + 1L,
        lend = jmap$end[i]
      ),
      tibble(
        lstart = jmap$start[i + 1L],
        lend = jmap$start[i + 1L] + JUNCTION_MASK_AA - 1L
      )
    )
  })

  danger_zones <- bind_rows(linker_intervals, mask_intervals)

  hits$junction_artefact <- map2_lgl(hits$start, hits$end, function(hs, he) {
    any(hs <= danger_zones$lend & he >= danger_zones$lstart)
  })

  hits
}

# ── 5. Map hits back to source isolates ───────────────────────────────────────

#' Map hit positions in the concatenated sequence back to source isolates.
#'
#' @param hits   Data frame with junction_artefact column (from flag_junction_hits).
#' @param jmap   $map data frame from build_concat().
#' @return hits with added columns:
#'   isolate, local_start, local_end
#'   (NA for junction artefacts)
map_hits_to_isolates <- function(hits, jmap) {
  if (nrow(hits) == 0L) {
    return(hits)
  }

  assign_isolate <- function(pos) {
    idx <- which(jmap$start <= pos & jmap$end >= pos)
    if (length(idx) == 0L) {
      return(list(isolate = NA_character_, local = NA_integer_))
    }
    list(
      isolate = jmap$isolate[idx[1L]],
      local = pos - jmap$start[idx[1L]] + 1L
    )
  }

  hits$isolate <- NA_character_
  hits$local_start <- NA_integer_
  hits$local_end <- NA_integer_

  clean <- !hits$junction_artefact
  if (any(clean)) {
    start_info <- map(hits$start[clean], assign_isolate)
    end_info <- map(hits$end[clean], assign_isolate)

    hits$isolate[clean] <- map_chr(start_info, "isolate")
    hits$local_start[clean] <- map_int(start_info, "local")
    hits$local_end[clean] <- map_int(end_info, "local")
  }

  hits
}

# ── 6. Conservation summary ───────────────────────────────────────────────────

#' Summarise motif hit frequency and conservation across isolates / genera.
#'
#' @param hits       Data frame from map_hits_to_isolates().
#' @param meta       Optional data frame with columns `isolate` and `genus`
#'                   for genus-level aggregation. If NULL, only isolate-level
#'                   stats are computed.
#' @param n_total    Total number of input isolates (for frequency calculation).
#' @return List with:
#'   $by_motif   : frequency table per ELM class
#'   $by_genus   : genus-level presence/absence (if meta supplied)
#'   $hit_detail : full cleaned hit table
summarise_conservation <- function(hits, meta = NULL, n_total) {
  clean_hits <- filter(hits, !junction_artefact)

  by_motif <- clean_hits %>%
    group_by(elm_identifier) %>%
    summarise(
      n_hits      = n(),
      n_isolates  = n_distinct(isolate),
      freq        = n_distinct(isolate) / n_total,
      .groups     = "drop"
    ) %>%
    arrange(desc(freq))

  by_genus <- NULL
  if (!is.null(meta)) {
    stopifnot("isolate" %in% names(meta), "genus" %in% names(meta))
    clean_hits <- left_join(clean_hits, meta, by = "isolate")

    by_genus <- clean_hits %>%
      group_by(elm_identifier, genus) %>%
      summarise(
        n_isolates = n_distinct(isolate),
        .groups    = "drop"
      ) %>%
      tidyr::pivot_wider(
        names_from = genus,
        values_from = n_isolates,
        values_fill = 0L
      )
  }

  list(
    by_motif   = by_motif,
    by_genus   = by_genus,
    hit_detail = clean_hits
  )
}

# ── 7. Master wrapper ─────────────────────────────────────────────────────────

#' Run the full ELM conservation pipeline.
#'
#' @param fragments  Named character vector. Names = isolate IDs, values = AA seqs.
#' @param taxon      ELM taxonomic context (default "Viruses").
#' @param meta       Optional data frame with `isolate` + `genus` for genus summary.
#' @param drop_artefacts  If TRUE (default), junction artefacts are removed from
#'                        the returned hit table. Set FALSE to inspect them.
#' @return List with $hits, $map, $summary (from summarise_conservation()).
run_elm_conservation <- function(fragments,
                                 taxon = "Viruses",
                                 meta = NULL,
                                 drop_artefacts = TRUE) {
  message("── Step 1/5: Building concatenated sequence ──────────────────────")
  cc <- build_concat(fragments)
  concat <- cc$concat
  jmap <- cc$map

  message("  Fragments  : ", nrow(jmap))
  message("  Total length: ", nchar(concat), " aa (incl. linkers)")
  message("  Linker length: ", LINKER_LEN, " X residues")

  message("── Step 2/5: Querying ELM ───────────────────────────────────────")
  message("  Taxon: ", taxon)
  elm_raw <- query_elm(concat, taxon = taxon)

  if (is.null(elm_raw)) {
    stop("ELM query returned no result. Check connectivity and taxon string.")
  }

  message("── Step 3/5: Parsing hits ────────────────────────────────────────")
  hits <- parse_elm_hits(elm_raw)
  message("  Raw hits: ", nrow(hits))

  message("── Step 4/5: Flagging junction artefacts ────────────────────────")
  hits <- flag_junction_hits(hits, jmap)
  n_artefacts <- sum(hits$junction_artefact)
  message(
    "  Junction artefacts flagged: ", n_artefacts,
    " (", round(100 * n_artefacts / max(nrow(hits), 1L), 1), "%)"
  )

  if (drop_artefacts) hits <- filter(hits, !junction_artefact)

  message("── Step 5/5: Mapping hits to isolates ───────────────────────────")
  hits <- map_hits_to_isolates(hits, jmap)

  summary <- summarise_conservation(
    hits    = hits,
    meta    = meta,
    n_total = nrow(jmap)
  )

  message("── Done ─────────────────────────────────────────────────────────")
  message("  Clean hits     : ", nrow(hits))
  message("  Unique motifs  : ", n_distinct(hits$elm_identifier))
  message(
    "  Isolates with hits: ",
    n_distinct(hits$isolate[!is.na(hits$isolate)]), " / ", nrow(jmap)
  )

  list(
    hits    = hits,
    map     = jmap,
    summary = summary,
    elm_raw = elm_raw # keep raw response for debugging
  )
}

# ── 8. Usage example ──────────────────────────────────────────────────────────

if (FALSE) {
  # fragments: named vector of ± 20 aa windows around your CI motif
  # (one entry per isolate, gaps already removed)
  fragments <- c(
    "PVY_isolate_01"    = "MASVKDTEMASVKDTEMSVKDTEMSVKDTEMSVKDTEM",
    "TuMV_isolate_01"   = "MDTVKSTEMASVKDTEMSVKDTEMSVKDTEMSVKDTEM",
    "SCMV_isolate_01"   = "MKTVKDTEMCSVKDTEMSVKDTEMSVKDTEMSVKDTEM"
    # ... all ~200 isolates
  )

  # Optional genus metadata for genus-level conservation summary
  meta <- data.frame(
    isolate = c("PVY_isolate_01", "TuMV_isolate_01", "SCMV_isolate_01"),
    genus   = c("Potyvirus", "Potyvirus", "Macluravirus")
  )

  result <- run_elm_conservation(
    fragments      = fragments,
    taxon          = "Viruses",
    meta           = meta,
    drop_artefacts = TRUE
  )

  # Top conserved motifs
  print(result$summary$by_motif)

  # Genus-level presence/absence matrix
  print(result$summary$by_genus)

  # Full hit detail
  View(result$hits)

  # Focus on tyrosine-based sorting motifs specifically
  tyrosine_motifs <- filter(
    result$hits,
    grepl("TRG_Yxx|LIG_AP4|TRG_DiLeu|LIG_SH2_CRK", elm_identifier)
  )
  print(tyrosine_motifs)
}
