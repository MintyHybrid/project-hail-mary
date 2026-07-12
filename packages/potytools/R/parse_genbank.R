# =============================================================================
# GenBank File Parser - Base R, no external dependencies
# =============================================================================
# Parses all .gb/.gbk/.genbank files in a folder.
#
# Output structure (list of records, one per file or per // record):
#   $name        : LOCUS name
#   $accession   : primary accession
#   $definition  : DEFINITION field (collapsed)
#   $mol_type    : molecule type from LOCUS (e.g. DNA, RNA, mRNA)
#   $topology    : "circular" or "linear"
#   $organism    : organism name (SOURCE > ORGANISM)
#   $length      : declared sequence length (bp)
#   $features    : list of feature entries, each with:
#       $type        : feature key (gene, CDS, rRNA, tRNA, misc_RNA, ...)
#       $location    : raw location string (join/complement preserved)
#       $strand      : 1 (forward) or -1 (reverse/complement)
#       $start       : minimum nucleotide position (1-based)
#       $end         : maximum nucleotide position (1-based)
#       $qualifiers  : named list of qualifier values; multi-value -> character vector
#                      /translation is stored as a clean string (protein seq)
#   $sequence    : full nucleotide sequence (lowercase string) or NA
#   $proteins    : named list of CDS protein sequences keyed by locus_tag or product
# =============================================================================

# ---- Main loader ------------------------------------------------------------

#' Load all GenBank files from a folder
#'
#' @param folder  Path to folder containing .gb / .gbk / .genbank files.
#' @param pattern Regex to match filenames (default catches .gb, .gbk, .genbank).
#' @param recursive Search subdirectories? (default FALSE)
#' @return Named list; names are file basenames without extension.
#'         Multi-record files get names like "filename__1", "filename__2".
#' @export
load_genbank_folder <- function(folder,
                                pattern = "\\.(gb|gbk|genbank)$",
                                recursive = FALSE) {
  files <- list.files(folder,
    pattern = pattern,
    full.names = TRUE,
    ignore.case = TRUE,
    recursive = recursive
  )

  if (length(files) == 0) {
    warning("No GenBank files found in: ", folder)
    return(list())
  }

  message("Found ", length(files), " GenBank file(s). Parsing...")

  all_records <- list()
  for (f in files) {
    base <- tools::file_path_sans_ext(basename(f))
    records <- parse_genbank_file(f)
    if (length(records) == 1) {
      all_records[[base]] <- records[[1]]
    } else {
      for (i in seq_along(records)) {
        key <- paste0(base, "__", i)
        all_records[[key]] <- records[[i]]
      }
    }
  }

  message("Loaded ", length(all_records), " record(s) total.")
  all_records
}

# ---- File-level parser ------------------------------------------------------

#' Parse a single GenBank file (may contain multiple // records)
#'
#' @param filepath Path to a .gb file.
#' @return List of records.
#' @export
parse_genbank_file <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)

  # Split on "//" terminators to get individual records
  terminators <- which(grepl("^//$", lines))
  if (length(terminators) == 0) {
    # No terminator - try parsing the whole file as one record
    return(list(.parse_record(lines)))
  }

  records <- list()
  start <- 1L
  for (term in terminators) {
    block <- lines[start:(term - 1L)]
    if (any(nchar(trimws(block)) > 0)) {
      records[[length(records) + 1L]] <- .parse_record(block)
    }
    start <- term + 1L
  }
  records
}

# ---- Record parser ----------------------------------------------------------

.parse_record <- function(lines) {
  # ---- Helper: collect a section that may span indented continuation lines --
  .collect_field <- function(start_keyword) {
    idx <- grep(paste0("^", start_keyword, "\\b"), lines)
    if (length(idx) == 0) {
      return(NA_character_)
    }
    i <- idx[1]
    text <- sub(paste0("^", start_keyword, "\\s*"), "", lines[i])
    i <- i + 1L
    while (i <= length(lines) &&
      grepl("^\\s+", lines[i]) &&
      !grepl("^[A-Z]{2,}\\s", lines[i])) {
      text <- paste(text, trimws(lines[i]))
      i <- i + 1L
    }
    trimws(text)
  }

  # ---- LOCUS ------------------------------------------------------------------
  locus_line <- lines[grep("^LOCUS\\b", lines)[1]]
  locus_tok <- strsplit(trimws(locus_line), "\\s+")[[1]]
  name <- if (length(locus_tok) >= 2) locus_tok[2] else NA_character_
  seq_length <- suppressWarnings(as.integer(locus_tok[3]))
  # mol_type: typically token 5 (DNA/RNA/mRNA etc.)
  mol_type <- NA_character_
  topology <- "linear"
  for (tok in locus_tok[-c(1, 2, 3)]) {
    if (grepl("circular|linear", tok, ignore.case = TRUE)) {
      topology <- tolower(tok)
    } else if (grepl("^(ss-|ds-|ms-)?(DNA|RNA|mRNA|rRNA|tRNA|cRNA|uRNA)$",
      tok,
      ignore.case = TRUE
    )) {
      mol_type <- toupper(tok)
    }
  }

  # ---- Scalar metadata fields ------------------------------------------------
  accession <- .collect_field("ACCESSION")
  accession <- strsplit(accession, "\\s+")[[1]][1] # first token only
  definition <- .collect_field("DEFINITION")
  version <- .collect_field("VERSION")

  # Organism: under SOURCE/ORGANISM
  org_idx <- grep("^  ORGANISM\\b", lines)
  organism <- if (length(org_idx) > 0) {
    trimws(sub("^\\s*ORGANISM\\s*", "", lines[org_idx[1]]))
  } else {
    NA_character_
  }

  # ---- Feature table ---------------------------------------------------------
  feat_start <- grep("^FEATURES\\b", lines)[1]
  orig_start <- grep("^ORIGIN\\b", lines)[1]

  features <- list()
  if (!is.na(feat_start)) {
    feat_end <- if (!is.na(orig_start)) orig_start - 1L else length(lines)
    feat_block <- lines[(feat_start + 1L):feat_end]
    features <- .parse_features(feat_block)
  }

  # ---- Nucleotide sequence ---------------------------------------------------
  sequence <- NA_character_
  if (!is.na(orig_start) && orig_start < length(lines)) {
    seq_lines <- lines[(orig_start + 1L):length(lines)]
    seq_lines <- seq_lines[!grepl("^\\s*$", seq_lines)]
    # Strip line numbers and spaces, keep only [acgturyswkmbdhvn] characters
    seq_clean <- gsub("[^acgturyswkmbdhvnACGTURYSWKMBDHVN]", "", seq_lines)
    sequence <- paste(seq_clean, collapse = "")
    sequence <- tolower(sequence)
  }

  # ---- Convenience: extract all protein sequences from CDS ------------------
  proteins <- .extract_proteins(features)

  list(
    name       = name,
    accession  = accession,
    definition = definition,
    version    = version,
    mol_type   = mol_type,
    topology   = topology,
    organism   = organism,
    length     = seq_length,
    features   = features,
    sequence   = sequence,
    proteins   = proteins
  )
}

# ---- Feature block parser ---------------------------------------------------

.parse_features <- function(feat_lines) {
  # Feature header lines: 5 leading spaces then a non-space character
  # Qualifier/continuation lines: 21 leading spaces
  # We detect feature starts by a change in indentation pattern.

  feat_start_idx <- grep("^     \\S", feat_lines)
  if (length(feat_start_idx) == 0) {
    return(list())
  }

  features <- vector("list", length(feat_start_idx))

  for (i in seq_along(feat_start_idx)) {
    s <- feat_start_idx[i]
    e <- if (i < length(feat_start_idx)) {
      feat_start_idx[i + 1L] - 1L
    } else {
      length(feat_lines)
    }
    blk <- feat_lines[s:e]

    # --- Feature type + location (may wrap onto next lines) ---
    header <- blk[1]
    feat_type <- trimws(substr(header, 1, 20)) # cols 1-20 (1-based)
    loc_part <- trimws(substr(header, 21, nchar(header)))

    # Location continuation: lines that start with 21 spaces but don't have /
    j <- 2L
    while (j <= length(blk) &&
      grepl("^                     [^/]", blk[j])) {
      loc_part <- paste0(loc_part, trimws(blk[j]))
      j <- j + 1L
    }
    location <- loc_part

    # --- Parse strand and numeric range from location string ---
    is_complement <- grepl("complement", location, ignore.case = TRUE)
    strand <- if (is_complement) -1L else 1L
    nums <- as.integer(regmatches(
      location,
      gregexpr("[0-9]+", location)
    )[[1]])
    feat_start_pos <- if (length(nums) > 0) min(nums) else NA_integer_
    feat_end_pos <- if (length(nums) > 0) max(nums) else NA_integer_

    # --- Qualifier lines ---
    qual_lines <- if (j <= length(blk)) blk[j:length(blk)] else character(0)
    qualifiers <- .parse_qualifiers(qual_lines)

    features[[i]] <- list(
      type       = feat_type,
      location   = location,
      strand     = strand,
      start      = feat_start_pos,
      end        = feat_end_pos,
      qualifiers = qualifiers
    )
  }

  features
}

# ---- Qualifier parser -------------------------------------------------------

.parse_qualifiers <- function(lines) {
  if (length(lines) == 0) {
    return(list())
  }

  # Glue continuation lines: a qualifier starts with /key or /key="
  # We walk through lines and accumulate each qualifier's text.
  chunks <- list()
  current <- NULL

  for (ln in lines) {
    stripped <- trimws(ln)
    if (startsWith(stripped, "/")) {
      if (!is.null(current)) chunks[[length(chunks) + 1L]] <- current
      current <- stripped
    } else {
      # Continuation of previous qualifier value
      current <- paste0(current, stripped)
    }
  }
  if (!is.null(current)) chunks[[length(chunks) + 1L]] <- current

  quals <- list()

  for (chunk in chunks) {
    if (!startsWith(chunk, "/")) next

    if (grepl("=", chunk)) {
      key <- sub("^/([A-Za-z_][A-Za-z0-9_]*)=.*", "\\1", chunk)
      val <- sub('^/[A-Za-z_][A-Za-z0-9_]*="?(.*?)"?$', "\\1", chunk, perl = TRUE)
      # Strip surrounding quotes explicitly
      if (startsWith(val, '"') && endsWith(val, '"')) {
        val <- substr(val, 2L, nchar(val) - 1L)
      }
      # Protein sequences have no internal whitespace - strip any artefacts
      if (key == "translation") val <- gsub("\\s", "", val)
    } else {
      # Flag qualifier: /pseudo, /partial, etc.
      key <- sub("^/([A-Za-z_][A-Za-z0-9_]*).*", "\\1", chunk)
      val <- TRUE
    }

    # Allow multi-value qualifiers (e.g. multiple /db_xref)
    if (key %in% names(quals)) {
      quals[[key]] <- c(quals[[key]], val)
    } else {
      quals[[key]] <- val
    }
  }

  quals
}

# ---- Protein extraction helper ----------------------------------------------

.extract_proteins <- function(features) {
  prots <- list()
  cds_features <- Filter(function(f) f$type == "CDS", features)
  if (length(cds_features) == 0) {
    return(prots)
  }

  for (f in cds_features) {
    trans <- f$qualifiers[["translation"]]
    if (is.null(trans) || !nzchar(trans)) next

    # Build a meaningful key: prefer locus_tag, then gene, then protein_id, then product
    key <- f$qualifiers[["locus_tag"]] %||%
      f$qualifiers[["gene"]] %||%
      f$qualifiers[["protein_id"]] %||%
      f$qualifiers[["product"]] %||%
      paste0("CDS_", f$start, "_", f$end)

    key <- key[1] # take first value if multi-valued
    # Make key unique if already present
    if (key %in% names(prots)) {
      key <- paste0(key, "_", f$start)
    }
    prots[[key]] <- trans
  }
  prots
}

# ---- Utility: null-coalescing operator --------------------------------------
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b


# =============================================================================
# Reference/CDS matching
# =============================================================================
# Matches a named reference sequence (e.g. a protein-alignment FASTA header) to
# the CDS feature in a parse_genbank_file() record that it came from, and
# extracts that CDS's coding sequence. Used to recover per-isolate CDS when all
# you have is a protein alignment and the genome flatfiles.

#' Extract GenBank-style accession numbers from a string
#'
#' @param x Character vector (e.g. FASTA headers) to search.
#' @return Character vector of matched accessions (e.g. `"AB011819"`,
#'   `"NC_001445.1"`); empty if none found.
#' @export
extract_accessions <- function(x) {
  unlist(regmatches(x, gregexpr("[A-Z]{1,2}_?[0-9]{5,8}(?:\\.[0-9]+)?", x)))
}

#' Does a string contain a GenBank-style accession?
#'
#' @param x Character vector to test.
#' @return Logical vector, same length as `x`.
#' @export
has_genbank_id <- function(x) {
  grepl("[A-Z]{1,2}_?[0-9]{5,8}", x)
}

.normalize_name <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

.aa_similarity <- function(x, y) {
  n <- min(nchar(x), nchar(y))
  if (n == 0) {
    return(0)
  }
  sum(strsplit(substr(x, 1, n), "")[[1]] == strsplit(substr(y, 1, n), "")[[1]]) / n
}

#' Frame-safe translation of a (possibly gapped) coding sequence
#'
#' Strips alignment gaps, trims to a whole number of codons, and translates.
#'
#' @param seq A `Biostrings::DNAString` (or coercible via `as.character()`).
#' @return An `Biostrings::AAString` translation.
#' @export
safe_translate <- function(seq) {
  s <- gsub("-", "", as.character(seq))
  trim_len <- nchar(s) - (nchar(s) %% 3)
  s <- substr(s, 1, trim_len)
  Biostrings::translate(Biostrings::DNAString(s))
}

#' Match a reference sequence to its CDS feature and extract the coding sequence
#'
#' Finds the CDS in a parsed GenBank `record` (as returned by one element of
#' [parse_genbank_file()]'s result) that corresponds to a named reference
#' sequence, trying in order: (1) an exact `protein_id` match against `name`,
#' (2) a normalised `product`-qualifier match, and (3) — if `ref_protein` is
#' supplied — the CDS whose translation is most similar to `ref_protein`,
#' accepted only above `min_similarity`.
#'
#' @param record A single parsed GenBank record (one element of the list
#'   returned by [parse_genbank_file()] or [load_genbank_folder()]).
#' @param name Character scalar. Reference sequence name to match, e.g. a
#'   protein-alignment FASTA header containing a `protein_id`.
#' @param ref_protein Character scalar. Reference amino-acid sequence (gaps
#'   allowed), used only for the similarity fallback. Optional.
#' @param min_similarity Numeric in `[0, 1]`. Minimum fractional amino-acid
#'   identity (over the shorter of the two sequences) required to accept a
#'   similarity-based match. Default `0.7`.
#' @return A `Biostrings::DNAString` with the matched CDS (reverse-complemented
#'   if the feature is on the minus strand), or `NULL` if no feature matched or
#'   `record$sequence` is unavailable.
#' @export
match_cds_to_reference <- function(record, name, ref_protein = NULL, min_similarity = 0.7) {
  cds_feats <- Filter(function(f) identical(f$type, "CDS"), record$features)
  if (length(cds_feats) == 0 || is.na(record$sequence)) {
    return(NULL)
  }

  genome <- Biostrings::DNAString(toupper(record$sequence))

  extract <- function(feat) {
    seq <- Biostrings::subseq(genome, start = feat$start, end = feat$end)
    if (identical(feat$strand, -1L)) seq <- Biostrings::reverseComplement(seq)
    seq
  }

  # 1) protein_id exact match against the reference name
  if (nzchar(name)) {
    for (f in cds_feats) {
      pid <- f$qualifiers[["protein_id"]]
      if (!is.null(pid) && any(grepl(pid[1], name, fixed = TRUE))) {
        return(extract(f))
      }
    }
  }

  # 2) normalised product-name match
  target <- .normalize_name(name)
  if (nzchar(target)) {
    for (f in cds_feats) {
      prod <- f$qualifiers[["product"]]
      if (!is.null(prod)) {
        prod_norm <- .normalize_name(prod[1])
        if (nzchar(prod_norm) && (identical(prod_norm, target) || grepl(target, prod_norm, fixed = TRUE))) {
          return(extract(f))
        }
      }
    }
  }

  # 3) translated-CDS similarity fallback
  if (!is.null(ref_protein) && nzchar(ref_protein)) {
    best_sim <- -Inf
    best_feat <- NULL
    for (f in cds_feats) {
      aa <- safe_translate(extract(f))
      sim <- .aa_similarity(as.character(aa), ref_protein)
      if (sim > best_sim) {
        best_sim <- sim
        best_feat <- f
      }
    }
    if (!is.null(best_feat) && best_sim >= min_similarity) {
      return(extract(best_feat))
    }
  }

  NULL
}


# =============================================================================
# Convenience accessors
# =============================================================================

#' Get a flat data.frame of all features from a parsed record list
#'
#' @param records  Output of load_genbank_folder()
#' @param record_name  Name of a single record, or NULL to return all.
#' @return data.frame with columns: record, type, location, strand, start, end,
#'         and one column per unique qualifier key (NAs where absent).
#' @export
features_as_df <- function(records, record_name = NULL) {
  if (!is.null(record_name)) records <- records[record_name]

  rows <- list()
  for (rn in names(records)) {
    rec <- records[[rn]]
    for (feat in rec$features) {
      base <- list(
        record   = rn,
        type     = feat$type,
        location = feat$location,
        strand   = feat$strand,
        start    = feat$start,
        end      = feat$end
      )
      qual_flat <- lapply(feat$qualifiers, function(v) {
        if (is.logical(v)) {
          return(as.character(v))
        }
        paste(v, collapse = "; ")
      })
      rows[[length(rows) + 1L]] <- c(base, qual_flat)
    }
  }

  if (length(rows) == 0) {
    return(data.frame())
  }

  # Combine rows with differing columns (fill missing with NA)
  all_cols <- unique(unlist(lapply(rows, names)))
  df_list <- lapply(rows, function(r) {
    missing <- setdiff(all_cols, names(r))
    r[missing] <- NA
    as.data.frame(r[all_cols], stringsAsFactors = FALSE)
  })
  do.call(rbind, df_list)
}

#' Get all protein sequences as a named character vector (FASTA-ready)
#'
#' @param records  Output of load_genbank_folder()
#' @return Named character vector of protein sequences.
#' @export
get_all_proteins <- function(records) {
  prots <- list()
  for (rn in names(records)) {
    ps <- records[[rn]]$proteins
    if (length(ps) == 0) next
    names(ps) <- paste0(rn, "|", names(ps))
    prots <- c(prots, ps)
  }
  unlist(prots)
}

#' Get all nucleotide sequences as a named character vector
#'
#' @param records  Output of load_genbank_folder()
#' @return Named character vector of nucleotide sequences.
#' @export
get_all_sequences <- function(records) {
  seqs <- sapply(records, `[[`, "sequence")
  seqs[!is.na(seqs)]
}

#' Write sequences to FASTA format
#'
#' @param seqs   Named character vector of sequences.
#' @param file   Output filepath.
#' @param width  Line wrap width (default 70).
#' @export
write_fasta <- function(seqs, file, width = 70L) {
  con <- file(file, "w")
  on.exit(close(con))
  for (nm in names(seqs)) {
    cat(">", nm, "\n", sep = "", file = con)
    s <- seqs[[nm]]
    n <- nchar(s)
    pos <- 1L
    while (pos <= n) {
      cat(substr(s, pos, min(pos + width - 1L, n)), "\n", sep = "", file = con)
      pos <- pos + width
    }
  }
  invisible(NULL)
}
