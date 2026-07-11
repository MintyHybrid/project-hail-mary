# Host Species Detection and Codon Table Retrieval for Potyviruses
# Detects host from isolate names, classifies as monocot/dicot, retrieves codon tables

# ==============================================================================
# HOST SPECIES DATABASES
# ==============================================================================

#' Common potyvirus hosts with monocot/dicot classification
#' Based on known potyvirus host ranges
POTYVIRUS_HOSTS <- list(
  # Common monocots
  monocot = c(
    "Freesia", "Pleione", "Polygonatnum", "Polygonatum", "blue squill",
    "Fritillary", "Paris", "saffron", "Zantedeschia", "Dasheen (taro)",
    "lily", "leek", "wild onion", "Narcissus", "scallion",
    "Japanese yam", "hyacinth", "konjac", "asparagus", "yam",
    "Miscanthus", "sorghum", "maize", "sugarcane", "johnsongrass",
    "Pennisetum", "Habenaria", "iris", "Gloriosa", "Clivia",
    "Dendrobium", "vanilla", "banana", "orchid", "Hippeastrum",
    "Cyrtanthus", "shallot", "onion", "canna", "johnson grass",
    "Pleioblastus", "snowdrop", "cocksfoot", "Costus"
  ),

  # Common dicots
  dicot = c(
    "Achyranthes", "peanut", "Basella", "Mirabilis", "beet", "Begonia",
    "Anemone", "Keunjorong", "Gomphocarpus", "watermelon", "soybean", "kudzu",
    "Wisteria", "Passiflora", "passionfruit", "Telosma", "yam bean", "bean",
    "Impatiens", "Hardenbergia", "passion fruit", "cowpea", "zucchini",
    "Dasheen", "Euphorbia", "Aconitum", "marigold", "turnip", "sweet potato",
    "jasmine", "Callistephus", "Apium", "celery", "Panax", "Ashitaba",
    "carrot", "potato", "tobacco", "tamarillo", "datura", "sunflower",
    "pokeweed", "clover", "Ruta (rue)", "hogweed", "periwinkle", "Thevetia",
    "lettuce", "papaya", "wild melon", "cucurbit", "chilli", "pepper",
    "eggplant", "tomato", "arracacha", "Barbacena", "Bidens", "mashua",
    "Brugmansia", "wild potato", "Verbena", "plum", "pecan", "Scorzonera",
    "endive", "Catharanthus", "noni", "Platycodon", "Thladiantha",
    "Daphne", "lupine", "pea", "Mediterranean ruda"
  )
)

# Codon usage tables for representative monocot and dicot species
# Data simplified from CoCoPUTS database / coRdon package
# Frequencies normalized to sum to 1.0

#' Monocot codon table (based on Zea mays / Triticum aestivum)
MONOCOT_CODON_TABLE <- c(
  TTT = 0.020, TTC = 0.028, TTA = 0.008, TTG = 0.020,
  TCT = 0.018, TCC = 0.021, TCA = 0.012, TCG = 0.015,
  TAT = 0.015, TAC = 0.022, TAA = 0.001, TAG = 0.001,
  TGT = 0.010, TGC = 0.014, TGA = 0.001, TGG = 0.015,
  CTT = 0.020, CTC = 0.026, CTA = 0.007, CTG = 0.028,
  CCT = 0.018, CCC = 0.020, CCA = 0.018, CCG = 0.012,
  CAT = 0.014, CAC = 0.020, CAA = 0.022, CAG = 0.032,
  CGT = 0.012, CGC = 0.014, CGA = 0.008, CGG = 0.010,
  ATT = 0.024, ATC = 0.028, ATA = 0.010, ATG = 0.028,
  ACT = 0.018, ACC = 0.024, ACA = 0.015, ACG = 0.015,
  AAT = 0.020, AAC = 0.028, AAA = 0.028, AAG = 0.042,
  AGT = 0.014, AGC = 0.018, AGA = 0.018, AGG = 0.020,
  GTT = 0.020, GTC = 0.022, GTA = 0.010, GTG = 0.030,
  GCT = 0.024, GCC = 0.026, GCA = 0.018, GCG = 0.014,
  GAT = 0.028, GAC = 0.032, GAA = 0.034, GAG = 0.040,
  GGT = 0.022, GGC = 0.026, GGA = 0.024, GGG = 0.020
)

#' Dicot codon table (based on Arabidopsis thaliana / Nicotiana tabacum)
DICOT_CODON_TABLE <- c(
  TTT = 0.022, TTC = 0.026, TTA = 0.010, TTG = 0.022,
  TCT = 0.020, TCC = 0.019, TCA = 0.015, TCG = 0.011,
  TAT = 0.017, TAC = 0.020, TAA = 0.001, TAG = 0.001,
  TGT = 0.012, TGC = 0.013, TGA = 0.001, TGG = 0.016,
  CTT = 0.024, CTC = 0.024, CTA = 0.009, CTG = 0.027,
  CCT = 0.020, CCC = 0.019, CCA = 0.021, CCG = 0.010,
  CAT = 0.016, CAC = 0.019, CAA = 0.027, CAG = 0.030,
  CGT = 0.011, CGC = 0.011, CGA = 0.007, CGG = 0.007,
  ATT = 0.027, ATC = 0.026, ATA = 0.013, ATG = 0.028,
  ACT = 0.020, ACC = 0.023, ACA = 0.018, ACG = 0.013,
  AAT = 0.024, AAC = 0.027, AAA = 0.033, AAG = 0.038,
  AGT = 0.016, AGC = 0.017, AGA = 0.024, AGG = 0.019,
  GTT = 0.024, GTC = 0.021, GTA = 0.012, GTG = 0.028,
  GCT = 0.026, GCC = 0.024, GCA = 0.021, GCG = 0.012,
  GAT = 0.033, GAC = 0.030, GAA = 0.039, GAG = 0.037,
  GGT = 0.025, GGC = 0.024, GGA = 0.027, GGG = 0.019
)

# Normalize to ensure sum = 1.0
MONOCOT_CODON_TABLE <- MONOCOT_CODON_TABLE / sum(MONOCOT_CODON_TABLE)
DICOT_CODON_TABLE <- DICOT_CODON_TABLE / sum(DICOT_CODON_TABLE)

# ==============================================================================
# HOST DETECTION FUNCTIONS
# ==============================================================================

#' Extract host species from isolate name
#' Common patterns: "VirusName_HostSpecies_IsolateID" or "VirusName-Host-ID"
#'
#' @param isolate_name Character string of isolate name
#' @return Character string of detected host (or NA)
#' @export
detect_host_from_name <- function(isolate_name) {
  # Common separators
  parts <- unlist(str_split(isolate_name, "[_\\-\\.]"))

  # Check each part against known hosts
  for (part in parts) {
    part_lower <- tolower(part)

    # Check monocots
    for (host in POTYVIRUS_HOSTS$monocot) {
      if (grepl(tolower(host), part_lower)) {
        return(host)
      }
    }

    # Check dicots
    for (host in POTYVIRUS_HOSTS$dicot) {
      if (grepl(tolower(host), part_lower)) {
        return(host)
      }
    }
  }

  return(NA)
}

#' Classify host as monocot or dicot
#'
#' @param host_name Character string of host name
#' @return "monocot", "dicot", or NA
#' @export
classify_host_type <- function(host_name) {
  if (is.na(host_name)) {
    return(NA)
  }

  host_lower <- tolower(host_name)

  # Check monocots
  for (host in POTYVIRUS_HOSTS$monocot) {
    if (grepl(tolower(host), host_lower)) {
      return("monocot")
    }
  }

  # Check dicots
  for (host in POTYVIRUS_HOSTS$dicot) {
    if (grepl(tolower(host), host_lower)) {
      return("dicot")
    }
  }

  return(NA)
}

#' Create host classification table for all isolates
#'
#' @param sequence_names Character vector of sequence names
#' @param manual_hosts Optional named vector for manual host assignment
#' @return Data frame with isolate, detected_host, host_type
#' @export
create_host_classification_table <- function(sequence_names,
                                             manual_hosts = NULL) {
  results <- data.frame(
    isolate = sequence_names,
    detected_host = NA_character_,
    host_type = NA_character_,
    source = NA_character_,
    stringsAsFactors = FALSE
  )

  # Try automatic detection
  for (i in seq_along(sequence_names)) {
    detected <- detect_host_from_name(sequence_names[i])
    if (!is.na(detected)) {
      results$detected_host[i] <- detected
      results$host_type[i] <- classify_host_type(detected)
      results$source[i] <- "auto_detected"
    }
  }

  # Apply manual overrides if provided
  if (!is.null(manual_hosts)) {
    for (isolate in names(manual_hosts)) {
      idx <- which(results$isolate == isolate)
      if (length(idx) > 0) {
        results$detected_host[idx] <- manual_hosts[isolate]
        results$host_type[idx] <- classify_host_type(manual_hosts[isolate])
        results$source[idx] <- "manual"
      }
    }
  }

  # Summary statistics
  cat("\nHost Classification Summary:\n")
  cat(sprintf("  Total isolates: %d\n", nrow(results)))
  cat(sprintf(
    "  Detected hosts: %d (%.1f%%)\n",
    sum(!is.na(results$detected_host)),
    100 * sum(!is.na(results$detected_host)) / nrow(results)
  ))
  cat(sprintf("  Monocots: %d\n", sum(results$host_type == "monocot", na.rm = TRUE)))
  cat(sprintf("  Dicots: %d\n", sum(results$host_type == "dicot", na.rm = TRUE)))
  cat(sprintf("  Unknown: %d\n", sum(is.na(results$host_type))))

  return(results)
}

# ==============================================================================
# CODON TABLE RETRIEVAL
# ==============================================================================

#' Get appropriate codon table for a host type
#'
#' @param host_type "monocot", "dicot", or NA
#' @return Named vector of codon frequencies
#' @export
get_codon_table_for_host_type <- function(host_type) {
  if (is.na(host_type)) {
    # Return average of monocot and dicot
    avg_table <- (MONOCOT_CODON_TABLE + DICOT_CODON_TABLE) / 2
    return(avg_table)
  } else if (host_type == "monocot") {
    return(MONOCOT_CODON_TABLE)
  } else if (host_type == "dicot") {
    return(DICOT_CODON_TABLE)
  } else {
    stop(sprintf("Unknown host type: %s", host_type))
  }
}

#' Create host-specific codon tables for each isolate
#'
#' @param host_classification Data frame from create_host_classification_table
#' @return Named list of codon tables (one per isolate)
#' @export
create_host_specific_codon_tables <- function(host_classification) {
  codon_tables <- list()

  for (i in seq_len(nrow(host_classification))) {
    isolate <- host_classification$isolate[i]
    host_type <- host_classification$host_type[i]

    codon_tables[[isolate]] <- get_codon_table_for_host_type(host_type)
  }

  return(codon_tables)
}

# ==============================================================================
# HOST-SPECIFIC ANALYSIS FUNCTIONS
# ==============================================================================

#' Calculate CAI (Codon Adaptation Index) for each isolate with its host
#'
#' @param sequences DNAStringSet or character vector
#' @param host_codon_tables Named list of host-specific codon tables
#' @return Data frame with isolate and CAI values
#' @export
calculate_host_specific_cai <- function(sequences, host_codon_tables) {
  if (inherits(sequences, "DNAStringSet")) {
    seq_names <- names(sequences)
    sequences <- as.character(sequences)
  } else {
    seq_names <- names(sequences)
    if (is.null(seq_names)) {
      seq_names <- paste0("seq_", seq_along(sequences))
    }
  }

  # Cache CAI weight vectors per host table - the weights depend only on the
  # host table, so compute each one once rather than per codon occurrence.
  weight_cache <- list()
  avg_table <- (MONOCOT_CODON_TABLE + DICOT_CODON_TABLE) / 2

  cai <- vapply(seq_along(sequences), function(i) {
    isolate <- seq_names[i]
    host_table <- host_codon_tables[[isolate]]
    if (is.null(host_table)) {
      isolate <- ".__avg__"
      host_table <- avg_table
    }

    w <- weight_cache[[isolate]]
    if (is.null(w)) {
      w <- .cai_weights(host_table)
      weight_cache[[isolate]] <<- w
    }

    codons <- extract_codons(sequences[i])
    if (length(codons) == 0) {
      return(NA_real_)
    }

    w_seq <- w[codons]
    w_seq <- w_seq[!is.na(w_seq)] # drop stop codons / unknowns
    if (length(w_seq) == 0) {
      return(NA_real_)
    }

    # CAI is the geometric mean of the per-codon weights
    exp(mean(log(pmax(w_seq, 1e-6))))
  }, numeric(1))

  data.frame(isolate = seq_names, CAI = cai, stringsAsFactors = FALSE)
}

#' Build a CAI weight lookup vector from a host codon-frequency table.
#'
#' The relative adaptiveness w(codon) = freq(codon) / max(freq of synonymous
#' codons). Computed once per host table and reused across all sequences.
#'
#' @param host_table Named numeric vector of codon frequencies (names = codons).
#' @return Named numeric vector of weights for every sense codon.
#' @keywords internal
.cai_weights <- function(host_table) {
  codons <- names(host_table)
  aas <- vapply(codons, function(cd) {
    aa <- GENETIC_CODE[[cd]]
    if (is.null(aa)) NA_character_ else aa
  }, character(1))

  keep <- !is.na(aas) & aas != "*"
  codons <- codons[keep]
  aas <- aas[keep]
  freqs <- host_table[codons]

  # Per-amino-acid maximum frequency, then w = freq / family_max (vectorized)
  fam_max <- tapply(freqs, aas, max, na.rm = TRUE)
  w <- freqs / fam_max[aas]
  w[!is.finite(w)] <- 1.0
  stats::setNames(as.numeric(w), codons)
}

#' Compare codon usage between isolates grouped by host type
#'
#' @param sequences DNAStringSet or character vector
#' @param host_classification Data frame with host classification
#' @return List with monocot_seqs, dicot_seqs, unknown_seqs
#' @export
group_sequences_by_host_type <- function(sequences, host_classification) {
  if (inherits(sequences, "DNAStringSet")) {
    seq_names <- names(sequences)
    sequences_char <- as.character(sequences)
  } else {
    seq_names <- names(sequences)
    sequences_char <- sequences
  }

  # Match sequences to classification
  monocot_idx <- which(host_classification$host_type == "monocot")
  dicot_idx <- which(host_classification$host_type == "dicot")
  unknown_idx <- which(is.na(host_classification$host_type))

  groups <- list(
    monocot = sequences_char[monocot_idx],
    dicot = sequences_char[dicot_idx],
    unknown = sequences_char[unknown_idx]
  )

  cat("\nSequence grouping by host type:\n")
  cat(sprintf("  Monocot-infecting: %d\n", length(groups$monocot)))
  cat(sprintf("  Dicot-infecting: %d\n", length(groups$dicot)))
  cat(sprintf("  Unknown host: %d\n", length(groups$unknown)))

  return(groups)
}

#' Analyze rare codons with host-specific context
#'
#' @param sequences DNAStringSet or character vector
#' @param host_classification Data frame with host classification
#' @param host_codon_tables Named list of host-specific tables
#' @return Data frame with rare codon analysis per isolate
analyze_host_specific_rare_codons <- function(sequences,
                                              host_classification,
                                              host_codon_tables) {
  if (inherits(sequences, "DNAStringSet")) {
    seq_names <- names(sequences)
    sequences <- as.character(sequences)
  } else {
    seq_names <- names(sequences)
  }

  avg_table <- (MONOCOT_CODON_TABLE + DICOT_CODON_TABLE) / 2
  optimal_cache <- list() # per host table: named vector aa -> optimal codon

  # Accumulate per-isolate data frames in a list, bind once at the end.
  per_isolate <- lapply(seq_along(sequences), function(i) {
    isolate <- seq_names[i]
    host_type <- host_classification$host_type[host_classification$isolate == isolate]
    host_type <- if (length(host_type) > 0) host_type[1] else NA

    key <- if (is.null(host_codon_tables[[isolate]])) ".__avg__" else isolate
    host_table <- if (key == ".__avg__") avg_table else host_codon_tables[[isolate]]

    # Precompute optimal codon per amino acid once per host table
    opt <- optimal_cache[[key]]
    if (is.null(opt)) {
      opt <- .optimal_codons(host_table)
      optimal_cache[[key]] <<- opt
    }

    codons <- unique(extract_codons(sequences[i]))
    aas <- vapply(codons, function(cd) {
      aa <- GENETIC_CODE[[cd]]
      if (is.null(aa)) NA_character_ else aa
    }, character(1))

    keep <- !is.na(aas) & aas != "*"
    codons <- codons[keep]
    aas <- aas[keep]
    if (length(codons) == 0) {
      return(NULL)
    }

    freqs <- host_table[codons]
    data.frame(
      isolate = isolate,
      host_type = host_type,
      codon = codons,
      amino_acid = aas,
      host_freq = as.numeric(freqs),
      is_rare = as.numeric(freqs) < 0.01,
      optimal_codon = opt[aas],
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, per_isolate)
}

#' Optimal (most frequent) codon per amino acid for a host table.
#'
#' @param host_table Named numeric vector of codon frequencies.
#' @return Named character vector mapping amino acid -> optimal codon.
#' @keywords internal
.optimal_codons <- function(host_table) {
  codons <- names(host_table)
  aas <- vapply(codons, function(cd) {
    aa <- GENETIC_CODE[[cd]]
    if (is.null(aa)) NA_character_ else aa
  }, character(1))
  keep <- !is.na(aas) & aas != "*"
  codons <- codons[keep]
  aas <- aas[keep]
  freqs <- host_table[codons]
  # For each amino acid, pick the codon with the highest frequency
  split_codons <- split(codons, aas)
  split_freqs <- split(as.numeric(freqs), aas)
  vapply(names(split_codons), function(a) {
    split_codons[[a]][which.max(split_freqs[[a]])]
  }, character(1))
}

# ==============================================================================
# INTERACTIVE HOST ASSIGNMENT
# ==============================================================================

#' Interactive function to manually assign hosts
#'
#' @param sequence_names Character vector of sequence names
#' @param auto_classification Data frame from create_host_classification_table
#' @return Updated classification data frame
interactive_host_assignment <- function(sequence_names, auto_classification) {
  cat("\nInteractive Host Assignment\n")
  cat("=" %s+% paste(rep("=", 60), collapse = "") %s+% "\n\n")

  # Show unclassified isolates
  unclassified <- auto_classification[is.na(auto_classification$host_type), ]

  if (nrow(unclassified) == 0) {
    cat("All isolates have been classified!\n")
    return(auto_classification)
  }

  cat(sprintf("Found %d unclassified isolates:\n", nrow(unclassified)))
  print(unclassified[, c("isolate")])

  cat("\nOptions:\n")
  cat("1. Assign all unknowns to monocot\n")
  cat("2. Assign all unknowns to dicot\n")
  cat("3. Use monocot/dicot average for unknowns\n")
  cat("4. Skip (leave as unknown)\n")

  choice <- readline(prompt = "Enter choice (1-4): ")

  if (choice == "1") {
    auto_classification$host_type[is.na(auto_classification$host_type)] <- "monocot"
    auto_classification$source[is.na(auto_classification$source)] <- "manual_batch"
  } else if (choice == "2") {
    auto_classification$host_type[is.na(auto_classification$host_type)] <- "dicot"
    auto_classification$source[is.na(auto_classification$source)] <- "manual_batch"
  } else if (choice == "3") {
    # Leave as NA, will use average table
    cat("Unknown hosts will use average monocot/dicot table\n")
  }

  return(auto_classification)
}

# ==============================================================================
# EXPORT FUNCTION
# ==============================================================================

#' Save host classification table to CSV for manual editing
#'
#' @param host_classification Data frame
#' @param filepath Output file path
#' @export
save_host_classification <- function(host_classification, filepath) {
  write.csv(host_classification, filepath, row.names = FALSE)
  cat(sprintf("\nSaved host classification to: %s\n", filepath))
  cat("You can manually edit this file and reload with:\n")
  cat(sprintf("  host_classification <- read.csv('%s')\n", filepath))
}

#' Load manually edited host classification
#'
#' @param filepath Path to CSV file
#' @return Data frame with host classification
#' @export
load_host_classification <- function(filepath) {
  df <- read.csv(filepath, stringsAsFactors = FALSE)

  # Validate
  required_cols <- c("isolate", "detected_host", "host_type")
  if (!all(required_cols %in% colnames(df))) {
    stop(sprintf(
      "CSV must contain columns: %s",
      paste(required_cols, collapse = ", ")
    ))
  }

  cat("\nLoaded host classification:\n")
  cat(sprintf("  Total isolates: %d\n", nrow(df)))
  cat(sprintf("  Monocots: %d\n", sum(df$host_type == "monocot", na.rm = TRUE)))
  cat(sprintf("  Dicots: %d\n", sum(df$host_type == "dicot", na.rm = TRUE)))
  cat(sprintf("  Unknown: %d\n", sum(is.na(df$host_type))))

  return(df)
}

# String concatenation
`%s+%` <- function(x, y) paste0(x, y)
