#!/usr/bin/env Rscript
# Prepare per-isolate CI motif-flanking region sequences + genus metadata for
# gget elm. Reads the LOCAL ICTV alignment (no live URL) and VMR table.
# Outputs (into data/elm_gget/):
#   region.fasta      - ungapped AA region per isolate (for gget elm)
#   region_meta.tsv   - name, virus_abbr, accession, genus, region_aligned

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(stringr)
  library(readxl)
  library(fuzzyjoin)
})

root <- "C:/Users/chris/R_projects/project-hail-mary"
out_dir <- file.path(root, "data", "elm_gget")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Read local alignment ------------------------------------------------
aln_path <- file.path(root, "data", "OPSR.Poty.Fig3.v16_align.txt")
protaln  <- readAAStringSet(aln_path, format = "fasta")

# --- 2. Anchor the CI motif region on the consensus -------------------------
cons    <- consensusString(protaln, ambiguityMap = "X")
m       <- matchPattern("VASYN", cons, max.mismatch = 1)
anchora <- start(m)[1]
anchorb <- end(m)[1]
region_aln <- subseq(protaln, anchora - 20, anchorb + 20)   # aligned (gapped)

# --- 3. Ungapped sequences for gget (drop '-') ------------------------------
region_ungapped <- AAStringSet(gsub("-", "", as.character(region_aln)))
names(region_ungapped) <- names(region_aln)
# gget/diamond dislike empty or ultra-short seqs; keep >= 6 aa
keep <- width(region_ungapped) >= 6
region_ungapped <- region_ungapped[keep]

# --- 4. Metadata: parse names, join ICTV VMR for genus ----------------------
nms <- names(region_aln)
tb <- tibble(nms) %>%
  mutate(virus_abbr = gsub("\\_.*", "", nms),
         accession  = gsub(".*\\_", "", nms))

ictv <- readxl::read_xlsx(file.path(root, "data", "VMR_MSL41.v1.20260320.xlsx"),
                          sheet = 2) %>%
  dplyr::filter(Family == "Potyviridae") %>%
  dplyr::select(where(~ any(!is.na(.))))

# --- 4a. Primary match: isolate accession within the VMR accession field ----
res <- fuzzyjoin::fuzzy_left_join(
  tb, ictv,
  by = c("accession" = "Virus GENBANK accession"),
  match_fun = ~ stringr::str_detect(.y, regex(.x, ignore_case = TRUE))
) %>%
  # a VMR row may match >1 way; keep first genus per input name
  dplyr::group_by(nms) %>% dplyr::slice(1) %>% dplyr::ungroup()

genus_by_acc <- res$Genus[match(nms, res$nms)]

# --- 4b. Fallback match: isolate abbreviation -> VMR abbreviation column -----
# Many added isolates are not VMR species exemplars, so their accession is
# absent; match their virus abbreviation instead. The VMR column can list
# several abbreviations per row (separated by ; , or /).
abbr_lookup <- ictv %>%
  dplyr::select(Genus, abbr = `Virus name abbreviation(s)`) %>%
  dplyr::filter(!is.na(abbr)) %>%
  tidyr::separate_rows(abbr, sep = "\\s*[;,/]\\s*") %>%
  dplyr::mutate(abbr = toupper(trimws(abbr))) %>%
  dplyr::filter(nzchar(abbr)) %>%
  dplyr::distinct(abbr, .keep_all = TRUE)

genus_by_abbr <- abbr_lookup$Genus[match(toupper(tb$virus_abbr), abbr_lookup$abbr)]

# Prefer accession match, fall back to abbreviation match
genus_col <- dplyr::coalesce(genus_by_acc, genus_by_abbr)

meta <- tibble(
  name          = nms,
  virus_abbr    = tb$virus_abbr,
  accession     = tb$accession,
  genus         = genus_col,
  genus_source  = dplyr::case_when(
    !is.na(genus_by_acc)  ~ "accession",
    !is.na(genus_by_abbr) ~ "abbreviation",
    TRUE                  ~ NA_character_),
  region_aligned = as.character(region_aln)
) %>% dplyr::filter(name %in% names(region_ungapped))

# --- 5. Write outputs -------------------------------------------------------
writeXStringSet(region_ungapped, file.path(out_dir, "region.fasta"))
readr::write_tsv(meta, file.path(out_dir, "region_meta.tsv"))

cat("Wrote", length(region_ungapped), "sequences to", file.path(out_dir, "region.fasta"), "\n")
cat("Anchor VASYN at consensus", anchora, "-", anchorb,
    "; region cols", anchora - 20, "-", anchorb + 20, "\n")
cat("Genera found:", paste(sort(unique(na.omit(meta$genus))), collapse = ", "), "\n")
cat("Genus source:", paste(names(table(meta$genus_source, useNA = "ifany")),
                           table(meta$genus_source, useNA = "ifany"),
                           sep = "=", collapse = "  "), "\n")
still_na <- meta$name[is.na(meta$genus)]
if (length(still_na)) cat("Still unclassified (", length(still_na), "): ",
                          paste(still_na, collapse = ", "), "\n", sep = "")
