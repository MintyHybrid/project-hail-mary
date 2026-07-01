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

res <- fuzzyjoin::fuzzy_left_join(
  tb, ictv,
  by = c("accession" = "Virus GENBANK accession"),
  match_fun = ~ stringr::str_detect(.y, regex(.x, ignore_case = TRUE))
)

genus_col <- if ("Genus" %in% names(res)) res$Genus else NA_character_
meta <- tibble(
  name          = nms,
  virus_abbr    = tb$virus_abbr,
  accession     = tb$accession,
  genus         = genus_col,
  region_aligned = as.character(region_aln)
) %>% dplyr::filter(name %in% names(region_ungapped))

# --- 5. Write outputs -------------------------------------------------------
writeXStringSet(region_ungapped, file.path(out_dir, "region.fasta"))
readr::write_tsv(meta, file.path(out_dir, "region_meta.tsv"))

cat("Wrote", length(region_ungapped), "sequences to", file.path(out_dir, "region.fasta"), "\n")
cat("Anchor VASYN at consensus", anchora, "-", anchorb,
    "; region cols", anchora - 20, "-", anchorb + 20, "\n")
cat("Genera found:", paste(sort(unique(na.omit(meta$genus))), collapse = ", "), "\n")
