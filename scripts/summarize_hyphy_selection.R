# Distills the raw HyPhy JSON/CSV exports in data/hyphy_selection/ (see the
# README there for provenance) into small derived tables the book can read
# quickly at render time, instead of loading multi-megabyte JSON per build.
#
# Run once after adding/changing anything in data/hyphy_selection/raw
# inputs; outputs are committed alongside the raw files.

library(jsonlite)
library(dplyr)

dir <- here::here("data", "hyphy_selection")

## ---- motif (18nt / 6 codon) — FEL + FUBAR side by side --------------------

fel  <- fromJSON(file.path(dir, "motif_18nt_FEL_result.json"))
fub  <- fromJSON(file.path(dir, "motif_18nt_FUBAR_result.json"))

fel_tbl <- as.data.frame(fel$MLE$content[[1]])[, 1:5]
names(fel_tbl) <- c("alpha_FEL", "beta_FEL", "alpha_eq_beta_FEL", "LRT_FEL", "p_FEL")

fub_tbl <- as.data.frame(fub$MLE$content[[1]])[, 1:6]
names(fub_tbl) <- c("alpha_FUBAR", "beta_FUBAR", "beta_minus_alpha_FUBAR",
                     "prob_purifying_FUBAR", "prob_positive_FUBAR", "bayes_factor_FUBAR")

motif_summary <- bind_cols(codon = seq_len(nrow(fel_tbl)), fel_tbl, fub_tbl) %>%
  mutate(dNdS_FEL = beta_FEL / alpha_FEL, .after = alpha_eq_beta_FEL)

write.csv(motif_summary, file.path(dir, "motif_18nt_sites_summary.csv"), row.names = FALSE)

## ---- large CI window (724 codons) — FUBAR ---------------------------------
## already exported as a tidy CSV by Datamonkey; just clean column names.

fubar_724 <- read.csv(file.path(dir, "CI_724codon_FUBAR_sites.csv"), check.names = FALSE)
names(fubar_724) <- c("site", "partition", "alpha", "beta", "beta_minus_alpha",
                       "prob_purifying", "prob_positive", "bayes_factor_positive")
fubar_724$dNdS <- fubar_724$beta / fubar_724$alpha
write.csv(fubar_724, file.path(dir, "CI_724codon_FUBAR_sites_clean.csv"), row.names = FALSE)

cat(sprintf(
  "FUBAR (724 codons): %d/%d sites purifying (P>0.9), %d/%d positive (P>0.9), median dN/dS = %.4f\n",
  sum(fubar_724$prob_purifying > 0.9), nrow(fubar_724),
  sum(fubar_724$prob_positive > 0.9), nrow(fubar_724),
  median(fubar_724$dNdS, na.rm = TRUE)))

## ---- large CI window (543 codons) — MEME ----------------------------------

meme <- fromJSON(file.path(dir, "CI_543codon_MEME_result.json"))
meme_hdr <- meme$MLE$headers[, 1]
meme_tbl <- as.data.frame(meme$MLE$content[[1]])
names(meme_tbl) <- c("alpha", "beta_neg", "p_neg", "beta_pos", "p_pos",
                      "LRT", "p_value", "n_branches_selected",
                      "total_branch_length", "MEME_logL", "FEL_logL",
                      "FEL_alpha", "FEL_beta")[seq_len(ncol(meme_tbl))]
meme_tbl$site <- seq_len(nrow(meme_tbl))
write.csv(meme_tbl, file.path(dir, "CI_543codon_MEME_sites_clean.csv"), row.names = FALSE)

cat(sprintf("MEME (543 codons): %d/%d sites with significant episodic diversifying selection (p<0.05)\n",
            sum(meme_tbl$p_value < 0.05, na.rm = TRUE), nrow(meme_tbl)))

## ---- locate the motif within the 724-codon window, and taxonomic scope ----
## The canonical 15nt catalytic core (GGCTAGCTACAACGA) is matched by exact
## string search against the labeled alignment to find its alignment-column
## codon site, independent of any assumption about how the two alignments
## (motif-only vs. 724-codon) were separately extracted.

nex_lines <- readLines(file.path(dir, "CI_724codon_alignment_labeled.nex"), warn = FALSE)
long_lines <- nex_lines[nchar(nex_lines) > 500]
taxlabels_block <- paste(long_lines[1], collapse = " ")
taxa <- gsub("'", "", regmatches(taxlabels_block, gregexpr("'[^']+'", taxlabels_block))[[1]])
seq_rows <- trimws(long_lines[-1])
stopifnot(length(taxa) == length(seq_rows))

core <- "GGCTAGCTACAACGA"
core_sites <- integer(0)
for (i in seq_along(seq_rows)) {
  p <- regexpr(core, seq_rows[i], fixed = TRUE)
  if (p[1] > 0) core_sites <- c(core_sites, (p[1] - 1) %/% 3 + 1)
}
motif_site <- as.integer(names(sort(table(core_sites), decreasing = TRUE))[1])
cat(sprintf("\nCanonical catalytic core found in %d/%d taxa, all at codon site %d\n",
            length(core_sites), length(taxa), motif_site))

motif_window <- (motif_site - 3):(motif_site + 3)
comparison <- fubar_724 %>%
  filter(site %in% motif_window) %>%
  mutate(percentile_vs_all_sites = sapply(dNdS, function(x) mean(fubar_724$dNdS <= x, na.rm = TRUE)))
write.csv(comparison, file.path(dir, "motif_window_vs_background_comparison.csv"), row.names = FALSE)
print(comparison[, c("site", "dNdS", "percentile_vs_all_sites", "prob_purifying")])

## Taxonomic scope: cross-reference against the project's own genus table.
meta <- readr::read_tsv(here::here("data", "elm_gget", "region_meta.tsv"), show_col_types = FALSE)
accessions <- sub(".*_", "", taxa)
matched <- meta[sapply(meta$accession, function(a) any(grepl(a, taxa, fixed = TRUE))), ]
cat(sprintf("\nTaxonomic scope: %d/%d taxa matched by accession against region_meta.tsv\n",
            nrow(matched), length(taxa)))
genus_tbl <- as.data.frame(table(matched$genus), stringsAsFactors = FALSE)
names(genus_tbl) <- c("genus", "n_taxa")
write.csv(genus_tbl, file.path(dir, "taxonomic_scope_by_genus.csv"), row.names = FALSE)
print(genus_tbl)
