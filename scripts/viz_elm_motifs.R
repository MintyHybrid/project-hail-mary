#!/usr/bin/env Rscript
# Visualize gget-elm regex motif hits across isolate CI region fragments.
# Produces:
#   elm_motif_grid.png       - motifs (y) x isolates (x), faceted & coloured by genus
#   elm_genus_specificity.png- genus x most genus-variable motifs, fill = prevalence
#   elm_hits_with_genus.tsv  - tidy per-hit table
#   elm_motif_by_genus.tsv   - motif x genus isolate counts
#   elm_motif_specificity.tsv- per-motif genus prevalences + specificity score

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(ggplot2); library(forcats)
})

root    <- "C:/Users/chris/R_projects/project-hail-mary"
dat_dir <- file.path(root, "data", "elm_gget")
MIN_N   <- 3    # min isolates per genus to be included in specificity analysis
TOP_M   <- 30   # number of most genus-variable motifs to show

regex <- readr::read_tsv(file.path(dat_dir, "elm_regex.tsv"), show_col_types = FALSE)
meta  <- readr::read_tsv(file.path(dat_dir, "region_meta.tsv"), show_col_types = FALSE) %>%
  dplyr::mutate(genus = dplyr::coalesce(genus, "Unclassified"))

hits <- regex %>%
  dplyr::left_join(dplyr::select(meta, query = name, virus_abbr, genus), by = "query")

# Presence/absence: one row per (isolate, motif)
pa <- hits %>%
  dplyr::distinct(query, virus_abbr, genus, ELMIdentifier, ELMType)

# Genus sizes, ordered largest first for a stable layout
genus_sizes <- meta %>% dplyr::count(genus, name = "n_iso") %>% dplyr::arrange(desc(n_iso))
isolate_order <- meta %>%
  dplyr::mutate(genus = factor(genus, levels = genus_sizes$genus)) %>%
  dplyr::arrange(genus, virus_abbr) %>% dplyr::pull(name)

# ---------------------------------------------------------------------------
# 1. Main grid: motif presence, faceted by genus (panels = separators + labels)
# ---------------------------------------------------------------------------
pa_plot <- pa %>%
  dplyr::mutate(
    query         = factor(query, levels = isolate_order),
    ELMIdentifier = forcats::fct_rev(factor(ELMIdentifier)),
    genus         = factor(genus, levels = genus_sizes$genus)
  )

n_iso   <- length(isolate_order)
n_motif <- dplyr::n_distinct(pa$ELMIdentifier)

p_grid <- ggplot(pa_plot, aes(x = query, y = ELMIdentifier, fill = genus)) +
  geom_tile() +
  facet_grid(cols = vars(genus), scales = "free_x", space = "free_x") +
  scale_fill_viridis_d(option = "turbo", guide = "none") +
  labs(
    title    = "ELM linear motifs across Potyviridae CI region fragments",
    subtitle = sprintf("%d isolates x %d motif classes (gget elm, regex mode); panels = genus",
                       n_iso, n_motif),
    x = "Isolate fragment (grouped by genus)", y = "ELM motif class"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    panel.grid      = element_blank(),
    panel.spacing.x = unit(1.5, "pt"),
    strip.text.x    = element_text(angle = 90, hjust = 0, size = 7),
    strip.background = element_rect(fill = "grey95", colour = NA)
  )

ggsave(file.path(dat_dir, "elm_motif_grid.png"), p_grid,
       width = 15, height = max(7, n_motif * 0.19), dpi = 200, limitsize = FALSE)

# ---------------------------------------------------------------------------
# 2. Genus specificity: prevalence of each motif within each genus
# ---------------------------------------------------------------------------
prevalence <- pa %>%
  dplyr::count(genus, ELMIdentifier, name = "n_pos") %>%
  dplyr::left_join(genus_sizes, by = "genus") %>%
  dplyr::mutate(prevalence = n_pos / n_iso)

# Full prevalence grid (fill absent combos with 0) for genera above size cutoff
big_genera <- genus_sizes %>% dplyr::filter(n_iso >= MIN_N) %>% dplyr::pull(genus)

prev_grid <- prevalence %>%
  dplyr::filter(genus %in% big_genera) %>%
  dplyr::select(genus, ELMIdentifier, prevalence) %>%
  tidyr::complete(genus, ELMIdentifier, fill = list(prevalence = 0))

# Specificity score: variance of prevalence across the qualifying genera
spec <- prev_grid %>%
  dplyr::group_by(ELMIdentifier) %>%
  dplyr::summarise(
    max_prev   = max(prevalence),
    spec_var   = stats::var(prevalence),
    top_genus  = big_genera[which.max(tapply(prevalence, genus, mean)[big_genera])],
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(spec_var))

top_motifs <- spec %>% dplyr::slice_head(n = TOP_M) %>% dplyr::pull(ELMIdentifier)

prev_top <- prev_grid %>%
  dplyr::filter(ELMIdentifier %in% top_motifs) %>%
  dplyr::mutate(
    genus         = factor(genus, levels = rev(big_genera)),
    ELMIdentifier = factor(ELMIdentifier, levels = spec$ELMIdentifier[spec$ELMIdentifier %in% top_motifs])
  )

p_spec <- ggplot(prev_top, aes(x = ELMIdentifier, y = genus, fill = prevalence)) +
  geom_tile(colour = "grey92") +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1),
                       name = "Prevalence\n(fraction of\nisolates)") +
  labs(
    title    = "Genus-characteristic ELM motifs",
    subtitle = sprintf("Top %d motifs by prevalence variance across genera (n >= %d isolates)",
                       TOP_M, MIN_N),
    x = "ELM motif class", y = "Genus"
  ) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        panel.grid = element_blank())

ggsave(file.path(dat_dir, "elm_genus_specificity.png"), p_spec,
       width = 12, height = 5.5, dpi = 200)

# ---------------------------------------------------------------------------
# 3. Tables
# ---------------------------------------------------------------------------
genus_summary <- pa %>%
  dplyr::count(ELMIdentifier, genus, name = "n_isolates") %>%
  tidyr::pivot_wider(names_from = genus, values_from = n_isolates, values_fill = 0)

spec_table <- prev_grid %>%
  tidyr::pivot_wider(names_from = genus, values_from = prevalence, values_fill = 0) %>%
  dplyr::left_join(dplyr::select(spec, ELMIdentifier, spec_var, top_genus, max_prev),
                   by = "ELMIdentifier") %>%
  dplyr::arrange(dplyr::desc(spec_var))

readr::write_tsv(hits, file.path(dat_dir, "elm_hits_with_genus.tsv"))
readr::write_tsv(genus_summary, file.path(dat_dir, "elm_motif_by_genus.tsv"))
readr::write_tsv(spec_table, file.path(dat_dir, "elm_motif_specificity.tsv"))

cat("Isolates:", n_iso, " Motif classes:", n_motif, " Total hits:", nrow(hits), "\n")
cat("Genera >=", MIN_N, "isolates:", paste(big_genera, collapse = ", "), "\n")
cat("Top genus-characteristic motifs:\n")
print(utils::head(spec_table[, c("ELMIdentifier", "top_genus", "max_prev", "spec_var")], 10))
