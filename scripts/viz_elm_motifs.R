#!/usr/bin/env Rscript
# Visualize gget-elm regex motif hits across isolate CI region fragments.
# Layout: ELM motifs on y, isolate fragments on x (sorted & coloured by genus).
# Outputs: data/elm_gget/elm_motif_grid.png and a tidy joined table.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
})

root    <- "C:/Users/chris/R_projects/project-hail-mary"
dat_dir <- file.path(root, "data", "elm_gget")

regex <- readr::read_tsv(file.path(dat_dir, "elm_regex.tsv"), show_col_types = FALSE)
meta  <- readr::read_tsv(file.path(dat_dir, "region_meta.tsv"), show_col_types = FALSE)

# Join genus onto each motif hit
hits <- regex %>%
  dplyr::left_join(meta %>% dplyr::select(query = name, virus_abbr, genus),
                   by = "query") %>%
  dplyr::mutate(genus = dplyr::coalesce(genus, "Unclassified"))

# Presence of each ELM motif per isolate (collapse multiple positions)
pa <- hits %>%
  dplyr::distinct(query, virus_abbr, genus, ELMIdentifier, ELMType) %>%
  dplyr::mutate(present = 1L)

# Order isolates by genus, then abbreviation; order motifs by type then id
isolate_order <- meta %>%
  dplyr::mutate(genus = dplyr::coalesce(genus, "Unclassified")) %>%
  dplyr::arrange(genus, virus_abbr) %>%
  dplyr::pull(name)

pa <- pa %>%
  dplyr::mutate(
    query        = factor(query, levels = isolate_order),
    ELMIdentifier = forcats::fct_rev(factor(ELMIdentifier)),
    genus        = factor(genus)
  )

# Genus colour strip data (one tile row beneath the grid)
genus_strip <- meta %>%
  dplyr::mutate(genus = factor(dplyr::coalesce(genus, "Unclassified")),
                query = factor(name, levels = isolate_order)) %>%
  dplyr::filter(!is.na(query))

n_iso   <- length(isolate_order)
n_motif <- length(unique(pa$ELMIdentifier))

# Main grid: motif presence, tiles coloured by genus
p_grid <- ggplot(pa, aes(x = query, y = ELMIdentifier, fill = genus)) +
  geom_tile(color = "grey90", linewidth = 0.1) +
  scale_fill_viridis_d(option = "turbo", name = "Genus") +
  labs(
    title    = "ELM linear motifs across Potyviridae CI region fragments",
    subtitle = sprintf("%d isolates × %d motif classes (gget elm, regex mode); sorted by genus",
                       n_iso, n_motif),
    x = "Isolate fragment (sorted by genus)",
    y = "ELM motif class"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid   = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(dat_dir, "elm_motif_grid.png"), p_grid,
       width = 14, height = max(6, n_motif * 0.22), dpi = 200, limitsize = FALSE)

# Genus-characteristic motifs: motifs whose hits concentrate in few genera
genus_summary <- hits %>%
  dplyr::distinct(query, genus, ELMIdentifier) %>%
  dplyr::count(ELMIdentifier, genus, name = "n_isolates") %>%
  tidyr::pivot_wider(names_from = genus, values_from = n_isolates, values_fill = 0)

readr::write_tsv(hits, file.path(dat_dir, "elm_hits_with_genus.tsv"))
readr::write_tsv(genus_summary, file.path(dat_dir, "elm_motif_by_genus.tsv"))

cat("Isolates:", n_iso, " Motif classes:", n_motif,
    " Total hits:", nrow(hits), "\n")
cat("Wrote elm_motif_grid.png, elm_hits_with_genus.tsv, elm_motif_by_genus.tsv\n")
