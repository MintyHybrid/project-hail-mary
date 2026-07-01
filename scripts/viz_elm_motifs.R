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

# Wilson score lower bound of a binomial proportion — shrinks toward 0.5 for
# small n, so small-genus prevalence estimates get an honest lower bound.
wilson_lower <- function(k, n, z = 1.96) {
  ifelse(n == 0, 0, {
    phat <- k / n
    denom <- 1 + z^2 / n
    centre <- phat + z^2 / (2 * n)
    margin <- z * sqrt((phat * (1 - phat) + z^2 / (4 * n)) / n)
    pmax(0, (centre - margin) / denom)
  })
}

# Small-genus confidence cutoffs (isolates in the genus driving the call)
CONF_LOW <- 5    # < 5 isolates -> low confidence
CONF_MOD <- 10   # 5-9 -> moderate; >= 10 -> high

# Completed count grid (n_pos + genus size) for qualifying genera
counts <- prevalence %>%
  dplyr::filter(genus %in% big_genera) %>%
  dplyr::select(genus, ELMIdentifier, n_pos, n_iso, prevalence) %>%
  tidyr::complete(genus, ELMIdentifier, fill = list(n_pos = 0, prevalence = 0)) %>%
  dplyr::left_join(genus_sizes, by = "genus", suffix = c("", ".g")) %>%
  dplyr::mutate(n_iso = dplyr::coalesce(n_iso, n_iso.g),
                wilson_lwr = wilson_lower(n_pos, n_iso)) %>%
  dplyr::select(-n_iso.g)

prev_grid <- counts %>% dplyr::select(genus, ELMIdentifier, prevalence)

# Specificity: variance of prevalence across genera, plus a confidence flag
# based on the size of the genus that drives each motif's signal.
spec <- counts %>%
  dplyr::group_by(ELMIdentifier) %>%
  dplyr::summarise(
    spec_var            = stats::var(prevalence),
    top_genus           = genus[which.max(prevalence)],
    top_n               = n_iso[which.max(prevalence)],
    top_prev            = max(prevalence),
    top_prev_wilson_lwr = wilson_lwr[which.max(prevalence)],
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    confidence = factor(dplyr::case_when(
      top_n >= CONF_MOD ~ "high",
      top_n >= CONF_LOW ~ "moderate",
      TRUE              ~ "low"),
      levels = c("low", "moderate", "high"))
  ) %>%
  dplyr::arrange(dplyr::desc(spec_var))

top_motifs <- spec %>% dplyr::slice_head(n = TOP_M) %>% dplyr::pull(ELMIdentifier)

# Genus axis labels annotated with n; genera < CONF_LOW isolates flagged with †
genus_label <- genus_sizes %>%
  dplyr::filter(genus %in% big_genera) %>%
  dplyr::mutate(label = sprintf("%s (n=%d)%s", genus, n_iso,
                                ifelse(n_iso < CONF_LOW, " †", "")))
lab_map <- stats::setNames(genus_label$label, genus_label$genus)

prev_top <- prev_grid %>%
  dplyr::filter(ELMIdentifier %in% top_motifs) %>%
  dplyr::mutate(
    genus_lab     = factor(lab_map[genus], levels = rev(lab_map[big_genera])),
    ELMIdentifier = factor(ELMIdentifier,
                           levels = spec$ELMIdentifier[spec$ELMIdentifier %in% top_motifs])
  )

# Mark cells whose "characteristic" genus is small (low confidence)
lowconf_motifs <- spec$ELMIdentifier[spec$confidence == "low"]
prev_top <- prev_top %>%
  dplyr::mutate(lowconf = ELMIdentifier %in% lowconf_motifs & prevalence > 0.5)

n_small <- sum(genus_sizes$genus %in% big_genera & genus_sizes$n_iso < CONF_LOW)

p_spec <- ggplot(prev_top, aes(x = ELMIdentifier, y = genus_lab, fill = prevalence)) +
  geom_tile(colour = "grey92") +
  # ring low-confidence high-prevalence cells so they read as tentative
  geom_tile(data = dplyr::filter(prev_top, lowconf),
            colour = "#2b8cbe", linewidth = 0.7, fill = NA) +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1),
                       name = "Prevalence\n(fraction of\nisolates)") +
  labs(
    title    = "Genus-characteristic ELM motifs",
    subtitle = sprintf("Top %d motifs by prevalence variance across genera (n >= %d isolates)",
                       TOP_M, MIN_N),
    caption  = paste0("† genus with <", CONF_LOW,
                      " isolates — prevalence is low-confidence; ",
                      "blue-ringed cells: motif called 'characteristic' from a small genus."),
    x = "ELM motif class", y = "Genus"
  ) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.caption = element_text(hjust = 0, size = 7, colour = "grey30"),
        panel.grid = element_blank())

ggsave(file.path(dat_dir, "elm_genus_specificity.png"), p_spec,
       width = 12, height = 5.8, dpi = 200)

# ---------------------------------------------------------------------------
# 3. Tables
# ---------------------------------------------------------------------------
genus_summary <- pa %>%
  dplyr::count(ELMIdentifier, genus, name = "n_isolates") %>%
  tidyr::pivot_wider(names_from = genus, values_from = n_isolates, values_fill = 0)

spec_table <- prev_grid %>%
  tidyr::pivot_wider(names_from = genus, values_from = prevalence, values_fill = 0) %>%
  dplyr::left_join(
    dplyr::select(spec, ELMIdentifier, top_genus, top_n, top_prev,
                  top_prev_wilson_lwr, confidence, spec_var),
    by = "ELMIdentifier") %>%
  dplyr::arrange(dplyr::desc(spec_var))

readr::write_tsv(hits, file.path(dat_dir, "elm_hits_with_genus.tsv"))
readr::write_tsv(genus_summary, file.path(dat_dir, "elm_motif_by_genus.tsv"))
readr::write_tsv(spec_table, file.path(dat_dir, "elm_motif_specificity.tsv"))

cat("Isolates:", n_iso, " Motif classes:", n_motif, " Total hits:", nrow(hits), "\n")
cat("Genera >=", MIN_N, "isolates:", paste(big_genera, collapse = ", "), "\n")
cat("Small genera (<", CONF_LOW, "isolates, low-confidence calls):",
    n_small, "\n")
cat("Top genus-characteristic motifs (with confidence):\n")
print(utils::head(spec_table[, c("ELMIdentifier", "top_genus", "top_n",
                                 "top_prev", "top_prev_wilson_lwr",
                                 "confidence", "spec_var")], 12))
