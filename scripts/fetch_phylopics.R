#!/usr/bin/env Rscript
# Fetch PhyloPic silhouettes (rOpenSci rphylopic) for the taxa in a saved
# BioNexus session and cache them as PNGs so notebook 05b can embed them in node
# hover tooltips and render offline.
#
# One silhouette PER NODE, resolved at the MOST SPECIFIC taxonomic rank PhyloPic
# has: species -> genus -> family -> order -> class. This makes each node's
# silhouette actually resemble the organism (a hoverfly vs a moth vs a beetle vs
# a tick), instead of a single generic class-level image shared by everything in
# a class.
#   data/phylopic/<sanitised node name>.png
# Network: run once.

suppressPackageStartupMessages({ library(rphylopic); library(png) })
root  <- "C:/Users/chris/R_projects/project-hail-mary"
sess  <- readRDS(file.path(root, "data", "bionexus_zucchini95.rds"))
outdir <- file.path(root, "data", "phylopic"); dir.create(outdir, showWarnings = FALSE)

# Only the nodes that survive into the displayed network (biotic layers, no
# orphans) need a silhouette — fetch those to keep the run short.
biotic <- sess$edges_df[sess$edges_df$layer != "cooccurrence", ]
keep   <- unique(c(biotic$from_name, biotic$to_name))
nd <- sess$nodes_df[sess$nodes_df$name %in% keep, ]
sanitize <- function(x) gsub("[^A-Za-z0-9]", "_", x)

first_uuid <- function(nm) {
  u <- tryCatch(suppressWarnings(rphylopic::get_uuid(name = nm, n = 1)),
                error = function(e) NA_character_)
  if (length(u) == 0 || is.na(u[1])) NA_character_ else u[1]
}

log <- character()
for (i in seq_len(nrow(nd))) {
  name <- nd$name[i]
  f <- file.path(outdir, paste0(sanitize(name), ".png"))
  if (file.exists(f)) next

  genus <- sub(" .*", "", name)   # binomial -> genus
  ranks <- c(species = name, genus = genus,
             family = nd$family[i], order = nd$order[i], class = nd$class[i])
  ranks <- ranks[!is.na(ranks) & nzchar(ranks)]
  ranks <- ranks[!duplicated(ranks)]

  uid <- NA_character_; matched <- NA_character_
  for (r in names(ranks)) {
    uid <- first_uuid(ranks[[r]])
    if (!is.na(uid)) { matched <- paste0(r, " (", ranks[[r]], ")"); break }
    Sys.sleep(0.2)
  }
  if (is.na(uid)) { log <- c(log, paste0(name, ": NONE")); next }

  img <- tryCatch(rphylopic::get_phylopic(uuid = uid, format = "raster"),
                  error = function(e) NULL)
  if (is.null(img)) { log <- c(log, paste0(name, ": fetch failed")); next }
  png::writePNG(img, target = f)
  log <- c(log, paste0(name, " -> ", matched))
  Sys.sleep(0.3)
}

writeLines(log, file.path(root, "data", "phylopic", "_resolution_log.txt"))
cat("Cached silhouettes:", length(list.files(outdir, "\\.png$")), "\n")
cat(paste(log, collapse = "\n"), "\n")
