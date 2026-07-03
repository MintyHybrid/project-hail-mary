#!/usr/bin/env Rscript
# Fetch PhyloPic silhouettes (rOpenSci rphylopic) for the taxonomic classes
# present in a saved BioNexus session, and cache them as PNGs so notebook 05b
# can embed them in node hover tooltips and render offline.
#   data/phylopic/<Class>.png
# Network: run once.

suppressPackageStartupMessages({ library(rphylopic); library(png) })
root  <- "C:/Users/chris/R_projects/project-hail-mary"
sess  <- readRDS(file.path(root, "data", "bionexus_turnip95.rds"))
outdir <- file.path(root, "data", "phylopic"); dir.create(outdir, showWarnings = FALSE)

taxa <- sort(unique(stats::na.omit(sess$nodes_df$class)))
taxa <- taxa[nzchar(taxa)]
message("Fetching PhyloPic silhouettes for classes: ", paste(taxa, collapse = ", "))

for (tx in taxa) {
  f <- file.path(outdir, paste0(gsub("[^A-Za-z0-9]", "_", tx), ".png"))
  if (file.exists(f)) next
  uid <- tryCatch(rphylopic::get_uuid(name = tx, n = 1), error = function(e) NA)
  if (is.na(uid)) { message("  no silhouette: ", tx); next }
  img <- tryCatch(rphylopic::get_phylopic(uuid = uid, format = "raster"),
                  error = function(e) NULL)
  if (is.null(img)) { message("  fetch failed: ", tx); next }
  # img is an RGBA raster array; write directly to PNG (no magick needed)
  png::writePNG(img, target = f)
  message("  saved ", basename(f))
  Sys.sleep(0.3)
}
cat("Cached silhouettes:", length(list.files(outdir, "\\.png$")), "\n")
