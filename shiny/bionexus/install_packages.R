# install_packages.R — Biodiversity Network Explorer v3
# Run this once before launching app.R

pkgs <- c(
  "shiny", "shinydashboard", "shinyWidgets", "shinycssloaders",
  "colourpicker",   # NEW in v3: interactive colour picker widget
  "visNetwork", "leaflet",
  "ape", "bipartite", "igraph",
  "ggplot2", "ggraph", "plotly",
  "dplyr", "tidyr", "purrr", "stringr",
  "httr", "jsonlite", "DT", "RColorBrewer", "viridis"
)

missing <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(missing) > 0) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All CRAN packages already installed.")
}

# BioConductor (optional, for richer phylogenetic plots)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!"ggtree" %in% installed.packages()[,"Package"])
  BiocManager::install("ggtree", update = FALSE)

message("\n✅ Ready. Launch with: shiny::runApp('app.R')")
