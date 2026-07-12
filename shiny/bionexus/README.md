# Biodiversity Network Explorer v2
## Live API Queries + Multilayer Network Analysis

---

## Quick Start

### 1. Install packages
```r
install.packages(c(
  "shiny", "shinydashboard", "shinyWidgets", "shinycssloaders",
  "visNetwork", "leaflet",
  "ape", "igraph", "bipartite",
  "ggplot2", "ggraph", "plotly",
  "dplyr", "tidyr", "purrr", "stringr",
  "httr", "jsonlite", "DT", "RColorBrewer", "viridis"
))

# Optional but recommended for richer phylogenetic trees:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
```

### 2. Run
```r
shiny::runApp("app.R")
```

---

## What's New in v2

### ① Dynamic Species Input
Paste any list of species names (one per line, comma- or semicolon-separated).
The app automatically:
- Resolves names against **GBIF taxonomy** → kingdom, class, order, family, functional group
- Queries **GLOBI** (Global Biotic Interactions) for known interactions
- Fetches **GBIF occurrence records** for the map
- Infers additional **heuristic edges** from taxonomy (pollinator→plant, herbivore→plant, co-occurrence within orders, parasite→host) for species not covered by GLOBI

### ② Multilayer Network Analysis (Tab ③)
The interaction network is decomposed into **up to 5 layers**:

| Layer | Type | Directionality |
|-------|------|----------------|
| Pollination | Insect → Plant | Directed |
| Herbivory / Host-plant | Herbivore → Plant | Directed |
| Parasitism | Parasite → Host | Directed |
| Co-occurrence | Any ↔ Any | Undirected |
| Predation | Predator → Prey | Directed |

**Metrics computed per layer:**
- Degree (in / out / total)
- Betweenness centrality (normalised)
- Closeness centrality (normalised)

**Cross-layer metrics:**
- **Multiplex Participation Coefficient (P)**: How evenly is a species distributed across layers?
  - P ≈ 1 → species active in all layers equally ("generalist" across interaction types)
  - P ≈ 0 → species active in only one layer ("specialist")
  - Formula: `P = (L/L-1) * (1 - Σ(k_il/k_i)²)`
- **Layer-Layer Edge Overlap (Jaccard)**: What fraction of node-pairs interact in both layers?
- **Species × Layer Presence Heatmap**: Binary presence of each species per layer
- **Aggregate Network Centrality**: Collapsed (all-layer) degree, betweenness, clustering coefficient

### ③ Bipartite Network Metrics (Tab ④)
For pollination and herbivory layers, the `bipartite` package computes:
- **Connectance**: fraction of possible links that are realised
- **Nestedness** (NODF): are specialist interactions a subset of generalist ones?
- **H2'**: functional complementarity of the network
- **Weighted nestedness**
- Per-species: degree, normalised degree, species strength

---

## Data Sources

| Source | What it provides | API endpoint |
|--------|-----------------|--------------|
| **GBIF** | Taxonomy, occurrence records | `api.gbif.org/v1/` |
| **GLOBI** | Known biotic interactions | `api.globalbioticinteractions.org/interaction` |
| **Heuristic rules** | Fallback edges from functional group classification | Built-in |

### GLOBI interaction type classification
| GLOBI keyword | Mapped layer |
|---------------|-------------|
| pollinates, flowersVisitedBy, visitedBy | Pollination |
| parasiteOf, parasitizedBy, cleptoparasiteOf | Parasitism |
| eats, preysOn, preyedUponBy | Predation |
| hostOf, hasHost, feeds on (plant) | Herbivory |
| interactsWith, coOccursWith | Co-occurrence |

---

## Functional Group Inference Logic
When GBIF returns taxonomy, the app assigns a functional group:

```
Plantae                  → Plant
Lepidoptera / Sphingidae → Pollinator
Lepidoptera (other)      → Herbivore (Lepidoptera)
Hymenoptera / Apidae etc → Pollinator (bee)
Hymenoptera (other)      → Herbivore / Parasite
Gobiesociformes          → Predator (fish)
Diptera                  → Pollinator / Predator (fly)
Aves                     → Pollinator / Predator (bird)
```

---

## Extending the App

**Add a new API source:**
Write a function following the pattern of `gbif_resolve()` or `globi_interactions()` —
return a `data.frame` with `from_name`, `to_name`, `layer`, `confidence`, `source` columns —
then add it to the main `observeEvent(input$run_query, {...})` block.

**Add a new network layer:**
Add an entry to `LAYERS` at the top of the script and update `classify_interaction()`.

**Custom interaction database:**
If you have your own interaction CSV, load it and `bind_rows()` it into `edges_df` before
passing to `compute_multilayer_metrics()`.

---

## Known Caveats
- GLOBI coverage is uneven: well-studied pollinators (bees, hawk-moths) have rich data; rare species may return nothing. Heuristic edges fill the gap but are lower confidence.
- GBIF name resolution uses the species suggestion API — verify accepted names in the status table.
- Bipartite metrics require ≥ 2 plants AND ≥ 2 animals with interactions in the chosen layer.
- Phylogenetic tree is topology-only (no branch lengths); family/order-level topology may not reflect the most current molecular systematics for all groups.
- Occurrence data may contain coordinate errors (geospatial cleaning is not performed).
