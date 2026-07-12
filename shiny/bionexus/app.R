# ============================================================
# Biodiversity Network Explorer v3.0
#
# NEW in v3:
#  - Multiple user-defined species groups with custom names
#    and colours; group membership drives node colouring,
#    phylogeny tips, and centrality plots
#  - Species appearing in multiple groups flagged as "Multiple"
#  - Colour-by toggle: Functional Group vs User Group
#  - Session save/load: entire analysis state (taxonomy,
#    edges, occurrences, metrics, tree) stored as a single
#    .rds file; reloading skips all API calls
#  - Session summary shown on load (timestamp, counts)
#
# REQUIRED PACKAGES — install once via install_packages.R
# ============================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinydashboardPlus)
  library(shinyWidgets)
  library(shinycssloaders)
  library(colourpicker)       # NEW: colour picker widget
  library(visNetwork)
  library(leaflet)
  library(ape)
  library(bipartite)          # load BEFORE igraph (namespace fix)
  library(igraph)
  library(ggplot2)
  library(ggraph)
  library(plotly)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(httr)
  library(jsonlite)
  library(DT)
  library(RColorBrewer)
  library(viridis)
})

APP_VERSION <- "3.0"

# ═══════════════════════════════════════════════════════════════
# UTILITY HELPERS
# ═══════════════════════════════════════════════════════════════

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0) return(b)
  if (length(a) == 1 && is.na(a))   return(b)
  a
}

coalesce_vec <- function(x, replacement) {
  ifelse(is.na(x) | is.null(x), replacement, x)
}

scalar_str <- function(x, fallback = NA_character_) {
  if (is.null(x) || length(x) == 0) return(fallback)
  v <- x[[1]]
  if (is.null(v) || length(v) == 0 || is.na(v[1])) return(fallback)
  as.character(v[1])
}

safe_color <- function(lookup, key, fallback = "#607d8b") {
  key <- as.character(key[1])
  if (is.null(key) || is.na(key) || !key %in% names(lookup)) return(fallback)
  lookup[[key]]
}

empty_plotly <- function(msg = "No data available") {
  plot_ly(type = "scatter", mode = "markers") %>%
    layout(
      paper_bgcolor = "#0d1117", plot_bgcolor = "#0d1117",
      xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
      annotations = list(list(
        text = msg, x = 0.5, y = 0.5,
        xref = "paper", yref = "paper",
        showarrow = FALSE, font = list(color = "#8b949e", size = 13)
      ))
    )
}

# Default colour palette for new groups
GROUP_PALETTE <- c("#f5a623","#3498db","#2ecc71","#e74c3c","#9b59b6",
                   "#1abc9c","#e67e22","#34495e","#f39c12","#16a085",
                   "#8e44ad","#27ae60","#d35400","#2980b9","#c0392b")

# ── LAYER DEFINITIONS ─────────────────────────────────────────
LAYERS <- list(
  pollination  = list(label="Pollination",      color="#f5a623", icon="🌼"),
  herbivory    = list(label="Herbivory / Host",  color="#e74c3c", icon="🌿"),
  parasitism   = list(label="Parasitism",        color="#9b59b6", icon="🪲"),
  cooccurrence = list(label="Co-occurrence",     color="#3498db", icon="📍"),
  predation    = list(label="Predation",         color="#e67e22", icon="🦁")
)
LAYER_COLORS <- sapply(LAYERS, `[[`, "color")

FUNC_COLORS <- c(
  "Plant"                              = "#2ecc71",
  "Pollinator"                         = "#f39c12",
  "Pollinator (bee)"                   = "#e67e22",
  "Pollinator / Predator (bird)"       = "#1abc9c",
  "Herbivore (Lepidoptera)"            = "#9b59b6",
  "Herbivore / Pollinator (beetle)"    = "#8e44ad",
  "Herbivore / Parasite (Hymenoptera)" = "#7f8c8d",
  "Predator (fish)"                    = "#e74c3c",
  "Predator (arachnid)"                = "#c0392b",
  "Predator"                           = "#e74c3c",
  "Insect (unknown role)"              = "#bdc3c7",
  "Pollinator / Predator (fly)"        = "#27ae60",
  "Consumer (mammal)"                  = "#795548",
  "Unknown"                            = "#607d8b"
)

# ═══════════════════════════════════════════════════════════════
# API FUNCTIONS
# ═══════════════════════════════════════════════════════════════
# ── Retry wrapper: up to `n` attempts with exponential backoff ──
retry_get <- function(url, query, timeout_sec = 20, n = 3) {
  cfg <- httr::config(http_version = 2)   # force HTTP/1.1 (libcurl constant = 2)
  for (attempt in seq_len(n)) {
    resp <- tryCatch(
      httr::GET(url, query = query, cfg, httr::timeout(timeout_sec)),
      error = function(e) {
        message(sprintf("  attempt %d/%d failed: %s", attempt, n, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(resp) && httr::status_code(resp) == 200) return(resp)
    if (!is.null(resp)) {
      sc <- httr::status_code(resp)
      message(sprintf("  attempt %d/%d — HTTP %d", attempt, n, sc))
      if (sc == 404) return(resp)       # don't retry true 404s
    }
    if (attempt < n) Sys.sleep(2 ^ attempt)   # 2s, 4s backoff
  }
  NULL
}

gbif_resolve <- function(name) {
  tryCatch({
    resp <- retry_get("https://api.gbif.org/v1/species/suggest",
                      list(q = name, limit = 1))
    if (is.null(resp) || httr::status_code(resp) != 200) return(NULL)
    dat <- jsonlite::fromJSON(rawToChar(resp$content), flatten = TRUE)
    if (is.data.frame(dat)) {
      if (nrow(dat) == 0) return(NULL)
      d <- as.list(dat[1, ])
    } else if (is.list(dat) && length(dat) > 0) {
      d <- if (is.list(dat[[1]])) dat[[1]] else as.list(dat)
    } else return(NULL)
    list(accepted_name = scalar_str(d$canonicalName, name),
         kingdom  = scalar_str(d$kingdom,  NA_character_),
         phylum   = scalar_str(d$phylum,   NA_character_),
         class    = scalar_str(d$class,    NA_character_),
         order    = scalar_str(d$order,    NA_character_),
         family   = scalar_str(d$family,   NA_character_),
         genus    = scalar_str(d$genus,    NA_character_),
         gbif_key = scalar_str(d$key,      NA_character_))
  }, error = function(e) { message("gbif_resolve '", name, "': ", conditionMessage(e)); NULL })
}

gbif_occurrences <- function(name, limit = 100) {
  tryCatch({
    resp <- retry_get("https://api.gbif.org/v1/occurrence/search",
                      list(scientificName = name, hasCoordinate = "true", limit = limit),
                      timeout_sec = 25)
    if (is.null(resp) || httr::status_code(resp) != 200) return(NULL)
    dat <- jsonlite::fromJSON(rawToChar(resp$content), flatten = TRUE)$results
    if (is.null(dat) || !is.data.frame(dat) || nrow(dat) == 0) return(NULL)
    cols <- intersect(c("decimalLatitude","decimalLongitude","country","year","stateProvince"),
                      names(dat))
    dat <- dat[, cols, drop = FALSE]
    dat$species <- name
    dat <- dat[!is.na(dat$decimalLatitude) & !is.na(dat$decimalLongitude), ]
    if (nrow(dat) == 0) return(NULL)
    dat
  }, error = function(e) { message("gbif_occ '", name, "': ", conditionMessage(e)); NULL })
}

globi_interactions <- function(name, limit = 100) {
  tryCatch({
    resp <- retry_get("https://api.globalbioticinteractions.org/interaction",
                      list(sourceTaxon = name, limit = limit, type = "json.v2"),
                      timeout_sec = 30)
    if (is.null(resp) || httr::status_code(resp) != 200) return(NULL)
    parsed <- jsonlite::fromJSON(rawToChar(resp$content), simplifyVector = TRUE, flatten = TRUE)
    
    if (!is.null(parsed$columns) && !is.null(parsed$data)) {
      cols <- parsed$columns; rows <- parsed$data
      if (length(rows) == 0) return(NULL)
      mat <- do.call(rbind, lapply(rows, function(r) {
        v <- unlist(r); length(v) <- length(cols); v
      }))
      df <- as.data.frame(mat, stringsAsFactors = FALSE)
      names(df) <- cols
    } else if (is.data.frame(parsed)) {
      df <- parsed
    } else return(NULL)
    
    rename_map <- c(
      "sourceTaxonName"     = "source_taxon_name",
      "sourcetaxonname"     = "source_taxon_name",
      "targetTaxonName"     = "target_taxon_name",
      "targettaxonname"     = "target_taxon_name",
      "interactionTypeName" = "interaction_type",
      "interactiontypename" = "interaction_type"
    )
    for (old in names(rename_map)) {
      new <- rename_map[[old]]
      if (old %in% names(df) && !new %in% names(df))
        names(df)[names(df) == old] <- new
    }
    needed <- c("source_taxon_name","target_taxon_name","interaction_type")
    if (!all(needed %in% names(df))) return(NULL)
    df <- df[, needed, drop = FALSE]
    df <- df[complete.cases(df), ]
    df <- df[nchar(trimws(df$source_taxon_name)) > 0 &
               nchar(trimws(df$target_taxon_name)) > 0, ]
    if (nrow(df) == 0) return(NULL)
    df$query_species <- name
    df
  }, error = function(e) { message("globi '", name, "': ", conditionMessage(e)); NULL })
}


classify_interaction <- function(int_type) {
  s <- tolower(as.character(int_type[1]))
  if (grepl("pollinat|flower|nectar|visit",   s)) return("pollination")
  if (grepl("parasit|cleptoparasit",          s)) return("parasitism")
  if (grepl("host|herbiv|feed|consume|gall|mine|bore", s)) return("herbivory")
  if (grepl("prey|predat|eat|hunt|kill",      s)) return("predation")
  "cooccurrence"
}

infer_functional_group <- function(kingdom, class, order, family) {
  kg <- as.character(kingdom[1]); cl <- as.character(class[1])
  or <- as.character(order[1]);   fa <- as.character(family[1])
  if (!is.na(kg) && kg=="Plantae") return("Plant")
  if (!is.na(or)) {
    if (or=="Lepidoptera") {
      if (!is.na(fa) && fa %in% c("Sphingidae","Hesperiidae","Pieridae")) return("Pollinator")
      return("Herbivore (Lepidoptera)")
    }
    if (or=="Hymenoptera") {
      if (!is.na(fa) && fa %in% c("Apidae","Halictidae","Megachilidae","Colletidae","Andrenidae"))
        return("Pollinator (bee)")
      if (!is.na(fa) && fa=="Vespidae") return("Predator")
      return("Herbivore / Parasite (Hymenoptera)")
    }
    if (or=="Diptera")    return("Pollinator / Predator (fly)")
    if (or=="Coleoptera") return("Herbivore / Pollinator (beetle)")
    if (or %in% c("Araneae","Scorpiones"))                          return("Predator (arachnid)")
    if (or %in% c("Gobiesociformes","Perciformes","Cypriniformes")) return("Predator (fish)")
  }
  if (!is.na(cl)) {
    if (cl=="Insecta")  return("Insect (unknown role)")
    if (cl=="Aves")     return("Pollinator / Predator (bird)")
    if (cl=="Mammalia") return("Consumer (mammal)")
  }
  "Unknown"
}

build_newick_from_taxonomy <- function(species_df) {
  df <- species_df[!is.na(species_df$order) & nchar(species_df$order)>0, ]
  if (nrow(df)<2) return(NULL)
  tip_names <- make.unique(gsub("[ /()]+","_", df$name))
  order_groups <- split(tip_names, df$order)
  order_subtrees <- sapply(names(order_groups), function(ord) {
    tips <- order_groups[[ord]]
    ord_safe <- gsub("[^A-Za-z0-9_]","_",ord)
    if (length(tips)==1) tips[1]
    else paste0("(",paste(tips,collapse=","),")",ord_safe)
  }, USE.NAMES=TRUE)
  order_kingdom <- df$kingdom[match(names(order_groups), df$order)]
  kingdom_groups <- split(names(order_groups), order_kingdom)
  kingdom_subtrees <- sapply(names(kingdom_groups), function(kng) {
    subs <- order_subtrees[kingdom_groups[[kng]]]
    kng_safe <- gsub("[^A-Za-z0-9_]","_",kng)
    if (length(subs)==1) subs[[1]]
    else paste0("(",paste(subs,collapse=","),")",kng_safe)
  })
  newick <- paste0("(",paste(kingdom_subtrees,collapse=","),")Life;")
  tryCatch(read.tree(text=newick), error=function(e) NULL)
}

# ═══════════════════════════════════════════════════════════════
# MULTILAYER NETWORK
# ═══════════════════════════════════════════════════════════════

build_layer_graphs <- function(edges_df, nodes_df) {
  layers <- unique(edges_df$layer)
  directed_layers <- c("pollination","herbivory","parasitism","predation")
  graphs <- lapply(setNames(layers,layers), function(lyr) {
    e <- edges_df[edges_df$layer==lyr, ]
    if (nrow(e)==0) return(NULL)
    is_dir <- lyr %in% directed_layers
    g <- tryCatch(
      graph_from_data_frame(
        d        = e[, c("from_name","to_name"), drop=FALSE],
        vertices = data.frame(name=nodes_df$name, stringsAsFactors=FALSE),
        directed = is_dir
      ),
      error=function(e){ message("build_layer_graphs '",lyr,"': ",conditionMessage(e)); NULL }
    )
    if (!is.null(g)) g$layer_directed <- is_dir
    g
  })
  Filter(Negate(is.null), graphs)
}

multiplex_participation <- function(layer_graphs, all_names) {
  L <- length(layer_graphs)
  if (L<2) return(setNames(rep(NA_real_,length(all_names)), all_names))
  deg_matrix <- matrix(0.0, nrow=length(all_names), ncol=L,
                       dimnames=list(all_names, names(layer_graphs)))
  for (lyr in names(layer_graphs)) {
    g <- layer_graphs[[lyr]]
    deg <- degree(g, mode="all")
    common <- intersect(names(deg), all_names)
    deg_matrix[common, lyr] <- as.numeric(deg[common])
  }
  ki <- rowSums(deg_matrix)
  ifelse(ki==0, 0, (L/(L-1))*(1-rowSums((deg_matrix/ifelse(ki==0,1,ki))^2)))
}

compute_multilayer_metrics <- function(edges_df, nodes_df) {
  layer_graphs <- build_layer_graphs(edges_df, nodes_df)
  if (length(layer_graphs)==0) { message("No layer graphs built"); return(NULL) }
  all_names <- nodes_df$name
  results   <- list()
  message(sprintf("Multilayer: %d layers, %d species", length(layer_graphs), length(all_names)))

  layer_metrics <- lapply(names(layer_graphs), function(lyr) {
    g  <- layer_graphs[[lyr]]; vn <- V(g)$name
    is_dir <- isTRUE(g$layer_directed)
    in_mode <- if (is_dir) "in" else "all"
    out_mode<- if (is_dir) "out" else "all"
    tryCatch(data.frame(
      species=vn, layer=lyr,
      degree     =as.numeric(degree(g,mode="all")[vn]),
      in_degree  =as.numeric(degree(g,mode=in_mode)[vn]),
      out_degree =as.numeric(degree(g,mode=out_mode)[vn]),
      betweenness=as.numeric(betweenness(g,normalized=TRUE)[vn]),
      closeness  =as.numeric(closeness(g,normalized=TRUE)[vn]),
      stringsAsFactors=FALSE),
    error=function(e){ message("Layer metrics '",lyr,"': ",conditionMessage(e)); NULL })
  })
  results$layer_metrics <- bind_rows(Filter(Negate(is.null), layer_metrics))

  agg_edges <- edges_df %>% select(from_name,to_name) %>% distinct() %>%
    filter(from_name %in% all_names, to_name %in% all_names)
  if (nrow(agg_edges)>0) {
    g_agg <- tryCatch(
      graph_from_data_frame(agg_edges,
        vertices=data.frame(name=all_names,stringsAsFactors=FALSE), directed=FALSE),
      error=function(e) NULL)
    if (!is.null(g_agg)) {
      vn2 <- V(g_agg)$name
      results$aggregate_degree      <- setNames(as.numeric(degree(g_agg)[vn2]), vn2)
      results$aggregate_betweenness <- setNames(as.numeric(betweenness(g_agg,normalized=TRUE)[vn2]), vn2)
      results$aggregate_clustering  <- setNames(
        as.numeric(transitivity(g_agg,type="local",isolates="zero"))[match(all_names,vn2)], all_names)
    }
  }

  results$participation <- multiplex_participation(layer_graphs, all_names)

  if (length(layer_graphs)>=2) {
    lnames <- names(layer_graphs)
    overlap_mat <- matrix(0.0,length(lnames),length(lnames),dimnames=list(lnames,lnames))
    make_edge_set <- function(g) {
      el <- as_edgelist(g,names=TRUE)
      if (nrow(el)==0) return(character(0))
      apply(el,1,function(r) paste(sort(r),collapse="~"))
    }
    edge_sets <- lapply(layer_graphs, make_edge_set)
    for (i in seq_along(lnames)) for (j in seq_along(lnames)) {
      if (i==j) { overlap_mat[i,j]<-1; next }
      si <- edge_sets[[lnames[i]]]; sj <- edge_sets[[lnames[j]]]
      if (!length(si)||!length(sj)) next
      u <- length(union(si,sj))
      overlap_mat[i,j] <- if (u==0) 0 else length(intersect(si,sj))/u
    }
    results$layer_overlap <- overlap_mat
  }

  layer_presence <- matrix(FALSE, nrow=length(all_names), ncol=length(layer_graphs),
                           dimnames=list(all_names, names(layer_graphs)))
  for (lyr in names(layer_graphs)) {
    g <- layer_graphs[[lyr]]
    active <- V(g)$name[degree(g,mode="all")>0]
    valid  <- intersect(active, all_names)
    if (length(valid)>0) layer_presence[valid,lyr] <- TRUE
  }
  results$layer_presence <- layer_presence

  agg_deg <- if (!is.null(results$aggregate_degree))
    results$aggregate_degree[all_names] else rep(0.0,length(all_names))
  agg_btw <- if (!is.null(results$aggregate_betweenness))
    results$aggregate_betweenness[all_names] else rep(0.0,length(all_names))

  results$species_summary <- data.frame(
    species             = all_names,
    n_layers_active     = rowSums(layer_presence),
    participation_coeff = as.numeric(results$participation[all_names]),
    agg_degree          = as.numeric(coalesce_vec(agg_deg,0)),
    agg_betweenness     = as.numeric(coalesce_vec(agg_btw,0)),
    stringsAsFactors=FALSE
  ) %>% left_join(nodes_df %>% select(name,functional_group,user_group,family,order),
                  by=c("species"="name"))

  results$layer_graphs <- layer_graphs
  results
}

compute_bipartite_metrics <- function(edges_df, nodes_df, layer_name) {
  e <- edges_df[edges_df$layer==layer_name, ]
  if (nrow(e)<2) return(NULL)
  plants  <- nodes_df$name[!is.na(nodes_df$kingdom) & nodes_df$kingdom=="Plantae"]
  animals <- nodes_df$name[!is.na(nodes_df$kingdom) & nodes_df$kingdom=="Animalia"]
  if (!length(plants)||!length(animals)) return(NULL)
  from_is_plant <- e$from_name %in% plants
  e$plant  <- ifelse(from_is_plant, e$from_name, e$to_name)
  e$animal <- ifelse(from_is_plant, e$to_name,   e$from_name)
  e <- e[e$plant %in% plants & e$animal %in% animals & e$plant!=e$animal, ]
  if (nrow(e)==0) return(NULL)
  web <- table(e$plant, e$animal)
  if (nrow(web)<2||ncol(web)<2) return(NULL)
  tryCatch({
    nm <- networklevel(web, index=c("connectance","links per species",
                                    "nestedness","H2","weighted nestedness"))
    sp <- specieslevel(web, index=c("degree","normalised degree","species strength"))
    list(network_metrics=nm, species_metrics=sp, web=web)
  }, error=function(e){ message("bipartite: ",conditionMessage(e)); NULL })
}

generate_heuristic_edges <- function(nodes_df) {
  edges <- list()
  plants      <- nodes_df$name[!is.na(nodes_df$kingdom) & nodes_df$kingdom=="Plantae"]
  pollinators <- nodes_df$name[grepl("Pollinator",nodes_df$functional_group)]
  herbivores  <- nodes_df$name[grepl("Herbivore", nodes_df$functional_group)]
  parasites   <- nodes_df$name[grepl("Parasite",  nodes_df$functional_group)]
  animals     <- nodes_df$name[!is.na(nodes_df$kingdom) & nodes_df$kingdom=="Animalia"]
  make_e <- function(from,to,layer,conf)
    data.frame(from_name=from,to_name=to,layer=layer,
               confidence=conf,source="Heuristic",stringsAsFactors=FALSE)
  for (p in pollinators) for (pl in plants)
    edges[[length(edges)+1]] <- make_e(p,pl,"pollination","Heuristic-Low")
  for (h in herbivores) for (pl in plants)
    edges[[length(edges)+1]] <- make_e(h,pl,"herbivory","Heuristic-Low")
  non_par <- animals[!animals %in% parasites]
  for (pa in parasites) for (ph in non_par)
    edges[[length(edges)+1]] <- make_e(pa,ph,"parasitism","Heuristic-Low")
  if (nrow(nodes_df)>1) {
    for (i in seq_len(nrow(nodes_df)-1)) for (j in (i+1):nrow(nodes_df)) {
      a <- nodes_df$name[i]; b <- nodes_df$name[j]
      oa <- nodes_df$order[i]; ob <- nodes_df$order[j]
      same_order <- !is.na(oa)&&!is.na(ob)&&oa==ob
      conf <- if (same_order) "Heuristic-Moderate" else "Heuristic-Low"
      edges[[length(edges)+1]] <- make_e(a,b,"cooccurrence",conf)
    }
  }
  if (!length(edges)) return(NULL)
  bind_rows(edges) %>% filter(from_name!=to_name) %>%
    distinct(from_name,to_name,layer,.keep_all=TRUE)
}

# ── Resolve groups list → flat species vector + user_group column ──
resolve_groups <- function(groups_list) {
  # groups_list: named list, each element = list(color, species=char vec)
  # Returns: data.frame with columns: name, user_group, user_group_color
  seen <- character(0)
  rows <- list()
  for (gname in names(groups_list)) {
    g   <- groups_list[[gname]]
    sps <- g$species
    col <- g$color %||% "#607d8b"
    for (sp in sps) {
      if (sp %in% seen) {
        # Mark as multiple — update existing row
        idx <- which(sapply(rows, function(r) r$name==sp))
        if (length(idx)>0) {
          rows[[idx[1]]]$user_group       <- "Multiple"
          rows[[idx[1]]]$user_group_color <- "#aaaaaa"
        }
      } else {
        seen <- c(seen, sp)
        rows[[length(rows)+1]] <- data.frame(
          name             = sp,
          user_group       = gname,
          user_group_color = col,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) return(data.frame(name=character(),
                                        user_group=character(),
                                        user_group_color=character(),
                                        stringsAsFactors=FALSE))
  bind_rows(rows)
}

# ═══════════════════════════════════════════════════════════════
# UI
# ═══════════════════════════════════════════════════════════════
ui <- dashboardPage(
  skin = "midnight",
  dashboardHeader(title=tags$span("🌿 Biodiversity Network v3"), titleWidth=310),

  dashboardSidebar(
    width=280,
    sidebarMenu(id="tabs",
      menuItem("① Species Groups",      tabName="input",     icon=icon("layer-group")),
      menuItem("② Network",             tabName="network",   icon=icon("project-diagram")),
      menuItem("③ Multilayer Analysis", tabName="multi",     icon=icon("chart-bar")),
      menuItem("④ Bipartite Metrics",   tabName="bipartite", icon=icon("sitemap")),
      menuItem("⑤ Phylogeny",           tabName="phylo",     icon=icon("tree")),
      menuItem("⑥ Map",                 tabName="map",       icon=icon("map")),
      menuItem("⑦ Species Info",        tabName="info",      icon=icon("info-circle"))
    ),
    hr(),
    # ── Save / Load (always visible) ──────────────────────────
    div(style="padding:0 14px;",
      h5("💾 Session", style="color:#8b949e; margin-bottom:6px;"),
      downloadButton("save_session", "Save session (.rds)",
                     style="width:100%;margin-bottom:6px;font-size:12px;"),
      fileInput("load_session", NULL, accept=".rds",
                placeholder="Load session (.rds)",
                buttonLabel="📂 Load",
                width="100%"),
      uiOutput("session_info_sidebar")
    ),
    hr(),
    # ── Network display controls ────────────────────────────────
    conditionalPanel("input.tabs == 'network'",
      div(style="padding:0 14px;",
        h5("Network Display", style="color:#8b949e;"),
        radioButtons("net_view", NULL,
          choices=c("All layers"="all","Single layer"="single"), selected="all"),
        conditionalPanel("input.net_view == 'single'",
          selectInput("net_layer","Layer:", choices=names(LAYERS))
        ),
        radioButtons("color_by", "Colour nodes by:",
          choices=c("Functional group"="functional_group",
                    "User group"="user_group"),
          selected="functional_group"),
        checkboxGroupInput("conf_filter","Show confidence:",
          choices=c("High","Moderate","Low","GLOBI",
                    "Heuristic-Moderate","Heuristic-Low"),
          selected=c("High","Moderate","Low","GLOBI","Heuristic-Moderate"))
      )
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background:#0d1117; }
      .box { background:#161b22; border-top:2px solid #30363d; color:#c9d1d9; }
      .box-header { color:#e6edf3; }
      .box-header .box-title { font-size:1rem; }
      body,label,.control-label,.checkbox,p,li { color:#c9d1d9 !important; }
      .skin-midnight .main-sidebar { background:#0d1117; }
      .skin-midnight .sidebar-menu>li>a { color:#8b949e; }
      .skin-midnight .sidebar-menu>li.active>a { color:#58a6ff; border-left:3px solid #58a6ff; }
      .layer-badge { display:inline-block; padding:3px 10px; border-radius:12px;
                     font-size:11px; font-weight:600; margin:2px; }
      .group-panel { background:#0d1117; border:1px solid #30363d; border-radius:8px;
                     padding:12px; margin-bottom:10px; position:relative; }
      .group-remove-btn { position:absolute; top:8px; right:8px; }
      .metric-card { background:#0d1117; border:1px solid #30363d; border-radius:8px;
                     padding:14px; margin-bottom:12px; text-align:center; }
      .metric-card .val { font-size:1.8rem; color:#58a6ff; font-weight:bold; }
      .metric-card .lbl { font-size:0.7rem; color:#8b949e; text-transform:uppercase;
                          letter-spacing:.08em; }
      .status-ok   { color:#3fb950; font-size:12px; }
      .status-wait { color:#d29922; font-size:12px; }
      .sp-status-tbl td { padding:4px 8px; font-size:12px; vertical-align:middle; }
      .session-info { background:#0d1117; border:1px solid #30363d; border-radius:6px;
                      padding:8px; font-size:11px; color:#8b949e; margin-top:6px; }
      /* colourpicker button sizing */
      .btn-colourpicker { height:30px !important; width:44px !important; border-radius:4px; }
    "))),

    tabItems(

      # ══════ TAB 1: SPECIES GROUPS ═════════════════════════════
      tabItem(tabName="input",
        fluidRow(
          box(width=8, title="Species Groups", status="primary", solidHeader=TRUE,
            p("Define one or more groups. Each group can represent a species list, study site, taxon set, etc.",
              style="color:#8b949e;font-size:13px;"),
            # Dynamic group panels inserted here
            div(id="group_container", uiOutput("group_panels_ui")),
            fluidRow(
              column(6,
                actionButton("add_group","➕ Add Group",
                             class="btn-info", style="width:100%;margin-top:4px;")),
              column(6,
                actionButton("clear_groups","🗑 Reset All",
                             class="btn-warning", style="width:100%;margin-top:4px;"))
            ),
            hr(),
            fluidRow(
              column(4,
                numericInput("gbif_occ_limit","GBIF occurrences/species:",
                             100,min=10,max=500,step=10)),
              column(4,
                numericInput("globi_limit","GLOBI interactions/species:",
                             50,min=10,max=200,step=10)),
              column(4,
                checkboxGroupInput("query_sources","Sources:",
                  choices=c("GBIF taxonomy"="gbif","GLOBI"="globi","Occurrences"="occ"),
                  selected=c("gbif","globi","occ")))
            ),
            actionButton("run_query","🔍 Query All Sources",
                         class="btn-primary btn-lg", width="100%")
          ),
          box(width=4, title="Query Status & Session", status="info", solidHeader=TRUE,
            withSpinner(uiOutput("query_status_ui"), color="#58a6ff"),
            hr(),
            h5("Resolved Taxonomy"),
            withSpinner(DTOutput("tax_preview"), color="#58a6ff")
          )
        ),
        fluidRow(
          box(width=12, title="GLOBI Raw Interactions", status="info", solidHeader=TRUE,
            p("Raw interactions from the Global Biotic Interactions database.",
              style="color:#8b949e;font-size:13px;"),
            withSpinner(DTOutput("globi_raw_table"), color="#58a6ff")
          )
        )
      ),

      # ══════ TAB 2: NETWORK ════════════════════════════════════
      tabItem(tabName="network",
        fluidRow(
          box(width=12, title="Interaction Network", status="primary", solidHeader=TRUE,
            uiOutput("layer_legend_ui"),
            uiOutput("group_legend_ui"),
            withSpinner(visNetworkOutput("main_network",height="580px"), color="#58a6ff")
          )
        ),
        fluidRow(
          box(width=5, title="Selected Entity", status="info", solidHeader=TRUE,
            uiOutput("net_detail_ui")),
          box(width=7, title="Per-Layer Degree", status="info", solidHeader=TRUE,
            withSpinner(plotlyOutput("layer_degree_chart",height="280px"), color="#58a6ff"))
        )
      ),

      # ══════ TAB 3: MULTILAYER ANALYSIS ════════════════════════
      tabItem(tabName="multi",
        fluidRow(uiOutput("multi_summary_cards")),
        fluidRow(
          box(width=6, title="Multiplex Participation Coefficient",
              status="primary", solidHeader=TRUE,
            p("P=1 → equally active across all layers; P=0 → active in only one.",
              style="color:#8b949e;font-size:12px;"),
            radioButtons("multi_color_by","Colour by:",
              choices=c("Functional group"="functional_group","User group"="user_group"),
              selected="functional_group", inline=TRUE),
            withSpinner(plotlyOutput("participation_plot",height="340px"), color="#58a6ff")
          ),
          box(width=6, title="Layer-Layer Edge Overlap (Jaccard)",
              status="primary", solidHeader=TRUE,
            withSpinner(plotlyOutput("layer_overlap_plot",height="340px"), color="#58a6ff"))
        ),
        fluidRow(
          box(width=6, title="Species × Layer Presence", status="info", solidHeader=TRUE,
            withSpinner(plotlyOutput("presence_heatmap",height="360px"), color="#58a6ff")),
          box(width=6, title="Aggregate Centrality", status="info", solidHeader=TRUE,
            p("Degree vs betweenness; bubble size = clustering coefficient.",
              style="color:#8b949e;font-size:12px;"),
            withSpinner(plotlyOutput("centrality_plot",height="360px"), color="#58a6ff"))
        ),
        fluidRow(
          box(width=12, title="Full Multilayer Metrics Table",
              status="info", solidHeader=TRUE,
            withSpinner(DTOutput("multi_table"), color="#58a6ff"))
        )
      ),

      # ══════ TAB 4: BIPARTITE ══════════════════════════════════
      tabItem(tabName="bipartite",
        fluidRow(
          box(width=3, title="Select Layer", status="primary", solidHeader=TRUE,
            selectInput("bp_layer","Layer:",
                        choices=c("pollination","herbivory","parasitism")),
            p("Requires ≥2 plants AND ≥2 animals with interactions.",
              style="color:#8b949e;font-size:12px;")),
          box(width=9, title="Network-Level Bipartite Metrics",
              status="info", solidHeader=TRUE,
            withSpinner(uiOutput("bp_net_metrics"), color="#58a6ff"))
        ),
        fluidRow(
          box(width=6, title="Interaction Web Matrix",
              status="primary", solidHeader=TRUE,
            withSpinner(plotlyOutput("bp_web_plot",height="380px"), color="#58a6ff")),
          box(width=6, title="Species-Level Metrics",
              status="info", solidHeader=TRUE,
            withSpinner(DTOutput("bp_species_table"), color="#58a6ff"))
        )
      ),

      # ══════ TAB 5: PHYLOGENY ══════════════════════════════════
      tabItem(tabName="phylo",
        fluidRow(
          box(width=12, title="Phylogenetic Tree (family/order topology)",
              status="primary", solidHeader=TRUE,
            fluidRow(
              column(6,
                p("Topology from GBIF taxonomy; branch lengths not to scale.",
                  style="color:#8b949e;font-size:13px;")),
              column(6,
                radioButtons("phylo_color_by","Colour tips by:",
                  choices=c("Functional group"="functional_group","User group"="user_group"),
                  selected="functional_group", inline=TRUE))
            ),
            withSpinner(plotlyOutput("phylo_plot",height="580px"), color="#58a6ff")
          )
        )
      ),

      # ══════ TAB 6: MAP ════════════════════════════════════════
      tabItem(tabName="map",
        fluidRow(
          box(width=12, title="Species Occurrences (GBIF)",
              status="primary", solidHeader=TRUE,
            fluidRow(
              column(9, uiOutput("map_sp_picker")),
              column(3, br(),
                actionButton("refresh_map","🗺 Refresh Map",
                             class="btn-primary",width="100%"))
            ),
            withSpinner(leafletOutput("occ_map",height="520px"), color="#58a6ff")
          )
        )
      ),

      # ══════ TAB 7: SPECIES INFO ═══════════════════════════════
      tabItem(tabName="info",
        fluidRow(
          box(width=4, title="Species Card", status="primary", solidHeader=TRUE,
            uiOutput("info_selector"), hr(), uiOutput("species_card_ui")),
          box(width=8, title="All Species Overview",
              status="info", solidHeader=TRUE,
            withSpinner(DTOutput("species_full_table"), color="#58a6ff"))
        )
      )
    )
  )
)

# ═══════════════════════════════════════════════════════════════
# SERVER
# ═══════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  # ── Reactive state ───────────────────────────────────────────
  rv <- reactiveValues(
    # Group management
    n_groups    = 1L,            # number of group panels currently shown
    group_ids   = 1L,            # running counter for unique IDs (never reused)
    active_ids  = 1L,            # IDs of currently visible panels
    # Data
    groups      = list(),        # resolved: named list(color, species)
    taxonomy    = list(),
    globi_raw   = NULL,
    occurrences = list(),
    query_log   = list(),
    edges_df    = NULL,
    nodes_df    = NULL,          # includes user_group, user_group_color
    multi_metrics= NULL,
    phylo_tree  = NULL,
    session_meta= NULL           # info about loaded session
  )

  # ── Helpers ─────────────────────────────────────────────────
  parse_species <- function(txt) {
    sp <- unlist(str_split(txt, "[,\n;]+"))
    unique(str_trim(sp)[nchar(str_trim(sp))>0])
  }

  # Build color lookup from nodes_df for a given color_by mode
  get_node_colors <- function(nodes_df, color_by="functional_group") {
    if (color_by=="user_group" && "user_group_color" %in% names(nodes_df)) {
      setNames(nodes_df$user_group_color, nodes_df$name)
    } else {
      sapply(nodes_df$functional_group, function(fg) safe_color(FUNC_COLORS, fg),
             USE.NAMES=FALSE) %>% setNames(nodes_df$name)
    }
  }

  # ═══════════════════════════════════════════════════════════
  # DYNAMIC GROUP UI
  # ═══════════════════════════════════════════════════════════

  # Build a single group panel UI
  make_group_panel <- function(i, default_name=NULL, default_color=NULL,
                                default_species=NULL) {
    gname  <- default_name    %||% paste("Group", i)
    gcol   <- default_color   %||% GROUP_PALETTE[(i-1) %% length(GROUP_PALETTE) + 1]
    gspecs <- default_species %||% ""

    div(id=paste0("group_panel_",i), class="group-panel",
      fluidRow(
        column(5,
          textInput(paste0("gname_",i), "Group name:", value=gname,
                    width="100%")
        ),
        column(3,
          tags$label("Colour:", class="control-label"),
          colourInput(paste0("gcol_",i), NULL, value=gcol,
                      showColour="background", width="100%")
        ),
        column(4,
          br(),
          actionButton(paste0("remove_group_",i), "✕ Remove",
                       class="btn-danger btn-sm group-remove-btn",
                       style="width:100%;")
        )
      ),
      textAreaInput(paste0("gspecies_",i), "Species (one per line):",
                    value=gspecs, rows=5, width="100%",
                    placeholder="Species name 1\nSpecies name 2\n...")
    )
  }

  output$group_panels_ui <- renderUI({
    ids <- rv$active_ids
    panels <- lapply(ids, function(i) make_group_panel(i))
    do.call(tagList, panels)
  })

  # Add group button
  observeEvent(input$add_group, {
    new_id <- rv$group_ids + 1L
    rv$group_ids  <- new_id
    rv$active_ids <- c(rv$active_ids, new_id)
    rv$n_groups   <- rv$n_groups + 1L
  })

  # Remove individual group — one observer per potential panel
  # Use a factory to avoid the loop-closure problem
  lapply(1:50, function(i) {
    observeEvent(input[[paste0("remove_group_",i)]], {
      if (length(rv$active_ids) <= 1) {
        showNotification("Keep at least one group.", type="warning")
        return()
      }
      rv$active_ids <- rv$active_ids[rv$active_ids != i]
      rv$n_groups   <- rv$n_groups - 1L
    }, ignoreNULL=TRUE, ignoreInit=TRUE)
  })

  # Reset all groups
  observeEvent(input$clear_groups, {
    rv$n_groups   <- 1L
    rv$group_ids  <- 1L
    rv$active_ids <- 1L
  })

  # Read current group definitions from inputs
  read_groups_from_input <- function() {
    ids <- rv$active_ids
    groups <- list()
    for (i in ids) {
      gname <- trimws(input[[paste0("gname_",i)]] %||% paste("Group",i))
      gcol  <- input[[paste0("gcol_",i)]]  %||% "#607d8b"
      gsps  <- parse_species(input[[paste0("gspecies_",i)]] %||% "")
      if (length(gsps)==0) next
      if (gname=="" || is.na(gname)) gname <- paste("Group",i)
      # Deduplicate group names
      if (gname %in% names(groups)) gname <- paste0(gname,"_",i)
      groups[[gname]] <- list(color=gcol, species=gsps)
    }
    groups
  }

  # ═══════════════════════════════════════════════════════════
  # SAVE / LOAD SESSION
  # ═══════════════════════════════════════════════════════════

  output$save_session <- downloadHandler(
    filename = function() {
      paste0("biodiversity_session_", format(Sys.time(),"%Y%m%d_%H%M%S"), ".rds")
    },
    content = function(file) {
      if (is.null(rv$nodes_df)) {
        showNotification("Nothing to save — run a query first.", type="error")
        return()
      }
      session_data <- list(
        version      = APP_VERSION,
        saved_at     = Sys.time(),
        groups       = rv$groups,
        taxonomy     = rv$taxonomy,
        globi_raw    = rv$globi_raw,
        occurrences  = rv$occurrences,
        query_log    = rv$query_log,
        edges_df     = rv$edges_df,
        nodes_df     = rv$nodes_df,
        multi_metrics= rv$multi_metrics,
        phylo_tree   = rv$phylo_tree
      )
      saveRDS(session_data, file)
      showNotification("✅ Session saved.", type="message", duration=4)
    }
  )

  observeEvent(input$load_session, {
    req(input$load_session)
    path <- input$load_session$datapath
    dat  <- tryCatch(readRDS(path), error=function(e){
      showNotification(paste("❌ Cannot read file:", conditionMessage(e)), type="error")
      NULL
    })
    if (is.null(dat)) return()

    # Validate
    required <- c("nodes_df","edges_df","groups","taxonomy")
    missing  <- setdiff(required, names(dat))
    if (length(missing)>0) {
      showNotification(paste("❌ Invalid session file. Missing:", paste(missing,collapse=", ")),
                       type="error"); return()
    }

    # Populate reactive values
    rv$groups       <- dat$groups       %||% list()
    rv$taxonomy     <- dat$taxonomy     %||% list()
    rv$globi_raw    <- dat$globi_raw
    rv$occurrences  <- dat$occurrences  %||% list()
    rv$query_log    <- dat$query_log    %||% list()
    rv$edges_df     <- dat$edges_df
    rv$nodes_df     <- dat$nodes_df
    rv$multi_metrics<- dat$multi_metrics
    rv$phylo_tree   <- dat$phylo_tree

    # Rebuild group panels from loaded groups
    grp_names <- names(rv$groups)
    if (length(grp_names)>0) {
      new_ids <- seq_along(grp_names)
      rv$active_ids <- new_ids
      rv$group_ids  <- max(new_ids)
      rv$n_groups   <- length(new_ids)
      # Update inputs once UI has rendered
      for (i in seq_along(grp_names)) {
        g <- rv$groups[[grp_names[i]]]
        updateTextInput(session,   paste0("gname_",i),    value=grp_names[i])
        updateColourInput(session, paste0("gcol_",i),     value=g$color%||%"#607d8b")
        updateTextAreaInput(session, paste0("gspecies_",i),
                            value=paste(g$species, collapse="\n"))
      }
    }

    n_sp  <- if (!is.null(rv$nodes_df)) nrow(rv$nodes_df) else 0
    n_eg  <- if (!is.null(rv$edges_df)) nrow(rv$edges_df) else 0
    n_grp <- length(rv$groups)
    saved <- if (!is.null(dat$saved_at)) format(dat$saved_at,"%Y-%m-%d %H:%M") else "?"
    rv$session_meta <- list(saved_at=saved, n_sp=n_sp, n_edges=n_eg, n_groups=n_grp,
                             version=dat$version%||%"?")

    showNotification(
      sprintf("✅ Session loaded: %d species, %d groups, %d edges.", n_sp, n_grp, n_eg),
      type="message", duration=6)
  })

  output$session_info_sidebar <- renderUI({
    m <- rv$session_meta
    if (is.null(m)) return(NULL)
    div(class="session-info",
      tags$b("Loaded session"), tags$br(),
      sprintf("Saved: %s", m$saved_at), tags$br(),
      sprintf("Species: %d | Groups: %d", m$n_sp, m$n_groups), tags$br(),
      sprintf("Edges: %d | v%s", m$n_edges, m$version)
    )
  })

  # ═══════════════════════════════════════════════════════════
  # MAIN QUERY
  # ═══════════════════════════════════════════════════════════
  observeEvent(input$run_query, {
    groups_raw <- read_groups_from_input()
    if (length(groups_raw)==0) {
      showNotification("No species entered in any group.", type="error"); return()
    }

    # Resolve to flat species + group membership table
    group_resolved <- resolve_groups(groups_raw)
    if (nrow(group_resolved)==0) {
      showNotification("No valid species found.", type="error"); return()
    }

    dupes <- group_resolved$name[group_resolved$user_group=="Multiple"]
    if (length(dupes)>0)
      showNotification(
        paste("⚠️ Species in multiple groups (shown in grey):",
              paste(dupes,collapse=", ")),
        type="warning", duration=8)

    species <- group_resolved$name
    rv$groups       <- groups_raw
    rv$query_log    <- setNames(lapply(species, function(s) "⏳ Queued"), species)
    rv$taxonomy     <- list()
    rv$globi_raw    <- NULL
    rv$occurrences  <- list()
    rv$edges_df     <- NULL
    rv$nodes_df     <- NULL
    rv$multi_metrics<- NULL
    rv$phylo_tree   <- NULL
    rv$session_meta <- NULL

    withProgress(message="Querying databases…", value=0, {
      n <- length(species)

      # ── GBIF taxonomy ──────────────────────────────────────
      if ("gbif" %in% input$query_sources) {
        for (i in seq_along(species)) {
          sp <- species[i]
          incProgress(0.3/n, detail=paste("GBIF:", sp))
          rv$query_log[[sp]] <- "🌐 GBIF…"
          tax <- gbif_resolve(sp)
          rv$taxonomy[[sp]] <- if (!is.null(tax)) {
            tax$functional_group <- infer_functional_group(
              tax$kingdom, tax$class, tax$order, tax$family)
            rv$query_log[[sp]] <- "✅ GBIF OK"
            tax
          } else {
            rv$query_log[[sp]] <- "⚠️ GBIF not found"
            list(accepted_name=sp, kingdom=NA, phylum=NA, class=NA,
                 order=NA, family=NA, genus=NA, gbif_key=NA,
                 functional_group="Unknown")
          }
          Sys.sleep(0.3)    # ← HERE, after the closing brace of the if/else
        }                   # ← this closes the for loop        }
      } else {
        for (sp in species) {
          rv$taxonomy[[sp]] <- list(accepted_name=sp,kingdom=NA,phylum=NA,
            class=NA,order=NA,family=NA,genus=NA,gbif_key=NA,functional_group="Unknown")
          rv$query_log[[sp]] <- "⏭ Skipped GBIF"
        }
      }

      # ── Build nodes_df ────────────────────────────────────
      nodes_df <- bind_rows(lapply(seq_along(species), function(i) {
        sp  <- species[i]
        tax <- rv$taxonomy[[sp]]
        data.frame(
          id               = i,
          name             = sp,
          kingdom          = tax$kingdom           %||% NA_character_,
          class            = tax$class             %||% NA_character_,
          order            = tax$order             %||% NA_character_,
          family           = tax$family            %||% NA_character_,
          functional_group = tax$functional_group  %||% "Unknown",
          gbif_key         = as.character(tax$gbif_key %||% NA_character_),
          stringsAsFactors = FALSE
        )
      })) %>%
        left_join(group_resolved, by="name")
      rv$nodes_df <- nodes_df

      # ── GLOBI ─────────────────────────────────────────────
      globi_all <- list()
      if ("globi" %in% input$query_sources) {
        for (i in seq_along(species)) {
          sp <- species[i]
          incProgress(0.4/n, detail=paste("GLOBI:", sp))
          rv$query_log[[sp]] <- paste(rv$query_log[[sp]], "| 🌐 GLOBI…")
          gi <- globi_interactions(sp, limit=input$globi_limit)
          if (!is.null(gi) && nrow(gi) > 0) {
            globi_all[[sp]] <- gi
            rv$query_log[[sp]] <- paste(rv$query_log[[sp]],
                                        sprintf("✅ %d interactions", nrow(gi)))
          } else {
            rv$query_log[[sp]] <- paste(rv$query_log[[sp]], "⚠️ No GLOBI")
          }
          Sys.sleep(0.5)    # ← HERE, after the closing brace of the if/else
        }                   # ← this closes the for loop
      }
      rv$globi_raw <- if (length(globi_all)>0) bind_rows(globi_all) else NULL

      # ── GBIF occurrences ──────────────────────────────────
      if ("occ" %in% input$query_sources) {
        for (i in seq_along(species)) {
          sp  <- species[i]
          incProgress(0.3/n, detail=paste("Occ:", sp))
          occ <- gbif_occurrences(sp, limit=input$gbif_occ_limit)
          if (!is.null(occ)) {
            rv$occurrences[[sp]] <- occ
            rv$query_log[[sp]] <- paste(rv$query_log[[sp]],
                                        sprintf("| 📍 %d occ.", nrow(occ)))
          }
          Sys.sleep(0.4)    # ← HERE, after the closing brace of the if block
        }                   # ← this closes the for loop
      }

      # ── Build edge table ──────────────────────────────────
      edge_list <- list()
      if (!is.null(rv$globi_raw) && nrow(rv$globi_raw)>0) {
        g <- rv$globi_raw
        if (all(c("source_taxon_name","target_taxon_name","interaction_type") %in% names(g))) {
          mask <- g$source_taxon_name %in% species | g$target_taxon_name %in% species
          g2   <- g[mask, ]
          if (nrow(g2)>0)
            edge_list[["globi"]] <- data.frame(
              from_name  = g2$source_taxon_name,
              to_name    = g2$target_taxon_name,
              layer      = sapply(g2$interaction_type, classify_interaction),
              confidence = "GLOBI", source="GLOBI", stringsAsFactors=FALSE)
        }
      }

      heur <- generate_heuristic_edges(nodes_df)
      if (!is.null(heur) && nrow(heur)>0) edge_list[["heuristic"]] <- heur

      all_edges <- if (length(edge_list)>0) bind_rows(edge_list) else NULL
      if (!is.null(all_edges) && nrow(all_edges)>0) {
        all_edges <- all_edges %>%
          filter(from_name %in% species, to_name %in% species, from_name!=to_name) %>%
          distinct(from_name, to_name, layer, .keep_all=TRUE)
        if (nrow(all_edges)==0) all_edges <- NULL
      }
      rv$edges_df <- all_edges

      # ── Phylogeny ─────────────────────────────────────────
      rv$phylo_tree <- build_newick_from_taxonomy(nodes_df)

      # ── Multilayer analysis ────────────────────────────────
      if (!is.null(all_edges) && nrow(all_edges)>0) {
        rv$multi_metrics <- tryCatch(
          compute_multilayer_metrics(all_edges, nodes_df),
          error=function(e){ message("multilayer: ",conditionMessage(e)); NULL })
      }
    })

    n_edges <- if (!is.null(rv$edges_df)) nrow(rv$edges_df) else 0
    showNotification(
      sprintf("Done: %d species across %d groups, %d interactions.",
              length(species), length(groups_raw), n_edges),
      type="message", duration=6)
  })

  # ═══════════════════════════════════════════════════════════
  # TAB 1: STATUS
  # ═══════════════════════════════════════════════════════════
  output$query_status_ui <- renderUI({
    if (!length(rv$query_log))
      return(p("Press '🔍 Query All Sources' or load a session.",
               style="color:#8b949e;"))
    rows <- lapply(names(rv$query_log), function(sp) {
      msg <- rv$query_log[[sp]]
      cls <- if (grepl("✅",msg,fixed=TRUE)) "status-ok" else "status-wait"
      grp <- if (!is.null(rv$nodes_df)) {
        nd <- rv$nodes_df[rv$nodes_df$name==sp, ]
        if (nrow(nd)>0) nd$user_group[1] else ""
      } else ""
      tags$tr(tags$td(tags$em(sp)), tags$td(grp), tags$td(tags$span(msg,class=cls)))
    })
    tags$table(class="sp-status-tbl",
      tags$thead(tags$tr(tags$th("Species"),tags$th("Group"),tags$th("Status"))),
      do.call(tags$tbody, rows))
  })

  output$tax_preview <- renderDT({
    req(rv$nodes_df)
    rv$nodes_df %>%
      select(name,user_group,kingdom,class,order,family,functional_group) %>%
      datatable(options=list(pageLength=8,dom="ftp"), rownames=FALSE,
                colnames=c("Species","Group","Kingdom","Class","Order","Family","Func. Group"))
  })

  output$globi_raw_table <- renderDT({
    req(rv$globi_raw)
    rv$globi_raw %>%
      select(any_of(c("source_taxon_name","interaction_type","target_taxon_name","query_species"))) %>%
      datatable(options=list(pageLength=10,dom="lftp",scrollX=TRUE), rownames=FALSE)
  })

  # ═══════════════════════════════════════════════════════════
  # TAB 2: NETWORK
  # ═══════════════════════════════════════════════════════════
  output$layer_legend_ui <- renderUI({
    layers_present <- if (!is.null(rv$edges_df) && nrow(rv$edges_df)>0)
      unique(rv$edges_df$layer) else character(0)
    div(style="padding:4px 0;",
      lapply(layers_present, function(l) {
        cfg <- LAYERS[[l]]; if (is.null(cfg)) cfg <- list(label=l,color="#888",icon="●")
        tags$span(class="layer-badge",
          style=sprintf("background:%s22;color:%s;border:1px solid %s;",
                        cfg$color,cfg$color,cfg$color),
          cfg$icon," ",cfg$label)
      }),
      tags$span("  ● animal  ■ plant",style="color:#555;font-size:11px;margin-left:10px;")
    )
  })

  output$group_legend_ui <- renderUI({
    req(rv$groups, length(rv$groups)>0)
    cb <- input$color_by %||% "functional_group"
    if (cb != "user_group") return(NULL)
    div(style="padding:4px 0;margin-top:2px;",
      tags$span("User groups: ", style="color:#8b949e;font-size:11px;"),
      lapply(names(rv$groups), function(gn) {
        col <- rv$groups[[gn]]$color %||% "#888"
        tags$span(class="layer-badge",
          style=sprintf("background:%s22;color:%s;border:1px solid %s;",col,col,col),
          gn)
      }),
      tags$span(class="layer-badge",
        style="background:#aaaaaa22;color:#aaaaaa;border:1px solid #aaaaaa;",
        "Multiple")
    )
  })

  get_vis_data <- reactive({
    req(rv$nodes_df, rv$edges_df)
    nodes  <- rv$nodes_df
    cb     <- input$color_by %||% "functional_group"
    colors <- get_node_colors(nodes, cb)

    nodes$color.background <- paste0(colors[nodes$name], "55")
    nodes$color.border      <- colors[nodes$name]
    nodes$color.highlight.background <- colors[nodes$name]
    nodes$shape      <- ifelse(!is.na(nodes$kingdom)&nodes$kingdom=="Plantae","square","dot")
    nodes$size       <- 22
    nodes$label      <- gsub(" ","\n",nodes$name)
    nodes$font.color <- "#ffffff"
    nodes$font.size  <- 12
    nodes$title <- paste0(
      "<b>",nodes$name,"</b><br>",
      "<i>",coalesce_vec(nodes$family,"?")," · ",coalesce_vec(nodes$order,"?"),"</i><br>",
      "Group: ",coalesce_vec(nodes$user_group,"—"),"<br>",
      nodes$functional_group,"<br>Kingdom: ",coalesce_vec(nodes$kingdom,"?")
    )

    edges <- rv$edges_df
    if (!is.null(input$conf_filter) && "confidence" %in% names(edges))
      edges <- edges[edges$confidence %in% input$conf_filter, ]
    if (!is.null(input$net_view) && input$net_view=="single" && !is.null(input$net_layer))
      edges <- edges[edges$layer==input$net_layer, ]

    name_to_id <- setNames(nodes$id, nodes$name)
    edges$from <- name_to_id[edges$from_name]
    edges$to   <- name_to_id[edges$to_name]
    edges <- edges[!is.na(edges$from) & !is.na(edges$to), ]

    edges$color  <- sapply(edges$layer,function(l){ cfg<-LAYERS[[l]]; if(is.null(cfg))"#888" else cfg$color })
    edges$dashes <- grepl("Low|Heuristic-Low", edges$confidence)
    edges$width  <- ifelse(grepl("GLOBI",edges$confidence),2.5,
                    ifelse(grepl("High|Moderate",edges$confidence),2,1.2))
    edges$arrows <- "to"
    src_col <- if ("source" %in% names(edges)) edges$source else rep("",nrow(edges))
    edges$title <- paste0("<b>",edges$layer,"</b><br>From: ",edges$from_name,
                          "<br>To: ",edges$to_name,"<br>Source: ",src_col)
    list(nodes=nodes, edges=edges)
  })

  output$main_network <- renderVisNetwork({
    vd <- get_vis_data()
    visNetwork(vd$nodes, vd$edges, background="#0d1117") %>%
      visNodes(shadow=list(enabled=TRUE,size=10), font=list(color="#ffffff",size=12)) %>%
      visEdges(smooth=list(type="dynamic")) %>%
      visOptions(highlightNearest=list(enabled=TRUE,degree=1,hover=TRUE),
                 selectedBy="user_group", nodesIdSelection=TRUE) %>%
      visPhysics(solver="forceAtlas2Based",
                 forceAtlas2Based=list(gravitationalConstant=-80,centralGravity=0.01,
                                       springLength=140,springConstant=0.08),
                 stabilization=list(iterations=300)) %>%
      visLayout(randomSeed=42) %>%
      visInteraction(hover=TRUE, navigationButtons=TRUE, keyboard=TRUE) %>%
      visEvents(
        selectNode  ="function(n){ Shiny.setInputValue('sel_node',n.nodes[0]); }",
        deselectNode="function(){  Shiny.setInputValue('sel_node',null); }"
      )
  })

  output$net_detail_ui <- renderUI({
    if (!is.null(input$sel_node) && length(input$sel_node)>0 && !is.na(input$sel_node)) {
      sp <- rv$nodes_df[rv$nodes_df$id==as.integer(input$sel_node), ]
      if (nrow(sp)==0) return(NULL)
      sp_name <- sp$name[1]
      ug      <- sp$user_group[1] %||% "—"
      ug_col  <- sp$user_group_color[1] %||% "#607d8b"
      fg_col  <- safe_color(FUNC_COLORS, sp$functional_group[1])
      active_layers <- if (!is.null(rv$edges_df))
        unique(rv$edges_df$layer[rv$edges_df$from_name==sp_name|rv$edges_df$to_name==sp_name])
      else character(0)
      tagList(
        tags$h4(tags$em(sp_name), style=paste0("color:",ug_col,";")),
        tags$table(style="font-size:13px;width:100%;",
          tags$tr(tags$td(tags$strong("User group")),
                  tags$td(tags$span(ug, style=paste0("color:",ug_col,";font-weight:bold;")))),
          tags$tr(tags$td(tags$strong("Kingdom")), tags$td(sp$kingdom[1])),
          tags$tr(tags$td(tags$strong("Order")),   tags$td(sp$order[1])),
          tags$tr(tags$td(tags$strong("Family")),  tags$td(sp$family[1])),
          tags$tr(tags$td(tags$strong("Function")),tags$td(sp$functional_group[1]))
        ),
        if (length(active_layers)>0) div(style="margin-top:6px;",
          lapply(active_layers, function(l) {
            cfg <- LAYERS[[l]]; if (is.null(cfg)) return(NULL)
            tags$span(class="layer-badge",
              style=sprintf("background:%s22;color:%s;border:1px solid %s;",
                            cfg$color,cfg$color,cfg$color),
              cfg$icon," ",cfg$label)
          })
        )
      )
    } else tags$p("Click a node for details.", style="color:#8b949e;")
  })

  output$layer_degree_chart <- renderPlotly({
    req(rv$multi_metrics)
    df <- rv$multi_metrics$layer_metrics %>%
      filter(degree>0) %>%
      mutate(layer_label=sapply(layer,function(l){cfg<-LAYERS[[l]];if(is.null(cfg))l else cfg$label}))
    if (nrow(df)==0) return(empty_plotly("No interactions with degree > 0"))
    label_colors <- setNames(sapply(names(LAYERS),function(l)LAYERS[[l]]$color),
                             sapply(names(LAYERS),function(l)LAYERS[[l]]$label))
    p <- ggplot(df,aes(x=reorder(species,degree),y=degree,fill=layer_label,
                       text=paste0(species,"\n",layer_label,": ",degree))) +
      geom_col(position="stack",width=0.7) +
      scale_fill_manual(values=label_colors) +
      coord_flip() + labs(x=NULL,y="Degree",fill="Layer") + theme_minimal(base_size=11) +
      theme(panel.background=element_rect(fill="#0d1117",color=NA),
            plot.background =element_rect(fill="#0d1117",color=NA),
            text=element_text(color="#c9d1d9"), axis.text=element_text(color="#c9d1d9",size=9),
            panel.grid=element_line(color="#21262d"),
            legend.background=element_rect(fill="#0d1117"))
    ggplotly(p,tooltip="text") %>%
      layout(paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
             font=list(color="#c9d1d9"),margin=list(l=120))
  })

  # ═══════════════════════════════════════════════════════════
  # TAB 3: MULTILAYER ANALYSIS
  # ═══════════════════════════════════════════════════════════

  # Helper: get color vector for a given color_by across all species in summary
  get_multi_colors <- function(summary_df, color_by) {
    if (color_by=="user_group" && "user_group" %in% names(summary_df) &&
        !is.null(rv$nodes_df)) {
      ug_colors <- setNames(rv$nodes_df$user_group_color, rv$nodes_df$name)
      sapply(summary_df$species, function(sp) ug_colors[sp] %||% "#607d8b")
    } else {
      sapply(summary_df$functional_group %||% rep("Unknown",nrow(summary_df)),
             function(fg) safe_color(FUNC_COLORS,fg))
    }
  }

  output$multi_summary_cards <- renderUI({
    req(rv$multi_metrics)
    m      <- rv$multi_metrics
    n_lay  <- length(m$layer_graphs)
    n_edg  <- if (!is.null(rv$edges_df)) nrow(rv$edges_df) else 0
    P      <- m$participation
    top_p  <- if (!is.null(P)&&any(!is.na(P))) names(which.max(P)) else "—"
    D      <- m$aggregate_degree
    hub    <- if (!is.null(D)&&length(D)>0) names(which.max(D)) else "—"
    n_grp  <- length(rv$groups)
    fluidRow(
      column(2, div(class="metric-card", div(class="val",n_lay),         div(class="lbl","Layers"))),
      column(2, div(class="metric-card", div(class="val",n_edg),         div(class="lbl","Edges"))),
      column(2, div(class="metric-card", div(class="val",n_grp),         div(class="lbl","User groups"))),
      column(3, div(class="metric-card", div(class="val",style="font-size:.9rem;",top_p), div(class="lbl","Top participation"))),
      column(3, div(class="metric-card", div(class="val",style="font-size:.9rem;",hub),   div(class="lbl","Aggregate hub")))
    )
  })

  output$participation_plot <- renderPlotly({
    req(rv$multi_metrics)
    P <- rv$multi_metrics$participation
    if (is.null(P)||all(is.na(P)))
      return(empty_plotly("Need ≥ 2 layers for participation coefficient"))
    cb  <- input$multi_color_by %||% "functional_group"
    df  <- data.frame(species=names(P), P=as.numeric(P), stringsAsFactors=FALSE) %>%
      left_join(rv$nodes_df %>% select(name,functional_group,user_group,user_group_color),
                by=c("species"="name")) %>% arrange(desc(P))
    col_vec <- get_multi_colors(df, cb)
    df$color_group <- if (cb=="user_group") coalesce_vec(df$user_group,"—")
                      else coalesce_vec(df$functional_group,"Unknown")
    # Build distinct color map for scale_fill_manual
    cmap <- setNames(col_vec, df$color_group)[!duplicated(df$color_group)]

    p <- ggplot(df,aes(x=reorder(species,P),y=P,fill=color_group,
                       text=paste0(species,"\nP = ",round(P,3)))) +
      geom_col(width=0.7) +
      geom_hline(yintercept=0.5,linetype="dashed",color="#58a6ff",alpha=0.5) +
      coord_flip() + scale_y_continuous(limits=c(0,1)) +
      scale_fill_manual(values=cmap,na.value="#607d8b") +
      labs(x=NULL,y="Participation Coefficient (P)",fill=if(cb=="user_group")"Group" else "Func. Group") +
      theme_minimal(base_size=11) +
      theme(panel.background=element_rect(fill="#0d1117",color=NA),
            plot.background =element_rect(fill="#0d1117",color=NA),
            text=element_text(color="#c9d1d9"),axis.text=element_text(color="#c9d1d9",size=9),
            panel.grid=element_line(color="#21262d"),
            legend.background=element_rect(fill="#0d1117"))
    ggplotly(p,tooltip="text") %>%
      layout(paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
             font=list(color="#c9d1d9"),margin=list(l=140))
  })

  output$layer_overlap_plot <- renderPlotly({
    req(rv$multi_metrics)
    mat <- rv$multi_metrics$layer_overlap
    if (is.null(mat)) return(empty_plotly("Need ≥ 2 layers for overlap matrix"))
    lbls <- sapply(rownames(mat),function(l){cfg<-LAYERS[[l]];if(is.null(cfg))l else cfg$label})
    plot_ly(z=mat,x=lbls,y=lbls,type="heatmap",colorscale="Blues",zmin=0,zmax=1,
            hovertemplate="<b>%{y} ↔ %{x}</b><br>Jaccard: %{z:.3f}<extra></extra>") %>%
      layout(paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
             font=list(color="#c9d1d9"),xaxis=list(tickangle=-30),
             margin=list(b=100,l=100))
  })

  output$presence_heatmap <- renderPlotly({
    req(rv$multi_metrics)
    mat <- rv$multi_metrics$layer_presence
    if (is.null(mat)) return(empty_plotly("No layer presence data"))
    lbls <- sapply(colnames(mat),function(l){cfg<-LAYERS[[l]];if(is.null(cfg))l else cfg$label})
    z <- matrix(as.integer(mat),nrow=nrow(mat),ncol=ncol(mat),dimnames=dimnames(mat))
    plot_ly(z=z,x=lbls,y=rownames(z),type="heatmap",
            colorscale=list(c(0,"#0d1117"),c(1,"#58a6ff")),
            showscale=FALSE,zmin=0,zmax=1,
            hovertemplate="<b>%{y}</b><br>%{x}: %{z}<extra></extra>") %>%
      layout(paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
             font=list(color="#c9d1d9"),margin=list(l=160,b=100),
             xaxis=list(tickangle=-30))
  })

  output$centrality_plot <- renderPlotly({
    req(rv$multi_metrics)
    m <- rv$multi_metrics
    if (is.null(m$aggregate_degree)) return(empty_plotly("No aggregate network data"))
    all_names <- names(m$aggregate_degree)
    df <- data.frame(
      species     = all_names,
      degree      = as.numeric(m$aggregate_degree),
      betweenness = as.numeric(m$aggregate_betweenness %||% rep(0,length(all_names))),
      clustering  = as.numeric(m$aggregate_clustering  %||% rep(0,length(all_names))),
      stringsAsFactors=FALSE
    ) %>%
      left_join(rv$nodes_df %>% select(name,functional_group,user_group,user_group_color),
                by=c("species"="name"))

    # Always colour by user group if available, otherwise functional
    has_ug <- "user_group" %in% names(df) && any(!is.na(df$user_group))
    col_vec <- if (has_ug && length(rv$groups)>0)
      setNames(coalesce_vec(df$user_group_color,"#607d8b"), df$species)
    else sapply(df$functional_group,function(fg)safe_color(FUNC_COLORS,fg))

    df$color <- col_vec
    df$color_group <- if (has_ug && length(rv$groups)>0)
      coalesce_vec(df$user_group,"—") else coalesce_vec(df$functional_group,"Unknown")

    plot_ly(df, x=~degree, y=~betweenness, size=~(clustering+0.1)*20,
            color=~color_group, colors=setNames(df$color,df$color_group),
            text=~paste0("<b>",species,"</b><br>Group: ",color_group,
                         "<br>Degree: ",degree,"<br>Betweenness: ",round(betweenness,3),
                         "<br>Clustering: ",round(clustering,3)),
            hoverinfo="text", type="scatter", mode="markers",
            marker=list(opacity=0.85,sizemode="diameter")) %>%
      layout(xaxis=list(title="Aggregate Degree",gridcolor="#21262d"),
             yaxis=list(title="Betweenness Centrality",gridcolor="#21262d"),
             paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
             font=list(color="#c9d1d9"),legend=list(bgcolor="#0d1117"))
  })

  output$multi_table <- renderDT({
    req(rv$multi_metrics)
    df <- rv$multi_metrics$species_summary %>%
      mutate(across(where(is.numeric),~round(.,4)))
    datatable(df,options=list(pageLength=10,dom="lftp",scrollX=TRUE),rownames=FALSE)
  })

  # ═══════════════════════════════════════════════════════════
  # TAB 4: BIPARTITE
  # ═══════════════════════════════════════════════════════════
  bp_data <- reactive({
    req(rv$edges_df, rv$nodes_df, input$bp_layer)
    compute_bipartite_metrics(rv$edges_df, rv$nodes_df, input$bp_layer)
  })

  output$bp_net_metrics <- renderUI({
    bpm <- bp_data()
    if (is.null(bpm))
      return(p("Insufficient data. Need ≥2 plants and ≥2 animals with interactions.",
               style="color:#8b949e;"))
    nm <- bpm$network_metrics
    div(style="display:flex;flex-wrap:wrap;gap:12px;",
      lapply(names(nm), function(n)
        div(class="metric-card",style="min-width:120px;",
          div(class="val",style="font-size:1.2rem;",round(nm[[n]],4)),
          div(class="lbl",n))))
  })

  output$bp_web_plot <- renderPlotly({
    bpm <- bp_data()
    if (is.null(bpm)) return(empty_plotly("Insufficient bipartite data"))
    mat <- as.matrix(bpm$web)
    plot_ly(z=mat,x=colnames(mat),y=rownames(mat),type="heatmap",colorscale="Viridis",
            hovertemplate="Plant: %{y}<br>Animal: %{x}<br>Interactions: %{z}<extra></extra>") %>%
      layout(xaxis=list(title="Animals",tickangle=-30),yaxis=list(title="Plants"),
             paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
             font=list(color="#c9d1d9"),margin=list(b=120,l=140))
  })

  output$bp_species_table <- renderDT({
    bpm <- bp_data()
    if (is.null(bpm)||is.null(bpm$species_metrics)) return(datatable(data.frame()))
    sm <- bpm$species_metrics
    df <- tryCatch({
      lower <- as.data.frame(sm[[1]],stringsAsFactors=FALSE)
      lower$species <- rownames(lower); lower$trophic <- "Plant"
      upper <- as.data.frame(sm[[2]],stringsAsFactors=FALSE)
      upper$species <- rownames(upper); upper$trophic <- "Animal"
      bind_rows(lower,upper)
    }, error=function(e) as.data.frame(sm,stringsAsFactors=FALSE))
    df %>% mutate(across(where(is.numeric),~round(.,3))) %>%
      datatable(options=list(pageLength=10,dom="ftp",scrollX=TRUE),rownames=FALSE)
  })

  # ═══════════════════════════════════════════════════════════
  # TAB 5: PHYLOGENY
  # ═══════════════════════════════════════════════════════════
  output$phylo_plot <- renderPlotly({
    req(rv$phylo_tree, rv$nodes_df)
    tree   <- rv$phylo_tree
    n_tips <- length(tree$tip.label)
    if (n_tips<2) return(empty_plotly("Need ≥2 species with known taxonomy"))
    H <- max(480, n_tips*40)
    pp <- tryCatch({
      ape::plot.phylo(tree,type="phylogram",direction="rightwards",
                      plot=FALSE,use.edge.length=FALSE)
      get("last_plot.phylo",envir=.PlotPhyloEnv)
    }, error=function(e) NULL)
    if (is.null(pp)) return(empty_plotly("Could not compute phylogenetic coordinates"))

    xc <- pp$xx; yc <- pp$yy
    edge_x<-c(); edge_y<-c(); vert_x<-c(); vert_y<-c()
    for (i in seq_len(nrow(tree$edge))) {
      pr<-tree$edge[i,1]; ch<-tree$edge[i,2]
      edge_x<-c(edge_x,xc[pr],xc[ch],NA); edge_y<-c(edge_y,yc[ch],yc[ch],NA)
    }
    for (nd in unique(tree$edge[,1])) {
      ch<-tree$edge[tree$edge[,1]==nd,2]
      vert_x<-c(vert_x,xc[nd],xc[nd],NA); vert_y<-c(vert_y,min(yc[ch]),max(yc[ch]),NA)
    }

    cb <- input$phylo_color_by %||% "functional_group"
    tip_df <- data.frame(tip_label=tree$tip.label,
                         x=xc[seq_len(n_tips)],y=yc[seq_len(n_tips)],
                         stringsAsFactors=FALSE) %>%
      mutate(display=gsub("_"," ",tip_label)) %>%
      left_join(rv$nodes_df %>% select(name,functional_group,family,
                                        user_group,user_group_color),
                by=c("display"="name"))

    tip_df$color <- if (cb=="user_group") {
      coalesce_vec(tip_df$user_group_color, "#607d8b")
    } else {
      sapply(coalesce_vec(tip_df$functional_group,"Unknown"),
             function(fg) safe_color(FUNC_COLORS,fg))
    }
    tip_df$legend_group <- if (cb=="user_group")
      coalesce_vec(tip_df$user_group,"—")
    else coalesce_vec(tip_df$functional_group,"Unknown")

    fig <- plot_ly(height=H) %>%
      add_trace(x=edge_x,y=edge_y,type="scatter",mode="lines",
                line=list(color="#3a5a7a",width=1.5),hoverinfo="none",
                showlegend=FALSE,name="") %>%
      add_trace(x=vert_x,y=vert_y,type="scatter",mode="lines",
                line=list(color="#3a5a7a",width=1.5),hoverinfo="none",
                showlegend=FALSE,name="") %>%
      add_trace(data=tip_df,x=~x,y=~y,type="scatter",mode="markers+text",
                marker=list(size=14,color=~color,line=list(color="#fff",width=1.5)),
                text=~display,textposition="middle right",
                textfont=list(color="#c9d1d9",size=11),
                hovertemplate=paste0("<b>%{text}</b><br>Group: ",
                  coalesce_vec(tip_df$legend_group,"?"),"<extra></extra>"),
                showlegend=FALSE)

    groups_shown <- unique(tip_df$legend_group); groups_shown <- groups_shown[!is.na(groups_shown)]
    for (g in groups_shown) {
      gcol <- tip_df$color[tip_df$legend_group==g][1] %||% "#607d8b"
      fig <- fig %>% add_trace(x=NA,y=NA,type="scatter",mode="markers",
                               marker=list(size=10,color=gcol),name=g,showlegend=TRUE)
    }
    fig %>% layout(paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
                   xaxis=list(visible=FALSE),yaxis=list(visible=FALSE),
                   margin=list(l=10,r=220,t=20,b=20),font=list(color="#c9d1d9"))
  })

  # ═══════════════════════════════════════════════════════════
  # TAB 6: MAP
  # ═══════════════════════════════════════════════════════════
  MAP_PAL <- c("#00f5c4","#f5a623","#ff6b6b","#9b59b6","#3498db",
               "#e67e22","#1abc9c","#e74c3c","#2ecc71","#f39c12",
               "#8e44ad","#27ae60","#d35400","#16a085","#c0392b")

  output$map_sp_picker <- renderUI({
    req(length(rv$occurrences)>0)
    pickerInput("map_species","Species:",choices=names(rv$occurrences),
                selected=names(rv$occurrences),multiple=TRUE,
                options=list(`actions-box`=TRUE,`live-search`=TRUE))
  })

  output$occ_map <- renderLeaflet({
    leaflet() %>% addProviderTiles(providers$CartoDB.DarkMatter) %>%
      setView(lng=10,lat=30,zoom=2)
  })

  observeEvent(input$refresh_map, {
    req(input$map_species, length(rv$occurrences)>0)
    leafletProxy("occ_map") %>% clearMarkers() %>% clearControls()
    selected <- input$map_species
    for (i in seq_along(selected)) {
      sp  <- selected[i]; occ <- rv$occurrences[[sp]]
      if (is.null(occ)||nrow(occ)==0) next
      col <- MAP_PAL[(i-1)%%length(MAP_PAL)+1]
      leafletProxy("occ_map") %>%
        addCircleMarkers(data=occ,lng=~decimalLongitude,lat=~decimalLatitude,
                         color=col,radius=4,stroke=FALSE,fillOpacity=0.75,
                         popup=paste0("<b>",sp,"</b>"))
    }
    leafletProxy("occ_map") %>%
      addLegend("bottomright",colors=MAP_PAL[seq_along(selected)],
                labels=selected,title="Species",opacity=0.9)
  })

  # ═══════════════════════════════════════════════════════════
  # TAB 7: SPECIES INFO
  # ═══════════════════════════════════════════════════════════
  output$info_selector <- renderUI({
    req(rv$nodes_df)
    selectInput("info_sp","Select species:",
                choices=setNames(rv$nodes_df$name,rv$nodes_df$name))
  })

  output$species_card_ui <- renderUI({
    req(rv$nodes_df, input$info_sp)
    sp <- rv$nodes_df[rv$nodes_df$name==input$info_sp, ]
    if (nrow(sp)==0) return(NULL)
    sp_name <- sp$name[1]
    fg      <- as.character(sp$functional_group[1]) %||% "Unknown"
    ug      <- as.character(sp$user_group[1]) %||% "—"
    ug_col  <- as.character(sp$user_group_color[1]) %||% "#607d8b"
    n_out   <- if (!is.null(rv$edges_df)) sum(rv$edges_df$from_name==sp_name) else 0
    n_in    <- if (!is.null(rv$edges_df)) sum(rv$edges_df$to_name==sp_name)   else 0
    tagList(div(style=paste0("border-left:4px solid ",ug_col,";padding:12px;margin-top:10px;"),
      tags$h4(tags$em(sp_name),style=paste0("color:",ug_col,";")),
      tags$table(style="font-size:13px;width:100%;",
        tags$tr(tags$td(tags$strong("User group")),
                tags$td(tags$b(ug,style=paste0("color:",ug_col,";")))),
        tags$tr(tags$td(tags$strong("Kingdom")),  tags$td(sp$kingdom[1])),
        tags$tr(tags$td(tags$strong("Class")),    tags$td(sp$class[1])),
        tags$tr(tags$td(tags$strong("Order")),    tags$td(sp$order[1])),
        tags$tr(tags$td(tags$strong("Family")),   tags$td(sp$family[1])),
        tags$tr(tags$td(tags$strong("Function")), tags$td(fg)),
        tags$tr(tags$td(tags$strong("Links out")),tags$td(n_out)),
        tags$tr(tags$td(tags$strong("Links in")), tags$td(n_in))
      ),
      div(style="margin-top:10px;display:flex;gap:8px;flex-wrap:wrap;",
        tags$a(href=paste0("https://www.gbif.org/species/search?q=",URLencode(sp_name)),
               target="_blank",class="btn btn-xs btn-info","GBIF ↗"),
        tags$a(href=paste0("https://en.wikipedia.org/wiki/",gsub(" ","_",sp_name)),
               target="_blank",class="btn btn-xs btn-default","Wikipedia ↗"),
        tags$a(href=paste0("https://www.globalbioticinteractions.org/?sourceTaxon=",URLencode(sp_name)),
               target="_blank",class="btn btn-xs btn-warning","GLOBI ↗")
      )
    ))
  })

  output$species_full_table <- renderDT({
    req(rv$nodes_df)
    df <- rv$nodes_df %>%
      select(name,user_group,kingdom,class,order,family,functional_group)
    if (!is.null(rv$multi_metrics)) {
      sm <- rv$multi_metrics$species_summary %>%
        select(species,n_layers_active,participation_coeff,agg_degree)
      df <- df %>% left_join(sm,by=c("name"="species"))
    }
    datatable(df,options=list(pageLength=12,dom="lftp",scrollX=TRUE),rownames=FALSE) %>%
      formatRound(columns=intersect(c("participation_coeff","agg_degree"),names(df)),digits=3)
  })

} # end server

shinyApp(ui=ui, server=server)
