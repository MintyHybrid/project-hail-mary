# Shared plot theme matching the book's Manuscript/Cyanotype visual theme
# (theme-light.scss at the repo root). Keeps figure chrome -- background,
# gridlines, fonts -- in sync with the site's HTML/CSS without hardcoding
# the palette separately in every notebook.

#' Cyanotype colour palette
#'
#' Named hex colours matching the book's `theme-light.scss` CSS custom
#' properties, for notebooks that want to opt a `scale_*_manual()` call into
#' the site palette rather than picking ad hoc hex values.
#'
#' @return Named character vector: `paper`, `ink`, `soft`, `accent`, `flag`,
#'   `line`.
#' @export
hailmary_palette <- function() {
  c(
    paper  = "#EFF1EC",
    ink    = "#171B1D",
    soft   = "#4B5459",
    accent = "#1F5C73",
    flag   = "#9C5B18",
    line   = "#C9CFC9"
  )
}

#' ggplot2 theme matching the book's Manuscript/Cyanotype visual theme
#'
#' A `theme_minimal()`-derived theme using the same paper background, ink
#' text colour, and hairline gridlines as `theme-light.scss`. Only themes
#' figure chrome (background, gridlines, axis/legend text, title styling) --
#' it does not touch data colours (`scale_colour_*`/`scale_fill_*`), which
#' stay whatever each plot's own encoding calls for.
#'
#' @param base_size Base font size in pt, passed to [ggplot2::theme_minimal()].
#' @return A ggplot2 theme object.
#' @export
theme_manuscript <- function(base_size = 11) {
  pal <- hailmary_palette()
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = pal[["paper"]], colour = NA),
      panel.background = ggplot2::element_rect(fill = pal[["paper"]], colour = NA),
      legend.background = ggplot2::element_rect(fill = pal[["paper"]], colour = NA),
      legend.key       = ggplot2::element_rect(fill = pal[["paper"]], colour = NA),
      panel.grid.major = ggplot2::element_line(colour = pal[["line"]], linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text   = ggplot2::element_text(colour = pal[["soft"]]),
      axis.title  = ggplot2::element_text(colour = pal[["ink"]]),
      plot.title  = ggplot2::element_text(colour = pal[["ink"]], face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = pal[["soft"]]),
      plot.caption  = ggplot2::element_text(colour = pal[["soft"]]),
      strip.background = ggplot2::element_rect(fill = pal[["line"]], colour = NA),
      strip.text  = ggplot2::element_text(colour = pal[["ink"]]),
      legend.text  = ggplot2::element_text(colour = pal[["ink"]]),
      legend.title = ggplot2::element_text(colour = pal[["ink"]])
    )
}

#' Background-only overlay for chrome-free ggplot themes
#'
#' A partial `theme()` covering just `plot.background` *and*
#' `panel.background` with the book's paper colour. Meant to be layered
#' after visualization-specific themes that intentionally hide axes/grid
#' (`ggplot2::theme_void()`, `ggtree::theme_tree()`/`theme_tree2()`,
#' `ggraph::theme_graph()`, `ggseqlogo::theme_logo()`) -- those set both
#' background elements to `element_blank()` (transparent) as part of their
#' own construction, so overriding `plot.background` alone still leaves the
#' PNG device's white background showing through the (larger) panel area.
#' Use `theme_manuscript()` instead for ordinary axis-bearing plots.
#'
#' @return A ggplot2 theme object (partial).
#' @export
theme_manuscript_bg <- function() {
  pal <- hailmary_palette()
  ggplot2::theme(
    plot.background  = ggplot2::element_rect(fill = pal[["paper"]], colour = NA),
    panel.background = ggplot2::element_rect(fill = pal[["paper"]], colour = NA)
  )
}
