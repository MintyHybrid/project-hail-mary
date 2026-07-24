# Shared plot theme matching the book's Manuscript/Cyanotype visual theme
# (theme-light.scss at the repo root). Keeps figure chrome -- background,
# gridlines, fonts -- in sync with the site's HTML/CSS without hardcoding
# the palette separately in every notebook.

#' Cyanotype colour palette
#'
#' Named hex colours matching the book's `theme-light.scss` CSS custom
#' properties, for notebooks that want to opt a `scale_*_manual()` call into
#' the site palette rather than picking ad hoc hex values. `text_strong` and
#' `text_soft` are a separate, deliberately neutral pair (see
#' [theme_manuscript()]) -- not simply `ink`/`soft` -- because static
#' R-rendered figures are baked once at render time and can't react to the
#' light/dark toggle the way the site's own CSS can; `ink`/`soft` are tuned
#' for the light paper background specifically and would be close to
#' illegible against the dark theme's paper colour.
#'
#' @return Named character vector: `paper`, `ink`, `soft`, `accent`, `flag`,
#'   `line`, `text_strong`, `text_soft`.
#' @export
hailmary_palette <- function() {
  c(
    paper  = "#EFF1EC",
    ink    = "#171B1D",
    soft   = "#4B5459",
    accent = "#1F5C73",
    flag   = "#9C5B18",
    line   = "#C9CFC9",
    text_strong = "#7A8286",
    text_soft   = "#979D9F"
  )
}

#' ggplot2 theme matching the book's Manuscript/Cyanotype visual theme
#'
#' A `theme_minimal()`-derived theme with a transparent background so figures
#' show whichever page theme (light or dark) is active, hairline gridlines,
#' and a neutral text palette chosen to stay legible against both the light
#' and dark paper colours (see [hailmary_palette()]) since a given render of
#' a figure can't itself react to the reader's theme toggle. Only themes
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
      plot.background  = ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.key       = ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.grid.major = ggplot2::element_line(colour = pal[["text_soft"]], linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text   = ggplot2::element_text(colour = pal[["text_soft"]]),
      axis.title  = ggplot2::element_text(colour = pal[["text_strong"]]),
      plot.title  = ggplot2::element_text(colour = pal[["text_strong"]], face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = pal[["text_soft"]]),
      plot.caption  = ggplot2::element_text(colour = pal[["text_soft"]]),
      strip.background = ggplot2::element_rect(fill = "transparent", colour = pal[["text_soft"]]),
      strip.text  = ggplot2::element_text(colour = pal[["text_strong"]]),
      legend.text  = ggplot2::element_text(colour = pal[["text_strong"]]),
      legend.title = ggplot2::element_text(colour = pal[["text_strong"]])
    )
}

#' Background-only overlay for chrome-free ggplot themes
#'
#' A partial `theme()` covering `plot.background`, `panel.background` and
#' the legend background/key, set to transparent so figures show whichever
#' page theme is active rather than a baked-in paper colour. Meant to be
#' layered after visualization-specific themes that intentionally hide
#' axes/grid (`ggplot2::theme_void()`, `ggtree::theme_tree()`/
#' `theme_tree2()`, `ggraph::theme_graph()`, `ggseqlogo::theme_logo()`) --
#' those set both background elements to `element_blank()` (transparent) as
#' part of their own construction, so overriding `plot.background` alone
#' still leaves the PNG device's own background (and, for plots with a
#' legend, the legend's own opaque background) showing through. Legend text
#' and title are also recoloured to the same neutral `text_strong` used by
#' [theme_manuscript()], since their default black would be illegible once
#' the legend background stops being paper-coloured on the dark theme.
#' Use `theme_manuscript()` instead for ordinary axis-bearing plots.
#'
#' @return A ggplot2 theme object (partial).
#' @export
theme_manuscript_bg <- function() {
  pal <- hailmary_palette()
  ggplot2::theme(
    plot.background  = ggplot2::element_rect(fill = "transparent", colour = NA),
    panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    legend.key       = ggplot2::element_rect(fill = "transparent", colour = NA),
    legend.text  = ggplot2::element_text(colour = pal[["text_strong"]]),
    legend.title = ggplot2::element_text(colour = pal[["text_strong"]])
  )
}
