#' @keywords internal
"_PACKAGE"

# Non-standard-evaluation column names used inside ggplot2 aes() and plotly
# formulae throughout the package. Declared here to avoid spurious
# "no visible binding for global variable" NOTEs from R CMD check / lintr.
utils::globalVariables(c(
  "PC", "PC1", "PC2", "PC3", "host_type", "isolate",
  "variance", "cumulative", "CAI"
))

# ── Null-coalescing operator ──────────────────────────────────────────────────

#' Null-coalescing operator
#' @param a Left-hand side value
#' @param b Default value returned when `a` is NULL or zero-length
#' @return `a` if non-NULL and non-empty, otherwise `b`
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# ── String concatenation ──────────────────────────────────────────────────────

#' String concatenation operator
#' @param x Character string
#' @param y Character string
#' @return `paste0(x, y)`
#' @keywords internal
`%s+%` <- function(x, y) paste0(x, y)
