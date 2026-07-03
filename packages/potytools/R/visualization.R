# Visualization helpers for codon usage analysis
# Ported from 2026_wsl/R/enhanced_visualization.R

#' Plot PCA biplot of RSCU codon usage, coloured by host type
#'
#' @param pca_results List from [perform_rscu_pca_improved()]; its `$data`
#'   already carries a `host_type` column.
#' @param host_classification Optional host classification data frame. If the
#'   PCA data lacks `host_type`, it is merged in from here by `isolate`.
#' @param interactive If `TRUE`, return a plotly object (requires plotly).
#' @return ggplot (or plotly) object.
#' @export
plot_pca_biplot <- function(pca_results, host_classification = NULL,
                            interactive = FALSE) {
  pca_data <- pca_results$data
  var_explained <- pca_results$variance

  if (!"host_type" %in% names(pca_data) && !is.null(host_classification)) {
    pca_data <- merge(pca_data,
      host_classification[, c("isolate", "host_type")],
      by = "isolate", all.x = TRUE
    )
  }

  p <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(
      ggplot2::aes(color = host_type, group = 1, text = isolate),
      size = 3, alpha = 0.7
    ) +
    ggplot2::scale_color_manual(
      values = c(monocot = "#2E86AB", dicot = "#A23B72", unknown = "#CCCCCC"),
      na.value = "#CCCCCC", name = "Host Type"
    ) +
    ggplot2::labs(
      title = "PCA of Codon Usage Patterns (RSCU)",
      x = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
      y = sprintf("PC2 (%.1f%% variance)", var_explained[2])
    ) +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.position = "right")

  if (interactive && requireNamespace("plotly", quietly = TRUE)) {
    p <- plotly::ggplotly(p, tooltip = c("text", "color"))
  }
  p
}

#' Scree plot of variance explained by each principal component
#'
#' @param pca_results List from [perform_rscu_pca_improved()].
#' @return ggplot object.
#' @export
plot_pca_scree <- function(pca_results) {
  var_data <- data.frame(
    PC         = paste0("PC", seq_along(pca_results$variance)),
    variance   = pca_results$variance,
    cumulative = cumsum(pca_results$variance)
  )
  var_data$PC <- factor(var_data$PC, levels = var_data$PC)
  n <- min(10, nrow(var_data))

  ggplot2::ggplot(var_data[seq_len(n), ], ggplot2::aes(x = PC, y = variance)) +
    ggplot2::geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.7) +
    ggplot2::geom_line(ggplot2::aes(y = cumulative, group = 1),
      color = "#D62828", linewidth = 1
    ) +
    ggplot2::geom_point(ggplot2::aes(y = cumulative),
      color = "#D62828", size = 2
    ) +
    ggplot2::scale_y_continuous(
      name = "Variance Explained (%)",
      sec.axis = ggplot2::sec_axis(~., name = "Cumulative Variance (%)")
    ) +
    ggplot2::labs(title = "PCA Scree Plot", x = "Principal Component") +
    cowplot::theme_cowplot() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Interactive 3D PCA scatter of RSCU codon usage
#'
#' Renders the first three principal components as a rotatable plotly 3D
#' scatter, coloured by host type. Requires the PCA to have been run with at
#' least three components (`perform_rscu_pca_improved()` returns `PC1..PC3`).
#'
#' @param pca_results List from [perform_rscu_pca_improved()]; its `$data`
#'   must contain `PC1`, `PC2`, `PC3`, `isolate` and (ideally) `host_type`.
#' @param host_classification Optional host classification data frame. If the
#'   PCA data lacks `host_type`, it is merged in from here by `isolate`.
#' @return A plotly object (or, if plotly is unavailable, the 2D biplot).
#' @export
plot_pca_3d <- function(pca_results, host_classification = NULL) {
  pca_data <- pca_results$data
  var_exp <- pca_results$variance

  if (!"host_type" %in% names(pca_data) && !is.null(host_classification)) {
    pca_data <- merge(pca_data,
      host_classification[, c("isolate", "host_type")],
      by = "isolate", all.x = TRUE
    )
  }

  if (!all(c("PC1", "PC2", "PC3") %in% names(pca_data))) {
    stop("pca_results$data must contain PC1, PC2 and PC3 columns")
  }

  if (!requireNamespace("plotly", quietly = TRUE)) {
    warning("plotly not available; returning 2D biplot instead")
    return(plot_pca_biplot(pca_results, host_classification))
  }

  plotly::plot_ly(
    pca_data,
    x = ~PC1, y = ~PC2, z = ~PC3,
    color = ~host_type,
    colors = c(monocot = "#2E86AB", dicot = "#A23B72", unknown = "#CCCCCC"),
    text = ~isolate,
    type = "scatter3d", mode = "markers",
    marker = list(size = 5)
  ) |>
    plotly::layout(
      title = "3D PCA of Codon Usage (RSCU)",
      scene = list(
        xaxis = list(title = sprintf("PC1 (%.1f%%)", var_exp[1])),
        yaxis = list(title = sprintf("PC2 (%.1f%%)", var_exp[2])),
        zaxis = list(title = sprintf("PC3 (%.1f%%)", var_exp[3]))
      )
    )
}

#' Violin/box plot of Codon Adaptation Index by host type
#'
#' @param cai_results Data frame from [calculate_host_specific_cai()]
#'   (columns `isolate`, `CAI`).
#' @param host_classification Data frame with columns `isolate`, `host_type`.
#' @return ggplot object.
#' @export
plot_cai_by_host <- function(cai_results, host_classification) {
  plot_data <- merge(cai_results,
    host_classification[, c("isolate", "host_type")],
    by = "isolate"
  )
  plot_data <- plot_data[!is.na(plot_data$host_type), ]

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = host_type, y = CAI, fill = host_type)
  ) +
    ggplot2::geom_violin(alpha = 0.7) +
    ggplot2::geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    ggplot2::scale_fill_manual(
      values = c(monocot = "#2E86AB", dicot = "#A23B72"), guide = "none"
    ) +
    ggplot2::labs(
      title = "Codon Adaptation Index (CAI) by Host Type",
      subtitle = "Higher CAI = better adaptation to host codon usage",
      x = "Host Type", y = "CAI"
    ) +
    cowplot::theme_cowplot()

  if (length(unique(plot_data$host_type)) == 2) {
    test_result <- stats::wilcox.test(CAI ~ host_type, data = plot_data)
    p <- p + ggplot2::annotate(
      "text",
      x = 1.5, y = max(plot_data$CAI) * 0.95,
      label = sprintf("Wilcoxon p = %.3f", test_result$p.value), size = 3
    )
  }
  p
}
