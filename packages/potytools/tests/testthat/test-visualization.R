make_pca_stub <- function(with_pc3 = TRUE) {
  d <- data.frame(
    isolate   = paste0("iso", 1:6),
    PC1       = rnorm(6), PC2 = rnorm(6),
    host_type = rep(c("monocot", "dicot"), 3),
    stringsAsFactors = FALSE
  )
  if (with_pc3) d$PC3 <- rnorm(6)
  list(data = d, variance = c(40, 25, 15, 10, 6, 4))
}

test_that("plot_pca_3d returns a plotly object when plotly is available", {
  skip_if_not_installed("plotly")
  p <- plot_pca_3d(make_pca_stub())
  expect_s3_class(p, "plotly")
})

test_that("plot_pca_3d errors when PC3 is missing", {
  skip_if_not_installed("plotly")
  expect_error(plot_pca_3d(make_pca_stub(with_pc3 = FALSE)),
               "PC1, PC2 and PC3")
})

test_that("plot_pca_biplot returns a ggplot", {
  p <- plot_pca_biplot(make_pca_stub())
  expect_s3_class(p, "ggplot")
})
