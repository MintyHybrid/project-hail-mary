test_that("classify_host_type returns monocot for lily", {
  expect_equal(classify_host_type("lily"), "monocot")
})

test_that("classify_host_type returns dicot for tobacco", {
  expect_equal(classify_host_type("tobacco"), "dicot")
})

test_that("classify_host_type returns NA for unknown host", {
  expect_true(is.na(classify_host_type("xenomorph")))
})

test_that("detect_host_from_name extracts lily from underscore-separated name", {
  host <- detect_host_from_name("FMV_lily_CHN01")
  expect_false(is.na(host))
  expect_true(grepl("lily", host, ignore.case = TRUE))
})

test_that("create_host_classification_table returns correct dimensions", {
  names <- c("PVY_tobacco_NL1", "WMV_cucumber_USA", "Unknown_seq")
  result <- create_host_classification_table(names)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3L)
  expect_true(all(c("isolate", "detected_host", "host_type") %in% names(result)))
})

test_that("get_codon_table_for_host_type returns 64-element vector", {
  mono <- get_codon_table_for_host_type("monocot")
  di   <- get_codon_table_for_host_type("dicot")
  expect_length(mono, 64L)
  expect_length(di,   64L)
  expect_equal(sum(mono), 1.0, tolerance = 1e-6)
  expect_equal(sum(di),   1.0, tolerance = 1e-6)
})
