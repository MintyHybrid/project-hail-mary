test_that("build_concat produces correct positions", {
  frags <- c(A = "MDFDY", B = "LKPTG")
  result <- build_concat(frags)
  expect_named(result, c("concat", "map"))
  expect_equal(result$map$isolate, c("A", "B"))
  # Start of A is 1; end = nchar("MDFDY") = 5
  expect_equal(result$map$start[1], 1L)
  expect_equal(result$map$end[1], 5L)
  # B starts after linker
  linker_len <- nchar(potytools:::LINKER)
  expect_equal(result$map$start[2], 5L + linker_len + 1L)
})
