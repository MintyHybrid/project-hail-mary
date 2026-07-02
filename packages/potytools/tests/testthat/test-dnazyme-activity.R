test_that("dnazyme_core_table has 15 positions and marks G14 essential", {
  tb <- dnazyme_core_table()
  expect_equal(nrow(tb), 15)
  expect_equal(paste(tb$base, collapse = ""), CANONICAL_1023_CORE)
  expect_equal(tb$tier[14], "essential")
  expect_true(all(tb$weight >= 0 & tb$weight <= 1))
})

test_that("canonical core scores as fully active", {
  res <- score_dnazyme_activity(c(canon = CANONICAL_1023_CORE))
  expect_equal(res$n_substitutions, 0L)
  expect_equal(res$weighted_identity, 1)
  expect_equal(res$activity_score, 1)
  expect_true(res$predicted_active)
})

test_that("substitution at the essential general base (G14) predicts inactive", {
  # canonical GGCTAGCTACAACGA -> change position 14 (G) to A
  mut <- CANONICAL_1023_CORE
  substr(mut, 14, 14) <- "A"
  res <- score_dnazyme_activity(c(g14 = mut))
  expect_equal(res$activity_score, 0)
  expect_false(res$predicted_active)
  expect_equal(res$critical_substitutions, "14")
})

test_that("a tolerated (modulatory) substitution keeps some predicted activity", {
  mut <- CANONICAL_1023_CORE
  substr(mut, 4, 4) <- "A"   # T4 -> A, weight 0.80
  res <- score_dnazyme_activity(c(t4 = mut))
  expect_gt(res$activity_score, 0)      # 1 - 0.80 = 0.20 retained
  expect_lt(res$weighted_identity, 1)
})

test_that("core window is located within a longer sequence", {
  longer <- paste0("TT", CANONICAL_1023_CORE, "AA")
  res <- score_dnazyme_activity(c(x = longer))
  expect_equal(res$core_window, CANONICAL_1023_CORE)
  expect_equal(res$n_substitutions, 0L)
})
