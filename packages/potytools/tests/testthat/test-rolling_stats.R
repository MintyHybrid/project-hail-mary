test_that("rolling_alignment_stats returns one row per window with expected cols", {
  set.seed(1)
  aln <- matrix(sample(c("A", "C", "G", "T"), 6 * 90, replace = TRUE),
                nrow = 6)
  rownames(aln) <- paste0("seq", 1:6)

  ras <- rolling_alignment_stats(aln, window = 30, step = 30)

  expect_s3_class(ras, "data.frame")
  expect_equal(nrow(ras), 3) # 90 sites / 30 step, non-overlapping
  expect_true(all(c("window", "start", "end", "GC", "AT", "GC3",
                    "CpG_OE", "gap_fraction", "entropy",
                    "GC_pos1", "GC_pos2", "GC_pos3") %in% names(ras)))
  expect_true(all(ras$GC >= 0 & ras$GC <= 1))
})

test_that("rolling_alignment_stats enforces codon-multiple window/step", {
  aln <- matrix("A", nrow = 2, ncol = 60)
  expect_error(rolling_alignment_stats(aln, window = 31, step = 30),
               "multiples of 3")
})

test_that("rolling_alignment_stats rejects non-matrix input", {
  expect_error(rolling_alignment_stats(list()), "character matrix")
})

test_that("per-sequence variant yields sequence x window rows", {
  set.seed(2)
  aln <- matrix(sample(c("A", "C", "G", "T"), 4 * 60, replace = TRUE),
                nrow = 4)
  rownames(aln) <- paste0("s", 1:4)

  ras_ps <- rolling_alignment_stats_per_sequence(aln, window = 30, step = 30)

  expect_equal(nrow(ras_ps), 4 * 2) # 4 seqs x 2 windows
  expect_true("sequence" %in% names(ras_ps))
  expect_setequal(unique(ras_ps$sequence), rownames(aln))
})

test_that("remove_gaps drops gap characters before metrics", {
  aln <- matrix(c(rep("-", 15), rep("G", 15)), nrow = 2, byrow = TRUE)
  # window covers all 15 sites of each row
  with_gaps <- rolling_alignment_stats(aln, window = 15, step = 15,
                                       codon_aligned = FALSE)
  no_gaps <- rolling_alignment_stats(aln, window = 15, step = 15,
                                     codon_aligned = FALSE, remove_gaps = TRUE)
  expect_gt(with_gaps$gap_fraction[1], 0)
  # GC computed over non-gap only should be >= GC over all (gaps count as 0)
  expect_gte(no_gaps$GC[1], with_gaps$GC[1])
})
