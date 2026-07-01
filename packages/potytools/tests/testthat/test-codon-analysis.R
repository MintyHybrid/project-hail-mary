test_that("extract_codons returns expected codons", {
  seq <- "ATGGACTTCGAC"
  codons <- extract_codons(seq)
  expect_equal(codons, c("ATG", "GAC", "TTC", "GAC"))
})

test_that("extract_codons handles frame offset", {
  seq <- "XXATGGACTTC"  # frame = 2 means skip 2 chars
  codons <- extract_codons(seq, frame = 2)
  expect_equal(codons, c("ATG", "GAC", "TTC"))
})

test_that("extract_codons drops codons with N", {
  seq <- "ATGNNN GAC"
  codons <- extract_codons(gsub(" ", "", seq))
  expect_false("NNN" %in% codons)
  expect_true("ATG" %in% codons)
})

test_that("calculate_rscu_correct returns expected structure", {
  seqs <- c("ATGATGATG")  # all Met; no synonymous codons → RSCU = 1
  result <- calculate_rscu_correct(seqs)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("codon", "amino_acid", "RSCU") %in% names(result)))
  expect_false(any(result$amino_acid == "*"))
})

test_that("calculate_gc_fast returns named vector with 4 values", {
  res <- calculate_gc_fast("ATGATGATG")
  expect_named(res, c("GC_overall", "GC1", "GC2", "GC3"))
  expect_true(all(res >= 0 & res <= 1))
})

test_that("calculate_enc_fast returns a numeric in (2, 61]", {
  seq <- paste(rep("ATGCGA", 10), collapse = "")
  enc <- calculate_enc_fast(seq)
  expect_true(is.numeric(enc))
  expect_true(enc > 2 && enc <= 61)
})

test_that("run_complete_analysis_optimized returns expected list keys", {
  seqs <- Biostrings::DNAStringSet(c(s1 = paste(rep("ATG", 20), collapse = ""),
                                     s2 = paste(rep("GAC", 20), collapse = "")))
  ref  <- Biostrings::DNAStringSet(c(r1 = paste(rep("TTC", 20), collapse = "")))
  result <- run_complete_analysis_optimized(seqs, ref, output_prefix = tempfile())
  expect_named(result,
    c("rscu_comparison", "usage_comparison", "motif_rscu", "ref_rscu",
      "enc_motif", "enc_ref", "gc_motif", "gc_ref"))
})
