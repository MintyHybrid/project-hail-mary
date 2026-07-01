test_that("extract_codons returns correct codons", {
  seq <- "ATGGCTTAG"
  codons <- extract_codons(seq)
  expect_equal(codons, c("ATG", "GCT", "TAG"))
})

test_that("extract_codons ignores trailing incomplete codon", {
  seq <- "ATGGCT"
  expect_length(extract_codons(seq), 2L)
})

test_that("extract_codons handles Ns by dropping them", {
  seq <- "ATGNCTTAG"
  codons <- extract_codons(seq)
  # NCT is dropped, ATG and TAG remain
  expect_true(all(nchar(codons) == 3))
  expect_false(any(grepl("N", codons)))
})

test_that("calculate_gc_content returns values between 0 and 1", {
  seq <- "ATGCATGC"
  gc <- calculate_gc_content(seq)
  expect_true(all(gc >= 0 & gc <= 1))
  expect_named(gc, c("GC_overall", "GC1.pos1", "GC2.pos2", "GC3.pos3"),
    ignore.order = TRUE)
})

test_that("calculate_gc_content: 100% GC sequence", {
  seq <- "GCGGCGGCG"
  gc <- calculate_gc_content(seq)
  expect_equal(unname(gc["GC_overall"]), 1.0)
})

test_that("calculate_enc returns value in [20, 61]", {
  seq <- paste(rep("ATG", 30), collapse = "")  # single codon repeated
  enc <- calculate_enc(seq)
  expect_true(is.numeric(enc))
  if (!is.na(enc)) {
    expect_true(enc >= 20 && enc <= 61)
  }
})

test_that("calculate_enc returns NA for very short sequences", {
  enc <- calculate_enc("ATG")
  expect_true(is.na(enc))
})

test_that("compare_codon_usage returns a data frame with required columns", {
  motif <- c("ATGGCTTAG", "ATGGCATAG")
  ref   <- c("ATGGCTTAG", "GCGGCTATG")
  result <- compare_codon_usage(motif, ref)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("codon", "amino_acid", "fold_change", "p_value", "p_adjusted")
  %in% names(result)))
})

test_that("GENETIC_CODE covers 64 codons including stops", {
  expect_length(GENETIC_CODE, 64L)
  expect_true("ATG" %in% names(GENETIC_CODE))
  expect_equal(GENETIC_CODE[["ATG"]], "M")
})
