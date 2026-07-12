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

test_that("codon_align_from_protein projects gaps and codons correctly", {
  prot_aln <- Biostrings::AAStringSet(c(a = "M-K", b = "MAK"))
  cds <- Biostrings::DNAStringSet(c(a = "ATGAAA", b = "ATGGCAAAA"))
  codon_aln <- codon_align_from_protein(prot_aln, cds)
  expect_equal(as.character(codon_aln[["a"]]), "ATG---AAA")
  expect_equal(as.character(codon_aln[["b"]]), "ATGGCAAAA")
})

test_that("codon_align_from_protein pads with NNN when CDS runs out early", {
  prot_aln <- Biostrings::AAStringSet(c(a = "MKV"))
  cds <- Biostrings::DNAStringSet(c(a = "ATGAAA")) # only 2 codons for 3 residues
  codon_aln <- codon_align_from_protein(prot_aln, cds)
  expect_equal(as.character(codon_aln[["a"]]), "ATGAAANNN")
})

test_that("codon_align_from_protein only projects names present in both", {
  prot_aln <- Biostrings::AAStringSet(c(a = "MK", b = "MK"))
  cds <- Biostrings::DNAStringSet(c(a = "ATGAAA"))
  codon_aln <- codon_align_from_protein(prot_aln, cds)
  expect_equal(names(codon_aln), "a")
})

test_that("safe_translate trims to a whole number of codons and translates", {
  aa <- safe_translate(Biostrings::DNAString("ATGGCTTA"))  # 8 nt -> trims to 6
  expect_equal(as.character(aa), "MA")
})

test_that("safe_translate strips alignment gaps before translating", {
  aa <- safe_translate(Biostrings::DNAString("ATG---GCT"))
  expect_equal(as.character(aa), "MA")
})

test_that("calculate_gc_content returns values between 0 and 1", {
  seq <- "ATGCATGC"
  gc <- calculate_gc_content(seq)
  expect_true(all(gc >= 0 & gc <= 1))
  expect_named(gc, c("GC_overall", "GC1", "GC2", "GC3"),
    ignore.order = TRUE
  )
})

test_that("calculate_gc_content: 100% GC sequence", {
  seq <- "GCGGCGGCG"
  gc <- calculate_gc_content(seq)
  expect_equal(unname(gc["GC_overall"]), 1.0)
})

test_that("calculate_enc returns value in [20, 61] for a realistic CDS", {
  # A sequence spanning several amino acids with synonymous codon variation,
  # so Wright's Nc is well defined.
  seq <- paste(rep(paste0(
    "GCTGCCGCAGCG", # Ala x4 (4-fold family, varied)
    "CTTCTCCTGCTA", # Leu x4 (6-fold family, varied)
    "GACGATGAAGAG", # Asp/Glu (2-fold families)
    "AAAAAGTTTTTC", # Lys/Phe (2-fold families)
    "TATTACATGTGG" #  Tyr (2-fold) + Met/Trp (single-codon)
  ), 4), collapse = "")
  enc <- calculate_enc(seq)
  expect_true(is.numeric(enc))
  expect_false(is.na(enc))
  expect_true(enc >= 20 && enc <= 61)
})

test_that("calculate_enc returns NA for degenerate (single amino-acid) input", {
  # All-Met sequence has no synonymous families, so Nc is undefined.
  expect_true(is.na(calculate_enc(paste(rep("ATG", 30), collapse = ""))))
})

test_that("calculate_enc returns NA for very short sequences", {
  enc <- calculate_enc("ATG")
  expect_true(is.na(enc))
})

test_that("compare_codon_usage returns a data frame with required columns", {
  motif <- c("ATGGCTTAG", "ATGGCATAG")
  ref <- c("ATGGCTTAG", "GCGGCTATG")
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

# --- merged from former test-codon-analysis.R: fast/optimized variants ---

test_that("extract_codons returns expected codons", {
  seq <- "ATGGACTTCGAC"
  codons <- extract_codons(seq)
  expect_equal(codons, c("ATG", "GAC", "TTC", "GAC"))
})

test_that("extract_codons handles frame offset", {
  seq <- "XXATGGACTTC" # frame = 2 means skip 2 chars
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
  seqs <- c("ATGATGATG") # all Met; no synonymous codons → RSCU = 1
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
  seqs <- Biostrings::DNAStringSet(c(
    s1 = paste(rep("ATG", 20), collapse = ""),
    s2 = paste(rep("GAC", 20), collapse = "")
  ))
  ref <- Biostrings::DNAStringSet(c(r1 = paste(rep("TTC", 20), collapse = "")))
  result <- run_complete_analysis_optimized(seqs, ref, output_prefix = tempfile())
  expect_named(
    result,
    c(
      "rscu_comparison", "usage_comparison", "motif_rscu", "ref_rscu",
      "enc_motif", "enc_ref", "gc_motif", "gc_ref"
    )
  )
})

test_that("calculate_gc_matrix_fast matches per-sequence calculate_gc_fast", {
  seqs <- Biostrings::DNAStringSet(c(
    a = "ATGGCTGACGGG", b = "GGGCCCTTTAAA", c = "ATATATCGCGCG"
  ))
  m <- calculate_gc_matrix_fast(seqs)
  expect_equal(dim(m), c(3L, 4L))
  expect_equal(colnames(m), c("GC_overall", "GC1", "GC2", "GC3"))
  # Row-by-row agreement with the scalar implementation
  for (i in seq_along(seqs)) {
    expect_equal(unname(m[i, ]),
      unname(calculate_gc_fast(as.character(seqs[i]))),
      tolerance = 1e-9)
  }
})

test_that("calculate_gc_matrix_fast handles a single sequence", {
  m <- calculate_gc_matrix_fast(Biostrings::DNAStringSet(c(x = "ATGCATGCATGC")))
  expect_equal(dim(m), c(1L, 4L))
  expect_true(all(m >= 0 & m <= 1))
})
