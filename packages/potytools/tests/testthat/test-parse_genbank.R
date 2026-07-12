test_that("parse_genbank_file returns a list with expected fields", {
  # Write a minimal GenBank record to a temp file
  gb_text <- c(
    "LOCUS       TEST_SEQ     21 bp    DNA     linear   VRL 01-JAN-2026",
    "DEFINITION  Test sequence for unit test.",
    "ACCESSION   TEST001",
    "VERSION     TEST001.1",
    "SOURCE      Potyvirus testus",
    "  ORGANISM  Potyvirus testus",
    "FEATURES             Location/Qualifiers",
    "     CDS             1..21",
    '                     /gene="test_gene"',
    '                     /translation="MGS"',
    "ORIGIN",
    "        1 atggcttcgt ag",
    "//"
  )
  tmp <- withr::local_tempfile(fileext = ".gb")
  writeLines(gb_text, tmp)

  records <- parse_genbank_file(tmp)
  expect_type(records, "list")
  expect_length(records, 1L)

  rec <- records[[1]]
  expect_equal(rec$accession, "TEST001")
  expect_equal(rec$mol_type, "DNA")
  expect_true(length(rec$features) > 0)
})

test_that("load_genbank_folder warns on empty folder", {
  tmp_dir <- withr::local_tempdir()
  expect_warning(load_genbank_folder(tmp_dir))
})

test_that("extract_accessions finds GenBank-style accessions", {
  expect_equal(extract_accessions("Isolate_AB011819 stuff"), "AB011819")
  expect_equal(extract_accessions("NC_001445.1 ref genome"), "NC_001445.1")
  expect_equal(extract_accessions("no accession here"), character(0))
})

test_that("has_genbank_id detects presence of an accession", {
  expect_true(has_genbank_id("Isolate_AB011819"))
  expect_false(has_genbank_id("plain_isolate_name"))
})

test_that("match_cds_to_reference matches by protein_id and extracts the CDS", {
  gb_text <- c(
    "LOCUS       TEST_SEQ     12 bp    DNA     linear   VRL 01-JAN-2026",
    "ACCESSION   TEST001",
    "FEATURES             Location/Qualifiers",
    "     CDS             1..12",
    '                     /protein_id="ABC12345.1"',
    '                     /product="polyprotein"',
    '                     /translation="MAS"',
    "ORIGIN",
    "        1 atggcttcgt ag",
    "//"
  )
  tmp <- withr::local_tempfile(fileext = ".gb")
  writeLines(gb_text, tmp)
  rec <- parse_genbank_file(tmp)[[1]]

  cds <- match_cds_to_reference(rec, name = "Isolate_ABC12345.1")
  expect_s4_class(cds, "DNAString")
  expect_equal(as.character(cds), "ATGGCTTCGTAG")
})

test_that("match_cds_to_reference falls back to product-name match", {
  gb_text <- c(
    "LOCUS       TEST_SEQ     12 bp    DNA     linear   VRL 01-JAN-2026",
    "ACCESSION   TEST001",
    "FEATURES             Location/Qualifiers",
    "     CDS             1..12",
    '                     /product="cylindrical inclusion protein"',
    '                     /translation="MAS"',
    "ORIGIN",
    "        1 atggcttcgt ag",
    "//"
  )
  tmp <- withr::local_tempfile(fileext = ".gb")
  writeLines(gb_text, tmp)
  rec <- parse_genbank_file(tmp)[[1]]

  cds <- match_cds_to_reference(rec, name = "cylindrical inclusion protein")
  expect_equal(as.character(cds), "ATGGCTTCGTAG")
})

test_that("match_cds_to_reference falls back to translation similarity", {
  gb_text <- c(
    "LOCUS       TEST_SEQ     12 bp    DNA     linear   VRL 01-JAN-2026",
    "ACCESSION   TEST001",
    "FEATURES             Location/Qualifiers",
    "     CDS             1..12",
    '                     /translation="MAS"',
    "ORIGIN",
    "        1 atggcttcgt ag",
    "//"
  )
  tmp <- withr::local_tempfile(fileext = ".gb")
  writeLines(gb_text, tmp)
  rec <- parse_genbank_file(tmp)[[1]]

  cds <- match_cds_to_reference(rec, name = "no_match_here", ref_protein = "MAS")
  expect_equal(as.character(cds), "ATGGCTTCGTAG")

  no_hit <- match_cds_to_reference(rec, name = "no_match_here", ref_protein = "WWWWWW")
  expect_null(no_hit)
})

test_that("match_cds_to_reference reverse-complements minus-strand CDS", {
  gb_text <- c(
    "LOCUS       TEST_SEQ     12 bp    DNA     linear   VRL 01-JAN-2026",
    "ACCESSION   TEST001",
    "FEATURES             Location/Qualifiers",
    "     CDS             complement(1..12)",
    '                     /protein_id="XYZ99999.1"',
    "ORIGIN",
    "        1 ctacgaagcc at",
    "//"
  )
  tmp <- withr::local_tempfile(fileext = ".gb")
  writeLines(gb_text, tmp)
  rec <- parse_genbank_file(tmp)[[1]]

  cds <- match_cds_to_reference(rec, name = "Isolate_XYZ99999.1")
  expect_equal(as.character(cds), as.character(Biostrings::reverseComplement(Biostrings::DNAString("CTACGAAGCCAT"))))
})

test_that("match_cds_to_reference returns NULL with no name/ref_protein match", {
  gb_text <- c(
    "LOCUS       TEST_SEQ     12 bp    DNA     linear   VRL 01-JAN-2026",
    "ACCESSION   TEST001",
    "FEATURES             Location/Qualifiers",
    "     CDS             1..12",
    '                     /product="polyprotein"',
    "ORIGIN",
    "        1 atggcttcgt ag",
    "//"
  )
  tmp <- withr::local_tempfile(fileext = ".gb")
  writeLines(gb_text, tmp)
  rec <- parse_genbank_file(tmp)[[1]]

  # empty/degenerate name must not spuriously match every CDS
  expect_null(match_cds_to_reference(rec, name = ""))
  expect_null(match_cds_to_reference(rec, name = "@@@"))
})

test_that("features_as_df returns a data frame", {
  gb_text <- c(
    "LOCUS       T1  9 bp  DNA linear VRL 01-JAN-2026",
    "ACCESSION   T1",
    "FEATURES             Location/Qualifiers",
    "     gene            1..9",
    '                     /gene="x"',
    "ORIGIN",
    "        1 atggcttag",
    "//"
  )
  tmp <- withr::local_tempfile(fileext = ".gb")
  writeLines(gb_text, tmp)
  records <- list(T1 = parse_genbank_file(tmp)[[1]])
  df <- features_as_df(records)
  expect_s3_class(df, "data.frame")
  expect_true("type" %in% names(df))
})
