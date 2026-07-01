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
