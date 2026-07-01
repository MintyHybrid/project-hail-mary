# Memoised wrappers around rentrez so repeated calls for the same accession
# within a session don't hit NCBI again.
.entrez_summary_memo <- memoise::memoise(rentrez::entrez_summary)
.entrez_fetch_memo   <- memoise::memoise(rentrez::entrez_fetch)

#' Fetch sequence fragments from NCBI for a table of accession/coordinate rows.
#'
#' @description Fetch several accessions from NCBI via rentrez using
#'   start/end coordinates. Results for each (accession, start, end, strand)
#'   combination are cached in-session via memoise so re-running the function
#'   does not issue duplicate network requests.
#'
#' @param df Data frame or tibble with columns: seqnames, start, end, strand
#' @param output_file Path to the output FASTA file (appended if it exists)
#' @export
fetch_custom_sequences <- function(df, output_file = "sequences.fasta") {
  if (file.exists(output_file)) file.remove(output_file)

  for (i in seq_len(nrow(df))) {
    acc         <- df$seqnames[i]
    start       <- df$start[i]
    end         <- df$end[i]
    s_char      <- df$strand[i]
    ncbi_strand <- ifelse(s_char == "+", 1, 2)

    tryCatch({
      summ     <- .entrez_summary_memo(db = "nuccore", id = acc)
      org_name <- summ$title

      raw_fasta <- .entrez_fetch_memo(
        db        = "nuccore",
        id        = acc,
        rettype   = "fasta",
        seq_start = start,
        seq_stop  = end,
        strand    = ncbi_strand
      )

      seq_lines     <- strsplit(raw_fasta, "\n")[[1]]
      sequence_body <- paste(seq_lines[-1], collapse = "\n")
      custom_header <- sprintf("> %s, %s|%s:%s-%s", org_name, acc, s_char, start, end)

      cat(custom_header, "\n", sequence_body, "\n", file = output_file, append = TRUE)
      message("Fetched: ", acc)
      Sys.sleep(0.4)   # NCBI rate limit: ~3 req/s

    }, error = function(e) {
      message("Error with ", acc, ": ", e$message)
    })
  }
}
