#!/usr/bin/env Rscript
# Broader strawberry BLAST (W3bc) — reproduce the two-flank hit distribution.
#
# The earlier fetch (fetch_strawberry_homologs.R) restricted hits to a few
# Rosoideae genera and returned only 3 RefSeq TCP7 mRNAs, which align to just
# ONE flank of the ~3 kb positive locus. A genomic BLAST across the whole
# Rosaceae recovers many subjects whose HSPs cluster at BOTH ends of the query
# with a central gap — i.e. the putative insert is the middle ~1.4 kb, flanked
# by conserved host sequence on both sides (matching the user's NCBI BLAST).
#
# Saves the PER-HSP table (one row per HSP, with query coordinates) so notebook
# 05 can render a faithful hit-distribution + per-position coverage figure
# offline:
#   data/strawberry_broad_blast_hsps.tsv
# Network + time heavy (remote BLAST ~2-4 min): run once; output cached.

suppressPackageStartupMessages({
  library(Biostrings); library(crul); library(dplyr); library(stringr)
  library(xml2)
})

root    <- "C:/Users/chris/R_projects/project-hail-mary"
out_tsv <- file.path(root, "data", "strawberry_broad_blast_hsps.tsv")

straw     <- readDNAStringSet(file.path(root, "data", "strawberryseqs.fasta"))
query_seq <- as.character(straw[[1]])
qlen      <- nchar(query_seq)
message("Query length: ", qlen, " bp")

cli <- crul::HttpClient$new(url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi")
# Broaden to the whole Rosaceae family (Malus, Prunus, Pyrus, Fragaria, Rosa ...)
entrez_q <- "Rosaceae[Organism]"
message("Submitting broad remote BLAST (Rosaceae) ...")
put <- cli$post(query = list(
  CMD = "Put", PROGRAM = "blastn", MEGABLAST = "on", DATABASE = "core_nt",
  QUERY = query_seq, ENTREZ_QUERY = entrez_q, HITLIST_SIZE = 250))
put_txt <- put$parse("UTF-8")
rid <- sub(".*RID = (\\S+).*", "\\1",
           regmatches(put_txt, regexpr("RID = \\S+", put_txt)))
stopifnot(nzchar(rid)); message("  RID = ", rid)

repeat {
  Sys.sleep(20)
  st <- cli$get(query = list(CMD = "Get", RID = rid,
                             FORMAT_OBJECT = "SearchInfo"))$parse("UTF-8")
  if (grepl("Status=READY", st)) break
  if (grepl("Status=FAILED|Status=UNKNOWN", st)) stop("BLAST job failed: ", rid)
  message("  ... waiting")
}

doc <- NULL
for (attempt in 1:10) {
  Sys.sleep(10)
  xml_txt <- cli$get(query = list(CMD = "Get", RID = rid,
                                  FORMAT_TYPE = "XML"))$parse("UTF-8")
  if (grepl("<Hit>", xml_txt, fixed = TRUE)) { doc <- read_xml(xml_txt); break }
  message("  results not formatted yet (attempt ", attempt, ") ...")
}
if (is.null(doc)) stop("No BLAST hits returned for RID ", rid)

hit_nodes <- xml_find_all(doc, ".//Hit")
gv <- function(node, xp) xml_text(xml_find_first(node, xp))
hsp <- bind_rows(lapply(hit_nodes, function(h) {
  acc   <- gv(h, "./Hit_accession")
  title <- gv(h, "./Hit_def")
  hsps  <- xml_find_all(h, ".//Hsp")
  bind_rows(lapply(hsps, function(x) {
    qs <- as.integer(gv(x, "./Hsp_query-from"))
    qe <- as.integer(gv(x, "./Hsp_query-to"))
    data.frame(
      subject  = acc, title = title,
      qstart   = min(qs, qe), qend = max(qs, qe),
      pident   = 100 * as.numeric(gv(x, "./Hsp_identity")) /
                       as.numeric(gv(x, "./Hsp_align-len")),
      length   = as.integer(gv(x, "./Hsp_align-len")),
      bitscore = as.numeric(gv(x, "./Hsp_bit-score")),
      evalue   = as.numeric(gv(x, "./Hsp_evalue")),
      stringsAsFactors = FALSE)
  }))
}))
message("  parsed ", nrow(hsp), " HSPs from ", length(hit_nodes), " subjects")

hsp$organism <- str_extract(hsp$title,
  "^[A-Z][a-z]+ [a-z]+(?: (?:subsp\\.|var\\.) [a-z]+)?")
attr(hsp, "qlen") <- qlen
readr::write_tsv(hsp, out_tsv)
message("Wrote ", nrow(hsp), " HSPs (qlen=", qlen, ") to ", out_tsv)
