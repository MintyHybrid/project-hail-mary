#!/usr/bin/env Rscript
# Retrieve the homologous ~3 kb strawberry locus from RELATED Fragaria genomes so
# the presence/absence of the long (~3 kb) viral-homologous fragment can be
# tested against the host phylogeny (single ancient insertion vs. repeated
# horizontal transfer, or lineage-specific loss).
#
# strawberryseqs.fasta holds only fragment-POSITIVE Fragaria loci (~3 kb each).
# Here we BLAST a representative ~3 kb positive locus against Fragaria genomes on
# NCBI, then for every other genome measure how much of the 3 kb query aligns
# (query coverage). A genome carrying the full fragment covers ~all of the query;
# one where the fragment is absent aligns only over the flanking host sequence
# (low/partial coverage, with a gap over the insert). We fetch the spanning
# homologous locus from each hit for downstream alignment.
#
# Network + time heavy (remote BLAST ~1-3 min): run once; outputs cached to
#   data/strawberry_homologs.fasta
#   data/strawberry_homolog_table.tsv
# and read by notebook 05 (which renders offline from the cache).

suppressPackageStartupMessages({
  library(Biostrings); library(crul); library(rentrez)
  library(dplyr); library(stringr); library(xml2)
})

root    <- "C:/Users/chris/R_projects/project-hail-mary"
out_fa  <- file.path(root, "data", "strawberry_homologs.fasta")
out_tsv <- file.path(root, "data", "strawberry_homolog_table.tsv")

# ---- 1. Representative ~3 kb query (a fragment-positive Fragaria locus) -------
straw     <- readDNAStringSet(file.path(root, "data", "strawberryseqs.fasta"))
# version-stripped accessions we already have (XML Hit_accession has no version)
have_acc  <- str_remove(str_extract(names(straw), "[A-Z]{2}[0-9]{6}\\.[0-9]"), "\\.[0-9]+$")
query_seq <- as.character(straw[[1]])
qlen      <- nchar(query_seq)
message("Query length: ", qlen, " bp")

# ---- 2. Submit remote BLAST (blastn) across Fragaria + close relatives -------
# Broadened beyond Fragaria to nearby Rosoideae genera so that assemblies
# lacking the ~3 kb fragment are also recovered (needed to score absence).
cli <- crul::HttpClient$new(url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi")
entrez_q <- paste("Fragaria[Organism] OR Potentilla[Organism] OR",
                  "Rosa[Organism] OR Duchesnea[Organism] OR Rubus[Organism]")
message("Submitting remote BLAST (Fragaria + relatives) ...")
put <- cli$post(query = list(
  CMD = "Put", PROGRAM = "blastn", MEGABLAST = "on", DATABASE = "core_nt",
  QUERY = query_seq, ENTREZ_QUERY = entrez_q, HITLIST_SIZE = 300))
put_txt <- put$parse("UTF-8")
rid <- sub(".*RID = (\\S+).*", "\\1",
           regmatches(put_txt, regexpr("RID = \\S+", put_txt)))
stopifnot(nzchar(rid)); message("  RID = ", rid)

# ---- 3. Poll until ready -----------------------------------------------------
repeat {
  Sys.sleep(20)
  st <- cli$get(query = list(CMD = "Get", RID = rid,
                             FORMAT_OBJECT = "SearchInfo"))$parse("UTF-8")
  if (grepl("Status=READY", st)) break
  if (grepl("Status=FAILED|Status=UNKNOWN", st)) stop("BLAST job failed: ", rid)
  message("  ... waiting")
}

# ---- 4. Retrieve HSPs as XML (the URL-API Tabular format is unreliable) ------
# Fetch FORMAT_TYPE=XML and parse <Hit>/<Hsp>; retry until real hits appear.
doc <- NULL
for (attempt in 1:8) {
  Sys.sleep(10)
  xml_txt <- cli$get(query = list(CMD = "Get", RID = rid,
                                  FORMAT_TYPE = "XML"))$parse("UTF-8")
  if (grepl("<Hit>", xml_txt, fixed = TRUE)) {
    doc <- xml2::read_xml(xml_txt); break
  }
  message("  results not formatted yet (attempt ", attempt, ") ...")
}
if (is.null(doc)) stop("No BLAST hits returned for RID ", rid,
                       " after retries (query may have no Fragaria hits).")

hit_nodes <- xml2::xml_find_all(doc, ".//Hit")
gv <- function(node, xp) xml2::xml_text(xml2::xml_find_first(node, xp))
hsp <- dplyr::bind_rows(lapply(hit_nodes, function(h) {
  acc <- gv(h, "./Hit_accession")
  hsps <- xml2::xml_find_all(h, ".//Hsp")
  dplyr::bind_rows(lapply(hsps, function(x) data.frame(
    subject  = acc,
    pident   = 100 * as.numeric(gv(x, "./Hsp_identity")) /
                     as.numeric(gv(x, "./Hsp_align-len")),
    length   = as.integer(gv(x, "./Hsp_align-len")),
    qstart   = as.integer(gv(x, "./Hsp_query-from")),
    qend     = as.integer(gv(x, "./Hsp_query-to")),
    sstart   = as.integer(gv(x, "./Hsp_hit-from")),
    send     = as.integer(gv(x, "./Hsp_hit-to")),
    bitscore = as.numeric(gv(x, "./Hsp_bit-score")),
    stringsAsFactors = FALSE)))
}))
message("  parsed ", nrow(hsp), " HSPs from ", length(hit_nodes), " hits")

# ---- 5. Per-subject query coverage (union of query intervals) ----------------
qcov <- function(qs, qe) {
  iv <- cbind(pmin(qs, qe), pmax(qs, qe))
  iv <- iv[order(iv[, 1]), , drop = FALSE]
  tot <- 0; cur_s <- iv[1, 1]; cur_e <- iv[1, 2]
  for (i in seq_len(nrow(iv))[-1]) {
    if (iv[i, 1] <= cur_e + 1) { cur_e <- max(cur_e, iv[i, 2]) }
    else { tot <- tot + (cur_e - cur_s + 1); cur_s <- iv[i, 1]; cur_e <- iv[i, 2] }
  }
  (tot + (cur_e - cur_s + 1)) / qlen
}

subj <- hsp %>%
  dplyr::mutate(subject = str_remove(subject, "\\.[0-9]+$")) %>%
  dplyr::filter(!subject %in% have_acc) %>%
  dplyr::group_by(subject) %>%
  dplyr::summarise(
    query_coverage = qcov(qstart, qend),
    span_start = min(pmin(sstart, send)),
    span_end   = max(pmax(sstart, send)),
    mean_pident = round(mean(pident), 1),
    minus_strand = stats::median(as.integer(send < sstart)) >= 0.5,
    .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(query_coverage))

# classify by how much of the 3 kb query is present
subj <- subj %>% dplyr::mutate(fragment = dplyr::case_when(
  query_coverage >= 0.70 ~ "present",
  query_coverage >= 0.25 ~ "partial",
  TRUE                   ~ "flank_only"))
message(sprintf("  %d related Fragaria loci: %d present / %d partial / %d flank-only",
                nrow(subj), sum(subj$fragment == "present"),
                sum(subj$fragment == "partial"), sum(subj$fragment == "flank_only")))

# ---- 6. Fetch each spanning locus (~3 kb) via efetch -------------------------
pad <- 300L; recs <- character(0); meta <- list()
for (i in seq_len(nrow(subj))) {
  acc <- subj$subject[i]
  s0  <- max(1, subj$span_start[i] - pad); s1 <- subj$span_end[i] + pad
  fa  <- tryCatch(rentrez::entrez_fetch(
           db = "nuccore", id = acc, rettype = "fasta",
           seq_start = s0, seq_stop = s1,
           strand = if (subj$minus_strand[i]) 2 else 1),
         error = function(e) NULL)
  if (is.null(fa)) next
  body <- paste(strsplit(fa, "\n")[[1]][-1], collapse = "")
  summ <- tryCatch(rentrez::entrez_summary("nuccore", acc), error = function(e) NULL)
  recs[acc] <- body
  meta[[acc]] <- data.frame(
    accession = acc, organism = if (!is.null(summ)) summ$title else acc,
    length = nchar(body), query_coverage = round(subj$query_coverage[i], 3),
    fragment = subj$fragment[i], mean_pident = subj$mean_pident[i],
    source = "blast_homolog", stringsAsFactors = FALSE)
  Sys.sleep(0.4)
}

homologs <- DNAStringSet(recs); names(homologs) <- names(recs)
writeXStringSet(homologs, out_fa)

meta_df <- dplyr::bind_rows(meta) %>%
  dplyr::mutate(species = stringr::str_extract(organism,
                  "Fragaria x ananassa|Fragaria [a-z]+"))
readr::write_tsv(meta_df, out_tsv)

cat("Wrote", length(homologs), "homologous loci to", out_fa, "\n")
print(dplyr::count(meta_df, species, fragment))
