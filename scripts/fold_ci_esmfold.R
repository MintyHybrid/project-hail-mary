#!/usr/bin/env Rscript
# Predict a 3D structure of the potyvirus CI helicase region for the book's
# interactive viewer. No experimental potyvirus CI structure exists in the PDB
# and AlphaFold DB has no model for the viral polyprotein accessions, so we fold
# a representative CI helicase domain de novo with the ESMFold API.
#
# Steps: take a well-covered isolate from the ICTV alignment, extract an
# ungapped CI-region domain (anchored on the conserved VASYN motif, widened to a
# foldable ~350 aa domain), submit to ESMFold, and cache:
#   data/CI_esmfold.pdb            (predicted structure)
#   data/CI_esmfold_motif.txt      (1-based residue range of the motif region)
# read by notebook 04b (renders offline from the cache).

suppressPackageStartupMessages({ library(Biostrings); library(crul) })
root    <- "C:/Users/chris/R_projects/project-hail-mary"
out_pdb <- file.path(root, "data", "CI_esmfold.pdb")
out_mot <- file.path(root, "data", "CI_esmfold_motif.txt")

aln  <- readAAStringSet(file.path(root, "data", "OPSR.Poty.Fig3.v16_align.txt"))
cons <- as.character(consensusString(aln, ambiguityMap = "X"))
m    <- Biostrings::matchPattern("VASYN", cons, max.mismatch = 1)
anchor_a <- start(m)[1]; anchor_b <- end(m)[1]

# ~350 aa CI domain window around the anchor (clamped to the alignment)
win_lo <- max(1, anchor_a - 150)
win_hi <- min(nchar(cons), anchor_b + 200)

# choose the isolate with the fewest gaps across the window (best-covered)
region_aln <- subseq(aln, win_lo, win_hi)
gapfrac <- vapply(as.character(region_aln),
                  function(s) mean(strsplit(s, "")[[1]] == "-"), numeric(1))
best <- names(region_aln)[which.min(gapfrac)]
domain <- gsub("-", "", as.character(region_aln[[best]]))
domain <- gsub("X", "", domain)
if (nchar(domain) > 400) domain <- substr(domain, 1, 400)   # ESMFold length cap
message(sprintf("Folding %s CI domain (%d aa) via ESMFold ...", best, nchar(domain)))

# motif region (±20 around VASYN) position within the extracted domain
motif_start <- max(1, (anchor_a - 20) - win_lo + 1)
motif_end   <- min(nchar(domain), (anchor_b + 20) - win_lo + 1)

cli  <- crul::HttpClient$new(url = "https://api.esmatlas.com")
resp <- cli$post(path = "foldSequence/v1/pdb/", body = domain,
                 headers = list("Content-Type" = "text/plain"))
pdb  <- resp$parse("UTF-8")
if (resp$status_code != 200 || !grepl("^ATOM|^HEADER|MODEL", pdb)) {
  stop("ESMFold request failed (HTTP ", resp$status_code, ").")
}
writeLines(pdb, out_pdb)
writeLines(c(paste("isolate", best),
             paste("domain_len", nchar(domain)),
             paste("motif_start", motif_start),
             paste("motif_end", motif_end)), out_mot)
cat("Wrote", out_pdb, "(", nchar(domain), "aa );",
    "motif residues", motif_start, "-", motif_end, "\n")
