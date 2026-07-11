# _targets.R — Project Hail Mary pipeline
#
# Run with: targets::tar_make()
# Visualise: targets::tar_visnetwork()
# Load a result: targets::tar_read(target_name)
#
# External tool requirements (must be on PATH):
#   mafft   — multiple sequence alignment
#   iqtree2 — maximum likelihood phylogenetics
#   hmmscan — protein domain annotation (HMMER suite)

library(targets)
library(tarchetypes)

# ── Package loading ────────────────────────────────────────────────────────────
tar_option_set(
  packages = c(
    "Biostrings", "dplyr", "stringr", "purrr", "rentrez",
    "ape", "seqinr", "potytools"
  ),
  format = "rds"
)

# ── Helper: wrap external CLI tool ────────────────────────────────────────────
run_mafft <- function(input_fasta, output_fasta, threads = 4L) {
  cmd <- sprintf("mafft --auto --thread %d %s > %s",
                 threads, shQuote(input_fasta), shQuote(output_fasta))
  stopifnot(system(cmd) == 0L)
  output_fasta
}

run_iqtree <- function(alignment_fasta, prefix, model = "GTR+G") {
  cmd <- sprintf("iqtree2 -s %s -m %s -pre %s -nt AUTO -bb 1000",
                 shQuote(alignment_fasta), model, shQuote(prefix))
  stopifnot(system(cmd) == 0L)
  paste0(prefix, ".treefile")
}

# ── Targets ───────────────────────────────────────────────────────────────────
list(

  # ── 0. Raw data inputs ──────────────────────────────────────────────────────
  tar_target(ictv_alignment_file,
             "data/OPSR.Poty.Fig3.v16_align.txt",
             format = "file"),

  tar_target(genbank_folder,
             "data/genbank/",
             format = "file"),

  # ── 1. Parse GenBank files ──────────────────────────────────────────────────
  tar_target(gb_records,
             load_genbank_folder(genbank_folder)),

  tar_target(all_proteins,
             get_all_proteins(gb_records)),

  tar_target(feature_table,
             features_as_df(gb_records)),

  # ── 2. Load ICTV alignment ──────────────────────────────────────────────────
  tar_target(ictv_alignment,
             Biostrings::readAAStringSet(ictv_alignment_file)),

  # ── 3. Extract CI-region flanks (± 20 aa) ───────────────────────────────────
  tar_target(ci_flanks,
             extract_motif_flanks(ictv_alignment,
                                  motif_pattern = "GSGKS",
                                  flank_aa = 20L)),

  # ── 4. NCBI sequence retrieval ──────────────────────────────────────────────
  tar_target(coord_table,
             read.csv("data/ci_coords.csv")),

  tar_target(fetched_seqs_file, {
    fetch_custom_sequences(coord_table,
                           output_file = "data/ci_cds_fetched.fasta")
    "data/ci_cds_fetched.fasta"
  }, format = "file"),

  tar_target(cds_seqs,
             Biostrings::readDNAStringSet(fetched_seqs_file)),

  # ── 5. Codon usage analysis ─────────────────────────────────────────────────
  tar_target(rscu_results,     calculate_rscu(cds_seqs)),
  tar_target(enc_results,      sapply(as.character(cds_seqs), calculate_enc)),
  tar_target(gc_results,       lapply(as.character(cds_seqs), calculate_gc_content)),

  # ── 6. Host classification ──────────────────────────────────────────────────
  tar_target(host_table,
             create_host_classification_table(names(cds_seqs))),

  tar_target(host_codon_tables,
             create_host_specific_codon_tables(host_table)),

  tar_target(cai_results,
             calculate_host_specific_cai(cds_seqs, host_codon_tables)),

  # Note: ELM motif search moved to the gget-based pipeline
  # (scripts/prep_elm_regions.R -> scripts/gget_elm_run.py); see notebook 04.

  # ── 8. Multiple sequence alignment (external: MAFFT) ────────────────────────
  tar_target(aligned_fasta, {
    run_mafft("data/ci_cds_fetched.fasta", "data/ci_aligned.fasta")
  }, format = "file"),

  # ── 9. Phylogenetic inference (external: IQ-TREE 2) ─────────────────────────
  tar_target(treefile, {
    run_iqtree("data/ci_aligned.fasta", prefix = "data/phylo/ci_tree")
  }, format = "file"),

  tar_target(phylo_tree,
             ape::read.tree(treefile)),

  # ── 10. Motif variant analysis ──────────────────────────────────────────────
  tar_target(motif_seqs,
             Biostrings::readDNAStringSet("data/motif_18nt.fasta")),

  tar_target(motif_variants,
             analyze_motif_variants(as.character(motif_seqs))),

  tar_target(codon_comparison,
             compare_codon_usage(motif_seqs, cds_seqs)),

  # ── 11. Render book ─────────────────────────────────────────────────────────
  tar_quarto(manuscript, path = ".")
)
