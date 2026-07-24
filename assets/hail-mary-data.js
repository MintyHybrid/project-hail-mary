/* ===========================================================================
   hail-mary-data.js — single source of truth for the site's navigation layer.
   Consumed by network-nav.js (landing graph, package graph) and workbench.js
   (per-chapter reasoning stepper).

   ---------------------------------------------------------------------------
   
   CHAPTER_STEPS below was seeded from the Claude Designer prototype and has
   since been AUTHOR-REVIEWED against the chapter text (workbench.js sets
   STEPS_REVIEWED = true, which drops the per-panel "draft summary" notice).
   The steps paraphrase rather than quote, so keep them in step with the
   chapters when the prose changes: the empirical figures must match index.qmd
   and the open questions (Ch.2 GARD screen, Ch.6 HGT origin, Ch.8/9 model
   results) must stay visibly open, per this book's hedging convention.
   Everything else here (chapters, edges, packages) is derived from the repo.
   ---------------------------------------------------------------------------

   Editing notes:
   - CHAPTERS.href is relative to the site root; network-nav.js rewrites it
     for pages rendered in notebooks/ subdirectories.
   - Each step is [headline, detail]. Plain mode shows the headline alone;
     expert mode shows the headline plus the detail.
   =========================================================================== */

(function (root) {
  "use strict";

  // --- chapters -----------------------------------------------------------
  // role: "origin" = where the whole thing started; "dest" = the synthesis.
  var CHAPTERS = [
    {
      id: "1", n: "01", title: "The DNAzyme-like motif discovery",
      short: "Discovery", href: "notebooks/01_discovery.html", role: "origin",
      expert: "A 10-23 DNAzyme catalytic core, BLASTed against the viral nucleotide collection, matches 16 potyvirus isolates across 6 species — and the 9+9 nt flanks turn up in eukaryotic genomes.",
      plain: "A hunch-driven database search turns up a strange shared snippet hiding inside several plant viruses."
    },
    {
      id: "2", n: "02", title: "Phylogenetics",
      short: "Phylogeny", href: "notebooks/02_phylogenetics.html",
      expert: "MAFFT alignment and codon-aware ML inference over the CI region, with codon-partitioned dN/dS across Potyviridae genera.",
      plain: "Building the virus family tree to see exactly where the snippet sits."
    },
    {
      id: "3", n: "03", title: "Codon usage & the DNAzyme locus",
      short: "Codon usage", href: "notebooks/03_codon_dnazyme.html",
      expert: "RSCU, ENC, GC3 and dinucleotide O/E around the DNAzyme locus, with a codon co-occurrence network and host-specific CAI.",
      plain: "Reading the fine print of the genetic code right around the snippet."
    },
    {
      id: "4", n: "04", title: "Protein motifs",
      short: "ELM motifs", href: "notebooks/04_protein_motifs.html",
      expert: "gget/ELM short-linear-motif scan of the translated CI region across Potyviridae genera, mapped back to genome coordinates.",
      plain: "Finding short functional 'address labels' hidden in the protein sequence."
    },
    {
      id: "5", n: "05", title: "The CI protein",
      short: "CI protein", href: "notebooks/05_potyvirus_motifs.html",
      expert: "CI structure prediction plus a curated motif catalogue and CI / HC-Pro / VPg–host interaction network.",
      plain: "Meet the protein — its 3D shape and the partners it talks to."
    },
    {
      id: "6", n: "06", title: "Horizontal gene transfer",
      short: "HGT", href: "notebooks/06_hgt.html",
      expert: "Phylogenetic incongruence testing and the lineage-specific Fragaria insert fused to the host TCP7 locus.",
      plain: "Did a virus sequence really end up inside a plant genome? Chasing the evidence."
    },
    {
      id: "7", n: "07", title: "Ecological networks",
      short: "Networks", href: "notebooks/07_networks.html",
      expert: "GBIF occurrence and GLOBI interaction records assembled into multilayer host–vector networks; the BioNexus app.",
      plain: "Mapping who-interacts-with-whom to see whether ecology can explain the pattern."
    },
    {
      id: "8", n: "08", title: "Immune simulation",
      short: "Immune model", href: "notebooks/08_immu_simmu.html",
      expert: "An interactive PTI/ETI trafficking simulation over AP-4 / ESCRT / autophagy receptor routes, perturbable live.",
      plain: "A living model of how a plant cell sorts its cargo while under attack."
    },
    {
      id: "9", n: "09", title: "AP-4 trafficking & the tyrosine motif",
      short: "Hypothesis", href: "notebooks/09_ap4_trafficking.html", role: "dest",
      expert: "Synthesis: the CI's tyrosine-based (YXXΦ) and acidic-dileucine sorting motifs as a route into host AP-mediated trafficking.",
      plain: "The big idea it all points to — hijacked cellular 'address labels'."
    },
    {
      id: "A", n: "A", title: "Data & code",
      short: "Data & code", href: "data-code.html",
      expert: "Accessions, the potytools package, and the targets pipeline — every result reproducible end to end.",
      plain: "Every dataset, script and citation, so you can reproduce it yourself."
    }
  ];

  // --- reasoning chain ----------------------------------------------------
  // Each edge is "this line of reasoning led to that one", not a citation.
  var EDGES = [
    ["1", "2"], ["1", "3"], ["1", "6"], ["1", "7"],
    ["2", "3"], ["2", "4"], ["4", "5"], ["5", "8"],
    ["5", "9"], ["8", "9"], ["6", "7"], ["7", "A"], ["3", "A"]
  ];

  var STEP_LABELS = ["Hypothesis", "Method", "Result", "Interpretation"];

  // --- per-chapter reasoning steps ---------------------------------------
  // ⚠ SEED TEXT — see the review warning at the top of this file.
  var CHAPTER_STEPS = {
    "1": [
      ["Could a catalytic-looking core sequence be conserved because it does something?",
       "A fragment matching the 10-23 deoxyribozyme catalytic core turned up while scanning viral sequences — specific enough to be worth chasing rather than dismissing as noise."],
      ["BLAST the core against the viral nucleotide collection; pull every hit; check the flanking sequence.",
       "BLASTn against NCBI, retaining hits above a strict identity/e-value threshold, then searching eukaryotic genomes for perfect matches to the 9 + 9 nt flanks of each isolate."],
      ["16 potyvirus isolates across 6 species carry the core, and their flanks match eukaryotic genomic loci.",
       "Not a single hit — a small, recurring family of matches spanning both viral and host sequence space, with 9 of the 16 isolates identical by sequence."],
      ["A shared motif spanning kingdoms is either deep homology, convergence, or transfer — worth chasing.",
       "The rest of the book exists to work out which of those it is."]
    ],
    "2": [
      ["If this locus is functionally constrained, it should leave phylogenetic and selection signal.",
       "Building a real tree lets us ask whether the motif region evolves like the rest of the genome, or differently."],
      ["MAFFT alignment of the CI region across Potyviridae; codon-aware ML tree; codon-partitioned dN/dS.",
       "A codon-aware phylogenetics pipeline, benchmarked against genome-wide topology."],
      ["Codon-partitioned dN/dS indicates purifying selection at the protein level, family-wide across 13 genera.",
       "The locus is constrained as protein, which is a different claim from the nucleotide motif itself being selected."],
      ["The signal needs a second lens — codon usage, not just tree shape.",
       "A GARD recombination screen has not yet been run, so recombination remains an unexcluded alternative explanation."]
    ],
    "3": [
      ["If the motif sits under distinct selective or mutational pressure, codon statistics should show it.",
       "Local codon usage tracks selection and genome context in ways whole-tree topology can miss."],
      ["RSCU, ENC, GC3 and dinucleotide observed/expected across the locus, plus a codon co-occurrence network.",
       "Position-resolved composition statistics computed per isolate and compared against genome background, with host-specific CAI reference tables."],
      ["The DNAzyme-like locus sits in a compositional outlier window relative to flanking coding sequence.",
       "A distinct local signature, consistent with either recent acquisition or strong local constraint — the two are not separated by composition alone."],
      ["Composition can't say why on its own — the next question is what the region actually encodes.",
       "Which moves the analysis from nucleotides to the protein the locus sits inside."]
    ],
    "4": [
      ["If this locus does something, it may carry a recognisable short linear motif, not just unusual base composition.",
       "Many trafficking and regulatory signals are encoded as short linear motifs (SLiMs), largely invisible to alignment alone."],
      ["gget/ELM query of the translated CI region across Potyviridae genera, with a known-domain overlay.",
       "Programmatic ELM database search, with hits mapped back onto genome and protein coordinates."],
      ["A recurring candidate sorting motif is present across genera, overlapping the CI helicase.",
       "Consistent presence across a virus family is a reasonable prior against the hit being spurious."],
      ["A motif implies a partner — so the CI protein itself becomes the next object of study.",
       "Motif presence alone says nothing about accessibility; that needs structure."]
    ],
    "5": [
      ["If the CI carries a real sorting motif, its structure and interactome should be consistent with that role.",
       "Structure and interaction context test whether the motif is surface-exposed and physically usable."],
      ["CI structure prediction plus a curated CI / HC-Pro / VPg–host interaction network.",
       "Predicted model cross-checked against published potyviral protein interaction studies."],
      ["The candidate tyrosine-based motif sits in an accessible position, alongside an acidic-dileucine AP-2 motif.",
       "Physically plausible: the motif is somewhere a sorting adaptor could actually reach it."],
      ["That raises a harder question — could the sequence, or its origin, be shared with the host?",
       "Which is the horizontal-transfer question, approached in Ch. 6."]
    ],
    "6": [
      ["Could the shared flank sequence reflect an actual transfer event rather than coincidence?",
       "Convergence is an unlikely way to produce near-identical multi-kb stretches, and transfer is at least testable."],
      ["Incongruence tests between host taxonomy and the viral-flank locus, plus targeted genome search.",
       "Comparing expected species topology against the gene tree for the candidate region, isolate by isolate."],
      ["A lineage-specific insert sits fused into the Fragaria (strawberry) TCP7 locus.",
       "A concrete, addressable integration event — though whether its origin is HGT specifically remains open."],
      ["If a viral sequence can land in a host genome, the ecological contact that allowed it matters.",
       "Which turns the question from sequence to ecology."]
    ],
    "7": [
      ["Transfer implies sustained ecological contact between virus, vector and host — is that visible in interaction data?",
       "Public occurrence and interaction data should reflect repeated contact opportunities if they exist."],
      ["GBIF occurrence and GLOBI interaction records assembled into multilayer host–vector networks (BioNexus).",
       "Layered network construction across geography, host range and known vector relationships."],
      ["Contact appears sporadic rather than continuous across the relevant host/vector neighbourhood.",
       "Ecological proximity is consistent with, but well short of proof of, a transfer route."],
      ["With a route at least plausible, the next question is what the motif does inside an infected cell.",
       "Moving from population-scale ecology to sub-cellular mechanism."]
    ],
    "8": [
      ["If the motif is a real sorting signal, perturbing host trafficking in a model should produce infection-like effects.",
       "A simulation lets the mechanism be probed without waiting on wet-lab timelines."],
      ["An interactive kinetic model of AP-4 / ESCRT / autophagy receptor trafficking, perturbable live.",
       "Rate relationships for receptor synthesis, sorting, recycling and degradation, exposed as a live simulation."],
      ["Blocking the AP-4 arm reproduces a receptor mis-sorting pattern.",
       "The simulated phenotype resembles what a hijacked sorting signal would be expected to cause — within the model's assumptions."],
      ["The mechanism now has a name: a virus-side tyrosine motif competing for host AP-4 recognition.",
       "Still a model result, not an experimental one."]
    ],
    "9": [
      ["The CI motif may act as an AP-4-competent tyrosine/dileucine sorting signal that reroutes host trafficking.",
       "Every prior chapter contributes a piece of evidence toward this single mechanistic hypothesis."],
      ["Synthesis across phylogenetics, motif detection, structure, HGT and simulation.",
       "No new data — this chapter integrates Ch. 1–8 into one explicitly testable model."],
      ["A parsimonious account connects a conserved motif, a real insert, and a reproducible modelled phenotype.",
       "The pieces fit without needing a separate explanation for each observation, which is suggestive rather than conclusive."],
      ["Next: wet-lab validation of AP-4 competition, and testing across additional Potyviridae hosts.",
       "Much of this chapter is predictive and labelled as such; it is a hypothesis, not a result."]
    ]
  };

  // --- package network ----------------------------------------------------
  // Derived from actual library()/:: usage across notebooks/, scripts/ and
  // shiny/, plus potytools' DESCRIPTION. Keep in sync when dependencies change.
  var PACKAGES = [
    { id: "rentrez",       lang: "R",      topic: "sequence retrieval", chapters: ["1"],           desc: "Programmatic NCBI E-utilities queries for isolate retrieval." },
    { id: "Biostrings",    lang: "R",      topic: "sequence retrieval", chapters: ["1", "3", "4"], desc: "Core containers and pattern matching for DNA/protein sequences." },
    { id: "BSgenome",      lang: "R",      topic: "sequence retrieval", chapters: ["6"],           desc: "Reference genome access for host flank searches." },
    { id: "plyranges",     lang: "R",      topic: "sequence retrieval", chapters: ["6"],           desc: "Tidy genomic range arithmetic over GenomicRanges." },

    { id: "msa",           lang: "R",      topic: "alignment",          chapters: ["1", "2"],      desc: "Multiple sequence alignment (ClustalW/MUSCLE/MAFFT wrappers)." },
    { id: "DECIPHER",      lang: "R",      topic: "alignment",          chapters: ["2"],           desc: "Alignment and sequence set curation for large collections." },

    { id: "ape",           lang: "R",      topic: "phylogenetics",      chapters: ["2", "6"],      desc: "Tree building, manipulation and comparative methods." },
    { id: "ggtree",        lang: "R",      topic: "phylogenetics",      chapters: ["2", "6"],      desc: "Grammar-of-graphics tree visualisation and annotation." },
    { id: "phangorn",      lang: "R",      topic: "phylogenetics",      chapters: ["2"],           desc: "Maximum likelihood and parsimony phylogenetic inference." },
    { id: "phytools",      lang: "R",      topic: "phylogenetics",      chapters: ["2"],           desc: "Comparative phylogenetics and ancestral state methods." },
    { id: "rotl",          lang: "R",      topic: "phylogenetics",      chapters: ["7"],           desc: "Open Tree of Life synthetic tree queries for host taxa." },

    { id: "potytools",     lang: "R",      topic: "codon usage",        chapters: ["1", "3", "6"], desc: "This project's own package: parsing, codon bias, host classification, motif flanks." },
    { id: "seqinr",        lang: "R",      topic: "codon usage",        chapters: ["3"],           desc: "Sequence I/O and codon/GC composition utilities." },

    { id: "gget",          lang: "Python", topic: "motif / structure",  chapters: ["4"],           desc: "Programmatic ELM linear-motif and UniProt queries." },
    { id: "pandas",        lang: "Python", topic: "motif / structure",  chapters: ["4"],           desc: "Tabular wrangling of ELM hit tables before they re-enter R." },
    { id: "ggseqlogo",     lang: "R",      topic: "motif / structure",  chapters: ["4", "6"],      desc: "Sequence logos for motif and flank conservation." },
    { id: "r3dmol",        lang: "R",      topic: "motif / structure",  chapters: ["5"],           desc: "Interactive 3D protein structure rendering." },

    { id: "igraph",        lang: "R",      topic: "networks",           chapters: ["5", "7", "A"], desc: "Graph construction, layout and network metrics." },
    { id: "visNetwork",    lang: "R",      topic: "networks",           chapters: ["5", "7"],      desc: "Interactive network widgets embedded in the chapters." },
    { id: "ggraph",        lang: "R",      topic: "networks",           chapters: ["7"],           desc: "Static publication-quality network layouts." },

    { id: "metacoder",     lang: "R",      topic: "ecology",            chapters: ["7"],           desc: "Taxonomic hierarchy visualisation for host ranges." },
    { id: "rphylopic",     lang: "R",      topic: "ecology",            chapters: ["7"],           desc: "Organism silhouettes for host/vector figures." },
    { id: "shiny",         lang: "R",      topic: "ecology",            chapters: ["7", "8"],      desc: "Reactive apps behind the BioNexus network explorer." },
    { id: "leaflet",       lang: "R",      topic: "ecology",            chapters: ["7"],           desc: "Geographic occurrence maps in the BioNexus app." },

    { id: "ggplot2",       lang: "R",      topic: "visualisation",      chapters: ["3", "6", "7"], desc: "The grammar of graphics underlying every static figure." },
    { id: "plotly",        lang: "R",      topic: "visualisation",      chapters: ["3", "8"],      desc: "Interactive versions of selected figures." },
    { id: "DT",            lang: "R",      topic: "visualisation",      chapters: ["5"],           desc: "Sortable, searchable data tables in the rendered book." },

    { id: "targets",       lang: "R",      topic: "reproducibility",    chapters: ["A"],           desc: "Pipeline orchestration so every expensive step rebuilds from raw data." },
    { id: "renv",          lang: "R",      topic: "reproducibility",    chapters: ["A"],           desc: "Locked per-project package library for exact reproducibility." },
    { id: "quarto",        lang: "R",      topic: "reproducibility",    chapters: ["A"],           desc: "Renders this book, with per-chapter freeze caching." }
  ];

  root.HailMary = {
    CHAPTERS: CHAPTERS,
    EDGES: EDGES,
    STEP_LABELS: STEP_LABELS,
    CHAPTER_STEPS: CHAPTER_STEPS,
    PACKAGES: PACKAGES
  };
})(window);
