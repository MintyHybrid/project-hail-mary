# Sequence Logo Generation
# Creates visual representation of sequence conservation


# ==============================================================================
# SEQUENCE LOGO FUNCTIONS
# ==============================================================================

#' Create DNA sequence logo from aligned sequences
#'
#' @param sequences DNAStringSet or character vector (must be aligned)
#' @param title Plot title
#' @param color_scheme Color scheme ("auto", "nucleotide", "base_pairing")
#' @return ggplot object
#' @export
create_dna_logo <- function(sequences,
                            title = "DNA Sequence Logo",
                            color_scheme = "auto") {
  if (inherits(sequences, "DNAStringSet")) {
    sequences <- as.character(sequences)
  }

  # Convert to uppercase
  sequences <- toupper(sequences)

  # Check alignment
  seq_lengths <- nchar(sequences)
  if (length(unique(seq_lengths)) > 1) {
    warning("Sequences have different lengths - not properly aligned")
  }

  # Create logo
  p <- ggseqlogo(sequences,
    method = "bits",
    seq_type = "dna"
  ) +
    labs(
      title = title,
      x = "Position",
      y = "Information Content (bits)"
    ) +
    theme_logo() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    theme_manuscript_bg()

  if (color_scheme == "base_pairing") {
    # Custom colors: A-T (red/blue), G-C (yellow/green)
    p <- p + scale_fill_manual(values = c(
      A = "#E41A1C", T = "#377EB8",
      G = "#FFD700", C = "#4DAF4A"
    ))
  }

  return(p)
}

#' Create protein sequence logo from translated sequences
#'
#' @param sequences AAStringSet or character vector (amino acids)
#' @param title Plot title
#' @param color_scheme Color scheme ("chemistry", "hydrophobicity", "auto")
#' @return ggplot object
#' @export
create_protein_logo <- function(sequences,
                                title = "Protein Sequence Logo",
                                color_scheme = "chemistry") {
  if (inherits(sequences, "AAStringSet")) {
    sequences <- as.character(sequences)
  }

  # Convert to uppercase
  sequences <- toupper(sequences)

  # Create logo
  p <- ggseqlogo(sequences,
    method = "bits",
    seq_type = "aa",
    col_scheme = color_scheme
  ) +
    labs(
      title = title,
      x = "Position",
      y = "Information Content (bits)"
    ) +
    theme_logo() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    theme_manuscript_bg()

  return(p)
}

#' Create combined DNA and protein logo plot
#'
#' @param dna_sequences DNAStringSet
#' @param title Main title
#' @return Combined ggplot object
#' @export
create_combined_logo <- function(dna_sequences,
                                 title = "Motif Sequence Logos") {
  # DNA logo
  p_dna <- create_dna_logo(dna_sequences,
    title = "DNA Sequence",
    color_scheme = "auto"
  )

  # Translate sequences
  aa_sequences <- translate(dna_sequences)

  # Protein logo
  p_protein <- create_protein_logo(aa_sequences,
    title = "Translated Protein",
    color_scheme = "chemistry"
  )

  # Combine plots
  combined <- plot_grid(
    p_dna,
    p_protein,
    ncol = 1,
    rel_heights = c(1, 1),
    labels = c("A", "B"),
    label_size = 14
  )

  # Add main title
  title_plot <- ggdraw() +
    draw_label(title, fontface = "bold", size = 16)

  final_plot <- plot_grid(
    title_plot,
    combined,
    ncol = 1,
    rel_heights = c(0.1, 1)
  ) + theme_manuscript_bg()

  return(final_plot)
}

# ==============================================================================
# POSITION-SPECIFIC CONSERVATION
# ==============================================================================

#' Calculate Shannon entropy for each position
#'
#' @param sequences Character vector or DNAStringSet
#' @return Vector of entropy values (lower = more conserved)
#' @export
calculate_position_entropy <- function(sequences) {
  if (inherits(sequences, "DNAStringSet")) {
    sequences <- as.character(sequences)
  }

  sequences <- toupper(sequences)
  seq_length <- nchar(sequences[1])

  entropy_values <- numeric(seq_length)

  for (pos in 1:seq_length) {
    bases_at_pos <- substr(sequences, pos, pos)
    counts <- table(bases_at_pos)
    freqs <- counts / sum(counts)

    # Shannon entropy
    H <- -sum(freqs * log2(freqs + 1e-10))
    entropy_values[pos] <- H
  }

  return(entropy_values)
}

#' Plot conservation across sequence
#'
#' @param sequences DNAStringSet or character vector
#' @param title Plot title
#' @return ggplot object
#' @export
plot_conservation <- function(sequences, title = "Sequence Conservation") {
  entropy_values <- calculate_position_entropy(sequences)

  # Convert entropy to conservation score (2 - H)
  conservation <- 2 - entropy_values

  plot_data <- data.frame(
    position = 1:length(conservation),
    conservation = conservation,
    entropy = entropy_values
  )

  p <- ggplot(plot_data, aes(x = position, y = conservation)) +
    geom_line(color = "#2E86AB", size = 1.2) +
    geom_point(color = "#2E86AB", size = 2) +
    geom_area(alpha = 0.3, fill = "#2E86AB") +
    scale_y_continuous(
      limits = c(0, 2),
      breaks = seq(0, 2, 0.5),
      sec.axis = sec_axis(~ 2 - ., name = "Entropy (bits)")
    ) +
    labs(
      title = title,
      x = "Position",
      y = "Conservation (bits)"
    ) +
    theme_manuscript() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )

  return(p)
}

# ==============================================================================
# SEQUENCE LOGO SUMMARY STATISTICS
# ==============================================================================

#' Calculate summary statistics for sequence logo
#'
#' @param sequences DNAStringSet or character vector
#' @return Data frame with position-wise statistics
#' @export
logo_statistics <- function(sequences) {
  if (inherits(sequences, "DNAStringSet")) {
    sequences <- as.character(sequences)
  }

  sequences <- toupper(sequences)
  seq_length <- nchar(sequences[1])

  results <- data.frame(
    position = 1:seq_length,
    entropy = numeric(seq_length),
    conservation = numeric(seq_length),
    most_common_base = character(seq_length),
    frequency = numeric(seq_length),
    stringsAsFactors = FALSE
  )

  for (pos in 1:seq_length) {
    bases_at_pos <- substr(sequences, pos, pos)
    counts <- table(bases_at_pos)
    freqs <- counts / sum(counts)

    # Entropy
    H <- -sum(freqs * log2(freqs + 1e-10))
    results$entropy[pos] <- H
    results$conservation[pos] <- 2 - H

    # Most common
    most_common_idx <- which.max(counts)
    results$most_common_base[pos] <- names(counts)[most_common_idx]
    results$frequency[pos] <- max(freqs)
  }

  return(results)
}

# ==============================================================================
# EXPORT FUNCTIONS
# ==============================================================================

#' Save sequence logo to file
#'
#' @param sequences DNAStringSet
#' @param output_file Output file path
#' @param width Plot width
#' @param height Plot height
#' @param dpi Resolution
#' @export
save_combined_logo <- function(sequences,
                               output_file = "sequence_logo.png",
                               width = 12,
                               height = 8,
                               dpi = 300) {
  p <- create_combined_logo(sequences)

  ggsave(output_file, p, width = width, height = height, dpi = dpi)

  cat(sprintf("Saved sequence logo to: %s\n", output_file))
}

#' Create combined DNA and protein logos
#'
#' @param dna_sequences DNAStringSet
#' @param protein_sequences AAStringSet or character vector
#' @param title_prefix Title prefix string
#' @return Combined ggplot
#' @export
create_combined_logos <- function(dna_sequences, protein_sequences,
                                  title_prefix = "Motif") {
  p_dna <- create_dna_logo(dna_sequences,
    title = paste0(title_prefix, " - DNA Sequence Logo")
  )
  p_protein <- create_protein_logo(protein_sequences,
    title = paste0(title_prefix, " - Protein Sequence Logo")
  )
  cowplot::plot_grid(p_dna, p_protein, ncol = 1, rel_heights = c(1, 1)) +
    theme_manuscript_bg()
}

#' Host-specific sequence logos (monocot vs dicot)
#'
#' @param sequences DNAStringSet or AAStringSet
#' @param host_classification Data frame with columns isolate, host_type
#' @param seq_type "dna" or "protein"
#' @return Combined ggplot or NULL
#' @export
create_host_specific_logos <- function(sequences, host_classification,
                                       seq_type = "dna") {
  seqs_char <- if (inherits(sequences, c("DNAStringSet", "AAStringSet"))) {
    as.character(sequences)
  } else {
    sequences
  }
  mono_idx <- which(host_classification$host_type == "monocot")
  dict_idx <- which(host_classification$host_type == "dicot")
  if (length(mono_idx) == 0 || length(dict_idx) == 0) {
    warning("Need both monocot and dicot sequences for comparison")
    return(NULL)
  }
  if (seq_type == "dna") {
    p1 <- create_dna_logo(seqs_char[mono_idx], title = "Monocot-infecting Isolates (DNA)")
    p2 <- create_dna_logo(seqs_char[dict_idx], title = "Dicot-infecting Isolates (DNA)")
  } else {
    p1 <- create_protein_logo(seqs_char[mono_idx], title = "Monocot-infecting Isolates (Protein)")
    p2 <- create_protein_logo(seqs_char[dict_idx], title = "Dicot-infecting Isolates (Protein)")
  }
  cowplot::plot_grid(p1, p2, ncol = 1, labels = c("A", "B")) +
    theme_manuscript_bg()
}

#' Information content per position
#'
#' @param sequences Character vector of equal-length sequences
#' @param alphabet Valid character set
#' @return Numeric vector of information content values
#' @export
calculate_information_content <- function(sequences,
                                          alphabet = c("A", "C", "G", "T")) {
  if (inherits(sequences, c("DNAStringSet", "AAStringSet", "RNAStringSet"))) {
    sequences <- as.character(sequences)
  }
  sequences <- toupper(sequences)
  seq_length <- nchar(sequences[1])
  n_seqs <- length(sequences)
  max_entropy <- log2(length(alphabet))
  ic_values <- numeric(seq_length)
  for (pos in seq_len(seq_length)) {
    chars <- substr(sequences, pos, pos)
    counts <- table(factor(chars, levels = alphabet))
    freqs <- counts[counts > 0] / n_seqs
    entropy <- -sum(freqs * log2(freqs))
    ic_values[pos] <- max_entropy - entropy
  }
  ic_values
}

#' Plot information content along sequence
#'
#' @param sequences Character vector
#' @param title Plot title
#' @param seq_type "DNA" or "protein"
#' @return ggplot object
#' @export
plot_information_content <- function(sequences, title = "Information Content",
                                     seq_type = "DNA") {
  if (inherits(sequences, c("DNAStringSet", "AAStringSet", "RNAStringSet"))) {
    sequences <- as.character(sequences)
  }
  sequences <- toupper(sequences)
  alphabet <- if (seq_type == "DNA") {
    c("A", "C", "G", "T")
  } else {
    c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  }
  ic_values <- calculate_information_content(sequences, alphabet)
  plot_data <- data.frame(
    position = seq_along(ic_values),
    information_content = ic_values
  )
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = position, y = information_content)) +
    ggplot2::geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.7) +
    ggplot2::geom_line(color = "#D62828", linewidth = 1) +
    ggplot2::labs(title = title, x = "Position", y = "Information Content (bits)") +
    theme_manuscript()
  if (seq_type == "DNA" && nchar(sequences[1]) == 18) {
    p <- p + ggplot2::geom_vline(
      xintercept = seq(3.5, 18, by = 3),
      linetype = "dotted", alpha = 0.5
    )
  }
  p
}

#' Calculate PSSM (Position-Specific Scoring Matrix)
#'
#' @param sequences DNAStringSet, AAStringSet, or character vector
#' @param background_freq Optional named background frequency vector
#' @return Log2-odds matrix (letters x positions)
#' @export
calculate_pssm <- function(sequences, background_freq = NULL) {
  if (inherits(sequences, "DNAStringSet")) {
    sequences <- as.character(sequences)
    alphabet <- c("A", "C", "G", "T")
  } else if (inherits(sequences, "AAStringSet")) {
    sequences <- as.character(sequences)
    alphabet <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  } else {
    alphabet <- c("A", "C", "G", "T")
  }
  sequences <- toupper(sequences)
  seq_length <- nchar(sequences[1])
  n_seqs <- length(sequences)
  freq_matrix <- matrix(0,
    nrow = length(alphabet), ncol = seq_length,
    dimnames = list(alphabet, paste0("Pos", seq_len(seq_length)))
  )
  for (pos in seq_len(seq_length)) {
    chars <- substr(sequences, pos, pos)
    counts <- table(factor(chars, levels = alphabet))
    freq_matrix[, pos] <- counts / n_seqs
  }
  freq_matrix <- (freq_matrix + 0.01)
  freq_matrix <- sweep(freq_matrix, 2, colSums(freq_matrix), "/")
  if (is.null(background_freq)) {
    background_freq <- rep(1 / length(alphabet), length(alphabet))
  }
  log2(freq_matrix / background_freq)
}

#' Plot PSSM as heatmap
#'
#' @param pssm PSSM matrix (letters x positions)
#' @param title Plot title
#' @return ggplot object
#' @export
plot_pssm_heatmap <- function(pssm, title = "Position-Specific Scoring Matrix") {
  pssm_df <- as.data.frame(as.table(pssm))
  colnames(pssm_df) <- c("Letter", "Position", "Score")
  ggplot2::ggplot(pssm_df, ggplot2::aes(x = Position, y = Letter, fill = Score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", Score)), size = 2.5) +
    ggplot2::scale_fill_gradient2(
      low = "#2E86AB", mid = "white", high = "#D62828",
      midpoint = 0, name = "Log2\nOdds"
    ) +
    ggplot2::labs(title = title, x = "Position", y = "Letter") +
    theme_manuscript() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}

#' Comprehensive sequence logo analysis
#'
#' @param dna_sequences DNAStringSet
#' @param host_classification Data frame with columns isolate, host_type
#' @param output_prefix Prefix for saved PNG files
#' @return Named list of ggplot objects
#' @export
create_comprehensive_logo_analysis <- function(dna_sequences, host_classification,
                                               output_prefix = "logo_analysis") {
  protein_sequences <- Biostrings::translate(dna_sequences,
    no.init.codon = TRUE,
    if.fuzzy.codon = "solve"
  )
  p_dna <- create_dna_logo(dna_sequences, title = "DNA Sequence Logo (All Isolates)")
  p_protein <- create_protein_logo(protein_sequences, title = "Protein Sequence Logo (All Isolates)")
  p_combined <- create_combined_logos(dna_sequences, protein_sequences, title_prefix = "18nt Motif")
  p_host_dna <- create_host_specific_logos(dna_sequences, host_classification, "dna")
  p_host_protein <- create_host_specific_logos(protein_sequences, host_classification, "protein")
  p_ic_dna <- plot_information_content(as.character(dna_sequences), "DNA Information Content", "DNA")
  p_ic_protein <- plot_information_content(as.character(protein_sequences), "Protein Information Content", "protein")
  pssm_dna <- calculate_pssm(dna_sequences)
  pssm_protein <- calculate_pssm(protein_sequences)
  p_pssm_dna <- plot_pssm_heatmap(pssm_dna, "DNA Position-Specific Scoring Matrix")
  p_pssm_protein <- plot_pssm_heatmap(pssm_protein, "Protein Position-Specific Scoring Matrix")
  list(
    dna_logo = p_dna, protein_logo = p_protein, combined_logo = p_combined,
    host_dna_logos = p_host_dna, host_protein_logos = p_host_protein,
    ic_dna = p_ic_dna, ic_protein = p_ic_protein,
    pssm_dna = p_pssm_dna, pssm_protein = p_pssm_protein
  )
}
