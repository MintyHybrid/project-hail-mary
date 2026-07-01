#' fetch custom sequences from rentrez
#' 
#' @description fetch several accessions from rentrez using start end coordinated provided in a dataframe
#' @param df dataframe or tibble with column names accession, start, end, strand
#' @param output_file path/to/outputfile.fasta

fetch_custom_sequences <- function(df, output_file = "sequences.fasta") {
  # Load table (expects columns: accession, start, end, strand)
  # df <- read.csv(csv_file, stringsAsFactors = FALSE)
  
  # Clear existing file
  if (file.exists(output_file)) file.remove(output_file)
  
  for (i in 1:nrow(df)) {
    acc    <- df$seqnames[i]
    start  <- df$start[i]
    end    <- df$end[i]
    s_char <- df$strand[i]
    
    # Map strand: '+' -> 1 (plus), '-' -> 2 (minus)
    ncbi_strand <- ifelse(s_char == "+", 1, 2)
    
    
    tryCatch({
      # 1. Get Organism Name
      summ <- entrez_summary(db="nuccore", id=acc)
      org_name <- summ$title
      
      # 2. Fetch Sequence Fragment
      raw_fasta <- entrez_fetch(db="nuccore", id=acc, 
                                rettype="fasta", 
                                seq_start=start, 
                                seq_stop=end, 
                                strand=ncbi_strand)
      
      # 3. Format Custom Header: "> organism, accession.strand:start-end"
      # Remove NCBI's default header (first line)
      seq_lines <- strsplit(raw_fasta, "\n")[[1]]
      sequence_body <- paste(seq_lines[-1], collapse = "\n")
      
      # custom_header <- sprintf("> %s, %s.%s:%s-%s", 
      custom_header <- sprintf("> %s, %s|%s:%s-%s",                 # swappped . for | to avoid confusion with accession versioning
                               org_name, acc, s_char, start, end)
      
      # 4. Write directly to file (append mode)
      cat(custom_header, "\n", sequence_body, "\n", 
          file = output_file, append = TRUE)
      
      message(paste("Fetched:", acc))
      
      # Rate limit: ~3 requests per second to avoid NCBI blocks
      Sys.sleep(0.4)
      
    }, error = function(e) {
      message(paste("Error with", acc, ":", e$message))
    })
    
  }
}
