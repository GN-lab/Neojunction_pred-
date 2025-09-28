#!/usr/bin/env Rscript
# Title: "Step 02: Protein-Coding Transcript Extraction ── SSNIP (Hartwig Edition) ───────────────────────────
# September 24, 2025 | Gaurav Raichand | The Insitute of Cancer Research
# Edits: Added flexible mitochondrial chromosome filtering (e.g., 'MT' or 'chrM'),
#        dynamic output filename with current date, logging for filtered counts,
#        and basic GTF validation for Hartwig compatibility.

suppressPackageStartupMessages({
  library(tidyverse)    # For data manipulation
  library(data.table)   # For fread
  library(stringr)      # For string extraction
})

# Step 0: Load GTF from config paths
input_dir <- Sys.getenv("INPUT_DIR")
setwd(input_dir)

filename_gtf <- Sys.getenv("GTF_FILE")
message("[INFO] Reading GTF file: ", basename(filename_gtf))

data_gtf <- as_tibble(data.table::fread(filename_gtf, sep = "\t", header = FALSE, skip = "#"))  # Skip comment lines if present

# Basic validation: Ensure at least 9 columns (standard GTF format)
if (ncol(data_gtf) < 9) {
  stop("[ERROR] GTF file must have at least 9 columns. Check your Hartwig-compatible GTF.")
}

# Step 1: Rename columns
colnames(data_gtf) <- paste0("X", 1:9)

# Step 2: Filter protein-coding transcripts (exclude mitochondrial variants)
# Flexible MT filter: Handles 'MT', 'chrM', 'chrMT' common in GRCh37/38 GTFs
prot_coding <- data_gtf %>% 
  filter(!X1 %in% c("MT", "chrM", "chrMT")) %>%  # New: Flexible mitochondrial exclusion
  filter(X3 == "transcript") %>%                 # Select transcript features
  filter(str_detect(X9, "protein_coding"))       # Select protein-coding (check if your GTF uses this tag)

# Validation: Warn if no protein-coding entries found
if (nrow(prot_coding) == 0) {
  warning("[WARN] No protein-coding transcripts found. Verify 'protein_coding' tag in GTF attributes.")
}

message("[INFO] Filtered to ", nrow(prot_coding), " protein-coding transcripts (excluding mitochondrial).")

# Step 3: Extract enst, symbol, ensg from X9
prot_coding <- prot_coding %>% 
  mutate(
    enst = str_extract(X9, 'transcript_id "[^"]+"') %>% str_remove_all('transcript_id "|"'),
    symbol = str_extract(X9, 'gene_name "[^"]+"') %>% str_remove_all('gene_name "|"'),
    ensg = str_extract(X9, 'gene_id "[^"]+"') %>% str_remove_all('gene_id "|"')
  )

# Step 4: Output to config-defined directory
output_dir <- Sys.getenv("OUTPUT_DIR")
setwd(output_dir)

# Dynamic filename with current date (for version tracking in Hartwig runs)
current_date <- format(Sys.Date(), "%Y%m%d")  # New: Uses today's date (e.g., 20250924)
filename_output <- paste0("GTF_Protein_Coding_Genes_", current_date, ".txt")
write.table(prot_coding, filename_output, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

message("[INFO] Written: ", filename_output)
