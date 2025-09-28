#!/usr/bin/env Rscript
# Title: "Step 03: TPM Expression Filtering ── SSNIP (Hartwig Edition) ───────────────────────────
# September 24, 2025 | Gaurav Raichand | The Insitute of Cancer Research
# Edits: Removed log2 TPM conversion (Hartwig Salmon outputs are direct TPM), integrated config env vars,
#        dynamic filenames/dates, flexible full sample ID matching, generalized to 'all' group,
#        added logging/checks for mismatches and NA handling.


suppressPackageStartupMessages({
  library(tidyverse)    # For data manipulation
  library(data.table)   # For fread
  library(ggsci)        # For plotting (if needed later)
})


# Step 0: Load data from config paths
input_dir <- Sys.getenv("INPUT_DIR")
output_dir <- Sys.getenv("OUTPUT_DIR")


# Input filenames (use env vars or dynamic from prior steps)
tpm_file <- Sys.getenv("TPM_FILE")  # e.g., transcript_tpm_matrix.tsv from Salmon
purity_thresh <- as.numeric(Sys.getenv("MIN_PURITY", unset = "0.60"))
purity_file <- file.path(output_dir, sprintf("Patient_List_Post_TumorPurity_Filter_%.2f.txt", purity_thresh))  # From Step 01
gtf_coding_file <- file.path(output_dir, paste0("GTF_Protein_Coding_Genes_", format(Sys.Date(), "%Y%m%d"), ".txt"))  # From Step 02


message("[INFO] Loading TPM matrix: ", tpm_file)
dataframe_tpm <- as_tibble(data.table::fread(tpm_file))


message("[INFO] Loading purity-filtered samples: ", purity_file)
dataframe_purity <- as_tibble(data.table::fread(purity_file)) %>% rename(case = sample_id)  # Standardize to 'case'


message("[INFO] Loading protein-coding GTF extract: ", gtf_coding_file)
dataframe_gtf <- as_tibble(data.table::fread(gtf_coding_file))


# Validation: Check for sample mismatches (use full IDs, no truncation)
shared_samples <- intersect(colnames(dataframe_tpm)[-c(1:2)], dataframe_purity$case)  # Exclude tx/gene_id columns
if (length(shared_samples) == 0) {
  stop("[ERROR] No matching samples between TPM matrix and purity file. Check ID formats (e.g., Hartwig full strings).")
}
message("[INFO] Found ", length(shared_samples), " matching samples out of ", nrow(dataframe_purity), " purity-filtered.")


# Step 1: Preprocess TPM data
# Filter for ENST transcripts (assuming 'tx' column is transcript_id)
dataframe_tpm_edit <- dataframe_tpm %>%
  rename(enst = tx) %>%  # Rename to match expected 'enst'
  filter(str_starts(enst, "ENST"))


# Retain only shared samples (TPM values already in linear scale – no conversion needed)
dataframe_tpm_edit <- dataframe_tpm_edit %>% select(enst, all_of(shared_samples))


# Handle NAs and negatives (set to 0)
dataframe_tpm_edit <- dataframe_tpm_edit %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.) | . < 0, 0, .)))


# Step 2: Filter by median TPM threshold
thresh <- as.numeric(Sys.getenv("MIN_TPM", unset = "10"))
message("[INFO] Applying median TPM threshold ≥ ", thresh)


# Define subgroups (generalized; default to 'all' – add Hartwig-specific if purity has metadata)
case_groups <- list(
  all = dataframe_purity$case
  # Example: Add if needed, e.g., UCEC = dataframe_purity$case[dataframe_purity$cohort == "UCEC"]
)


tpm_median_final <- list()
transcript_passed <- list()


for (group in names(case_groups)) {
  samples <- intersect(case_groups[[group]], colnames(dataframe_tpm_edit)[-1])
  message("[INFO] Processing group '", group, "' with ", length(samples), " samples")
  
  df_sub <- dataframe_tpm_edit %>% select(enst, all_of(samples))
  median_df <- df_sub %>% rowwise() %>% mutate(median = median(c_across(-enst), na.rm = TRUE)) %>% ungroup() %>% select(enst, median)
  colnames(median_df)[2] <- group
  
  passed <- median_df %>% filter(!!sym(group) >= thresh) %>% pull(enst)
  
  transcript_passed[[group]] <- passed
  tpm_median_final[[group]] <- median_df
}


# Merge results
all_transcripts <- unique(unlist(transcript_passed))
tpm_merged <- reduce(tpm_median_final, full_join, by = "enst")


# Filter final TPM table
dataframe_tpm_final <- dataframe_tpm_edit %>% filter(enst %in% all_transcripts)


# Step 3: Generate summaries (ensure protein-coding match)
summary_raw <- dataframe_gtf %>%
  select(enst, symbol, ensg) %>%
  right_join(tpm_merged %>% mutate(enst = str_sub(enst, 1, 15)), by = "enst") %>%
  arrange(desc(all))


summary_passed <- summary_raw %>% filter(str_sub(enst, 1, 15) %in% str_sub(all_transcripts, 1, 15))


# Validation: Check if all filtered transcripts are protein-coding
non_coding <- anti_join(summary_passed %>% select(enst), dataframe_gtf %>% select(enst), by = "enst")
if (nrow(non_coding) > 0) {
  message("[WARN] ", nrow(non_coding), " filtered transcripts not found in protein-coding GTF.")
}


gtf_passed <- dataframe_gtf %>%
  semi_join(summary_passed, by = "enst") %>%
  select(chr = X1, source = X2, start = X4, end = X5, strand = X7, enst, symbol, ensg)


# Step 4: Output
setwd(output_dir)
current_date <- format(Sys.Date(), "%Y%m%d")


write_tsv(dataframe_tpm_final, paste0("TPM_Filter", thresh, "_ProteinCodingTx_", current_date, ".tsv"))
write_tsv(summary_raw, paste0("TPM_Summary_Table_Complete_", current_date, ".tsv"))
write_tsv(summary_passed, paste0("TPM_Summary_Table_Filter", thresh, "_", current_date, ".tsv"))
write_tsv(gtf_passed, paste0("GTF_ProteinCoding_Filter", thresh, "_", current_date, ".tsv"))


message("[INFO] Outputs written with threshold ", thresh, " and date ", current_date, ".")
