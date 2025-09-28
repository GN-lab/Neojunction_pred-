#!/usr/bin/env Rscript
# Title: "Step 5b - Further Filter SJ for Unique Protein-Coding Transcripts with GTF and a PSR > 10"
# September 24, 2025 | Gaurav Raichand | The Institute of Cancer Research
# Last updated: 10-23-2020

# Part A - Purpose: Filter splicing junctions for:
#             1. Non-annotated
#             2. Read count > 10
# Part B - Purpose: Filter splicing junctions for:
#             3. Protein Coding
#             4. Uniquely Expressed
#             5. PSR > 10 in at least one Disease Group

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(tidyverse)
library(readxl)
library(ggsci)
library(data.table)  # For faster reading

input_dir <- Sys.getenv("INPUT_DIR")
output_dir <- Sys.getenv("OUTPUT_DIR")

# Directory for SJ RNA-seq data (from config, e.g., STAR_SJ_DIR for all Hartwig samples)
directory_sj <- Sys.getenv("STAR_SJ_DIR")

# Load non-annotated candidates from previous step (dynamic filename - using file.path for reliability)
current_date <- format(Sys.Date(), "%Y%m%d")
filename_dataframe_sj.combined <- paste0("SJ_List_NonAnnotated_Candidates_Protein_Coding_", current_date, ".tsv")
full_path_sj <- file.path(output_dir, filename_dataframe_sj.combined)

# Check if file exists before loading
if (!file.exists(full_path_sj)) {
  stop("[ERROR] File does not exist: ", full_path_sj, ". Check if previous script ran successfully and output_dir is correct.")
}
dataframe_sj.combined <- fread(full_path_sj, colClasses = list(character = "chr"))

# Load filtered GTF from Step 3 (assuming it's in output_dir based on ls)
tpm_thresh <- as.numeric(Sys.getenv("MIN_TPM", unset = "10"))
filename_gtf <- paste0("GTF_ProteinCoding_Filter", tpm_thresh, "_", current_date, ".tsv")
full_path_gtf <- file.path(output_dir, filename_gtf)
if (!file.exists(full_path_gtf)) {
  stop("[ERROR] File does not exist: ", full_path_gtf)
}
dataframe_gtf <- fread(full_path_gtf, colClasses = list(character = "chr"))

# Load purity-filtered patient list from Step 1 (assuming it's in output_dir based on ls) - fixed formatting to "0.60"
purity_thresh <- as.numeric(Sys.getenv("MIN_PURITY", unset = "0.60"))
filename_tcga <- paste0("Patient_List_Post_TumorPurity_Filter_", format(purity_thresh, nsmall = 2), ".txt")
full_path_tcga <- file.path(output_dir, filename_tcga)
if (!file.exists(full_path_tcga)) {
  stop("[ERROR] File does not exist: ", full_path_tcga)
}
dataframe_tcga <- fread(full_path_tcga) %>% rename(case = sample_id)  # Standardize

###########################################################################
#  Step 1: Filter for SJ in Protein Coding Transcripts --------------------
###########################################################################
# The previous list of splicing junctions generated are filtered for:
#         1. Are non-annotated
#         2. Have a read count >= 10
# The goal of this step is to further filter for splicing junctions for:
#         3. Within protein-coding transcript regions

start.time <- proc.time()

# Initialize as empty to avoid errors if no rows qualify
dataframe_sj.filtered <- tibble(junc.id = character(), chr = character(), strand = character(), int.start = numeric(), int.end = numeric())

if (nrow(dataframe_sj.combined) == 0) {
  message("[WARNING] No non-annotated candidates to filter for protein-coding.")
} else {
  for (i in 1:nrow(dataframe_sj.combined)) {
    # Progress Bar
    message("[PROGRESS] Filtering: ", round((i / nrow(dataframe_sj.combined)) * 100, 2), "%")
    
    # Iterate through each row to work with them individually
    dataframe_sj.combined.i <- dataframe_sj.combined %>% 
      slice(i)
    
    CHR <- dataframe_sj.combined.i %>% pull(chr)
    STRAND <- dataframe_sj.combined.i %>% pull(strand)
    START <- dataframe_sj.combined.i %>% pull(int.start)
    END <- dataframe_sj.combined.i %>% pull(int.end)
    
    # 1. Generate a new variable RETAIN that is 1 for retain or 0 for do NOT retain
    RETAIN <- dataframe_gtf %>% 
      filter(chr == CHR & strand == STRAND) %>%
      filter(start < START & END < end) %>%
      nrow()
    
    # 2. Append only if RETAIN >= 1
    if (RETAIN >= 1) {
      dataframe_sj.filtered <- bind_rows(dataframe_sj.filtered, dataframe_sj.combined.i)
    }
  }
}

RUNTIME <- proc.time() - start.time
print(RUNTIME)

# 3. Retain splicing junctions with RETAIN >= 1
dataframe_sj.combined.retain <- dataframe_sj.filtered %>% 
  select(junc.id, chr, strand, int.start, int.end)

###########################################################################
#  Step 2: Output SJ Data (Filtered for TPM and PC) -----------------------
###########################################################################

filename_output.sj.filtered <- paste0("SJ_List_NonAnnotated_Candidates_Protein_Coding_", current_date, ".tsv")
full_path_output_sj <- file.path(output_dir, filename_output.sj.filtered)
write_tsv(dataframe_sj.combined.retain, full_path_output_sj, na = "NA", col_names = TRUE, quote_escape = "double")

###########################################################################
#  Step 3: Generate SJ Count Table from sj.out.tab with the SJ List -------
###########################################################################
# From the list of splicing junctions generated in Step 1, this step generates the "count" 
# by overlaying the list with the sample SJ files

# Reload the output we just wrote (for consistency)
dataframe_sj.combined.retain <- fread(full_path_output_sj, colClasses = list(character = "chr"))

# 1. Initiate the count dataframe by taking only the junction coordinates 
dataframe_sj.count <- dataframe_sj.combined.retain %>% 
  select(junc.id)

start.time <- proc.time()

# 2. Filter for the count data using the junction IDs and cross analyzing it with the sample SJ files
setwd(directory_sj)  # Single directory for Hartwig
all_files <- list.files(pattern = "\\.SJ\\.out\\.tab$")  # Assume standard STAR naming
sample_ids <- gsub("\\.SJ\\.out\\.tab$", "", all_files)
matching_indices <- sample_ids %in% dataframe_tcga$case
list_files <- all_files[matching_indices]

if (length(list_files) == 0) {
  message("[WARNING] No matching SJ files found in directory: ", directory_sj)
} else {
  for (j in 1:length(list_files)) {
    in.f <- list_files[j]
    sample_id <- gsub("\\.SJ\\.out\\.tab$", "", in.f)  # Extract clean sample ID for column name
    dataframe_sj.combined.j <- fread(in.f, col.names = c("chr", "int.start", "int.end", "strand", "int.motif", "annot", "n.uniq.map", "n.mult.map", "max.spl.overhang"),
                                     colClasses = list(character = "chr"))
    
    # 3. Retain junctions that are uniquely expressed
    dataframe_sj.count <- dataframe_sj.count %>% 
      left_join(dataframe_sj.combined.j %>% 
                  mutate(strand = ifelse(strand == 2, "-", ifelse(strand == 1, "+", "undefined"))) %>% 
                  filter(strand != "undefined") %>%
                  mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end)) %>% 
                  select(junc.id, n.uniq.map), by = "junc.id") %>% 
      mutate(n.uniq.map = ifelse(is.na(n.uniq.map), 0, n.uniq.map))
    
    colnames(dataframe_sj.count)[ncol(dataframe_sj.count)] <- sample_id  # Use clean sample ID
    
    message("[PROGRESS] Counting: ", round((j / length(list_files)) * 100, 4), "%")
  }
}

RUNTIME <- proc.time() - start.time
print(RUNTIME)

###########################################################################
#  Step 4: Filter for Junctions with Read Count > 10 per Disease Group ----
###########################################################################
# Filter for junctions that have a minimum Read Count of 10 in at least MIN_SAMPLES_PCT of each group

# 1. Manipulate the dataframe to have all the counts for each junction+sample in one column
dataframe_sj.count.warp <- dataframe_sj.count %>% 
  gather("sample.id", "count", 2:ncol(dataframe_sj.count)) %>% 
  mutate(judge = ifelse(count < 10, 0, 1)) %>% 
  select(-count)

# 2. Filter for minimum read count >= 10 in at least MIN_SAMPLES_PCT within each group
# Generalized to "all" group (add more if needed)
case_groups <- list(all = dataframe_tcga$case)

dataframe_sj.psr <- NULL  # Initialize as NULL

for (i in 1:length(case_groups)) {
  message("[PROGRESS] Group: ", i)
  CASES <- case_groups[[i]]
  
  dataframe_sj.count.case <- dataframe_sj.count.warp %>% 
    filter(sample.id %in% CASES) %>% 
    spread(sample.id, judge)
  
  dataframe_sj.psr.i <- dataframe_sj.count.case %>% 
    mutate(n.pos = rowSums(select(., -junc.id), na.rm = TRUE)) %>% 
    mutate(pos.sample.rate = round(n.pos / length(CASES), 3)) %>% 
    select(junc.id, pos.sample.rate)
  
  colnames(dataframe_sj.psr.i)[2] <- paste0("all_n", length(CASES))  # Generalized
  
  if (i == 1) {
    dataframe_sj.psr <- dataframe_sj.psr.i
  } else {
    dataframe_sj.psr <- dataframe_sj.psr %>% 
      bind_cols(dataframe_sj.psr.i %>% select(-junc.id))
  }
}

# 3. Determine the junctions to retain by whether max PSR >= MIN_SAMPLES_PCT
min_pct <- as.numeric(Sys.getenv("MIN_SAMPLES_PCT", unset = "0.10"))
dataframe_sj.psr <- dataframe_sj.psr %>% 
  mutate(max = apply(select(., -junc.id), 1, max, na.rm = TRUE)) %>% 
  mutate(retain = ifelse(max >= min_pct, 1, 0))

dataframe_sj.psr.retain <- dataframe_sj.psr %>% 
  filter(retain == 1) %>% 
  select(-retain)

# Generate the Junction Count List for retained
dataframe_sj.count.retain <- dataframe_sj.count %>% 
  semi_join(dataframe_sj.psr.retain, by = "junc.id")

###########################################################################
#  Step 5: Output Data ----------------------------------------------------
###########################################################################

filename_sj.retained.psr <- paste0("SJ_PSR_NonAnnotated_Candidates_Protein_Coding_", current_date, ".tsv")
full_path_psr <- file.path(output_dir, filename_sj.retained.psr)
write_tsv(dataframe_sj.psr.retain, full_path_psr, na = "NA", col_names = TRUE, quote_escape = "double")

filename_sj.retained.count <- paste0("SJ_Count_NonAnnotated_Candidates_Protein_Coding_", current_date, ".tsv")
full_path_count <- file.path(output_dir, filename_sj.retained.count)
write_tsv(dataframe_sj.count.retain, full_path_count, na = "NA", col_names = TRUE, quote_escape = "double")
