#!/usr/bin/env Rscript
# Title: "Step 5 - Select for Non-annotated Junctions with a Read Count > 10"
# September 24, 2025 | Gaurav Raichand | The Institute of Cancer Research

# Part A - Purpose: Filter splicing junctions for two things:
#             1. Non-annotated
#             2. Read count > 10
# Pipeline: 1. Filter the folder directory for patients/samples in the purity list generated in Step 1 (Tumor Purity > 0.60)
#           2. From each file in the directory, select for junctions that have a read count > 10
#           3. Combine the all junctions from each file into a comprehensive list of Junctions with Read Count > 10

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(readxl)
library(tidyverse)
library(ggsci)
library(data.table)  # For faster file reading

input_dir <- Sys.getenv("INPUT_DIR")
output_dir <- Sys.getenv("OUTPUT_DIR")

# Directory for SJ RNA-seq data (from config, e.g., STAR_SJ_DIR for all Hartwig samples)
directory_sj <- Sys.getenv("STAR_SJ_DIR")

# Load Files --------------------------------------------------------------

# Splicing Junction Coordinates (Raw)
filename_sj <- "sjdbList.fromGTF.out.tab"
filename_sj_full <- file.path(input_dir, filename_sj)
message("[INFO] Using SJ file: ", filename_sj_full)
dataframe_sj <- fread(filename_sj_full, header = FALSE)
setnames(dataframe_sj, old = 1:4, new = c("chr", "int.start", "int.end", "strand"))

# Patient List (From Step 1) Samples with a Tumor Purity >= MIN_PURITY
purity_thresh <- as.numeric(Sys.getenv("MIN_PURITY", unset = "0.60"))
filename_tcga <- file.path(output_dir, paste0("Patient_List_Post_TumorPurity_Filter_", format(purity_thresh, nsmall = 2), ".txt"))
dataframe_tcga <- fread(filename_tcga, colClasses = list(character = "sample_id")) %>% rename(case = sample_id)  # Standardize

# GTF File - (From Step 3) Protein-Coding Transcripts with Median TPM >= MIN_TPM
tpm_thresh <- as.numeric(Sys.getenv("MIN_TPM", unset = "10"))
filename_gtf <- file.path(output_dir, paste0("GTF_ProteinCoding_Filter", tpm_thresh, "_", format(Sys.Date(), "%Y%m%d"), ".tsv"))
dataframe_gtf <- fread(filename_gtf, colClasses = list(character = "chr"))

head(dataframe_tcga, 3) %>% print()

###########################################################################
#  Step 1: Rename the Columns of the SJ Coordinates Table -----------------
###########################################################################

dataframe_sj <- dataframe_sj %>% 
  mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end))

###########################################################################
#  Step 2: Generate the List of Non-Annotated SJ --------------------------
###########################################################################
# Goal: This step generates a list of splicing junctions that:
#         1. Are non-annotated
#         2. Have a read count >= MIN_READ_COUNT

start.time <- proc.time()

# Create an empty list to hold all of the "non-annotated" junction counts
list_sj <- list()

# Initialize dataframe_combined as empty to avoid "not found" error if no files
dataframe_combined <- tibble(chr = character(), int.start = numeric(), int.end = numeric(), strand = character(), junc.id = character())

# Set to SJ directory (single for Hartwig; add more to list if needed)
groups_dirs <- list(directory_sj)

for (i in seq_along(groups_dirs)) {
  setwd(groups_dirs[[i]])
  
  # 1. Find all potential SJ files and extract sample IDs for matching
  all_files <- list.files(pattern = "\\.SJ\\.out\\.tab$")  # Assume standard STAR naming
  sample_ids <- gsub("\\.SJ\\.out\\.tab$", "", all_files)
  matching_indices <- sample_ids %in% dataframe_tcga$case
  list_rnaseq.files <- all_files[matching_indices]
  
  if (length(list_rnaseq.files) == 0) {
    message("[WARNING] No matching SJ files found in directory: ", groups_dirs[[i]])
    next
  }
  
  for (j in seq_along(list_rnaseq.files)) {
    filename_j <- list_rnaseq.files[j]
    dataframe_j <- fread(filename_j, col.names = c("chr", "int.start", "int.end", "strand", "int.motif", "annot", "n.uniq.map", "n.mult.map", "max.spl.overhang"),
                         colClasses = list(character = "chr"))
    
    # 2. Remove splicing junctions not in the raw reference splicing junction file
    list_sj[[length(list_sj) + 1]] <- dataframe_j %>% 
      mutate(strand = ifelse(strand == 2, "-", ifelse(strand == 1, "+", "undefined"))) %>% 
      filter(strand != "undefined") %>%
      filter(chr != "M") %>% 
      mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end)) %>% 
      anti_join(dataframe_sj, by = "junc.id")
    
    # Extract junctions with supportive read counts >= MIN_READ_COUNT
    min_read <- as.numeric(Sys.getenv("MIN_READ_COUNT", unset = "10"))
    dataframe_j.read10 <- list_sj[[length(list_sj)]] %>%
      filter(n.uniq.map >= min_read) %>%
      select(chr, int.start, int.end, strand, junc.id)
    
    # Combine the junctions (append to existing)
    dataframe_combined <- dataframe_combined %>% 
      bind_rows(dataframe_j.read10) %>% 
      distinct(junc.id, .keep_all = TRUE)
    
    # Progress Bar
    message("[PROGRESS] Group ", i, ": ", round((j / length(list_rnaseq.files)) * 100, 2), "%")
  }
}

# 3. Edit and finalize the combined SJ dataframe (safe even if empty)
dataframe_combined <- dataframe_combined %>% 
  arrange(chr, int.start, int.end)

RUNTIME <- proc.time() - start.time 
print(RUNTIME)

###########################################################################
#  Step 3: Output Data ----------------------------------------------------
###########################################################################

setwd(output_dir)

filename_output.sj.filtered <- paste0("SJ_List_NonAnnotated_Candidates_Protein_Coding_", format(Sys.Date(), "%Y%m%d"), ".tsv")
write_tsv(dataframe_combined, filename_output.sj.filtered, na = "NA", col_names = TRUE, quote_escape = "double")

# Part B - Purpose: Further filter for protein-coding and PSR >= MIN_SAMPLES_PCT

# Load the combined non-annotated candidates
setwd(output_dir)
filename_sj.combined <- filename_output.sj.filtered  # Fixed typo
dataframe_sj.combined <- fread(filename_sj.combined, colClasses = list(character = "chr"))

###########################################################################
#  Step 4: Filter for Splicing Junctions in Protein Coding Transcripts ----
###########################################################################

start.time <- proc.time()

# Initialize as empty to avoid errors if no rows qualify
dataframe_sj.filtered <- tibble(chr = character(), int.start = numeric(), int.end = numeric(), strand = character(), junc.id = character())

if (nrow(dataframe_sj.combined) == 0) {
  message("[WARNING] No non-annotated candidates to filter for protein-coding.")
} else {
  for (i in seq_len(nrow(dataframe_sj.combined))) {
    
    # Progress Bar
    message("[PROGRESS] Filtering: ", round((i / nrow(dataframe_sj.combined)) * 100, 2), "%")
    
    # Iterate through each row to work with them individually
    dataframe_sj.combined.i <- dataframe_sj.combined[i]
    
    CHR <- dataframe_sj.combined.i$chr
    STRAND <- dataframe_sj.combined.i$strand
    START <- dataframe_sj.combined.i$int.start
    END <- dataframe_sj.combined.i$int.end
    
    # Generate a new variable RETAIN that is 1 for retain or 0 for do NOT retain
    RETAIN <- dataframe_gtf %>% 
      filter(chr == CHR & strand == STRAND) %>%
      filter(start < START & END < end) %>%
      nrow()
    
    # Append only if RETAIN >= 1
    if (RETAIN >= 1) {
      dataframe_sj.filtered <- bind_rows(dataframe_sj.filtered, dataframe_sj.combined.i)
    }
  }
}

RUNTIME <- proc.time() - start.time 
print(RUNTIME)

###########################################################################
#  Step 5: Filter for Junctions with Read Count > 10 per Disease Group ----
###########################################################################

# 1. Initiate the count dataframe by taking only the junction coordinates 
dataframe_sj.count <- dataframe_sj.filtered %>% 
  select(junc.id)

# 2. Filter for the count data using the junction IDs and cross analyzing it with the sample SJ files
start.time <- proc.time()

# Use single SJ directory for Hartwig
setwd(directory_sj)
all_files <- list.files(pattern = "\\.SJ\\.out\\.tab$")
sample_ids <- gsub("\\.SJ\\.out\\.tab$", "", all_files)
matching_indices <- sample_ids %in% dataframe_tcga$case
list_rnaseq.files <- all_files[matching_indices]

if (length(list_rnaseq.files) == 0) {
  message("[WARNING] No matching SJ files for count filtering.")
} else {
  for (i in seq_along(list_rnaseq.files)) {
    in.f <- list_rnaseq.files[i]
    sample_id <- gsub("\\.SJ\\.out\\.tab$", "", in.f)  # Extract clean sample ID for column name
    sj.i <- fread(in.f, col.names = c("chr", "int.start", "int.end", "strand", "int.motif", "annot", "n.uniq.map", "n.mult.map", "max.spl.overhang"),
                  colClasses = list(character = "chr")) %>% 
      mutate(strand = ifelse(strand == 2, "-", ifelse(strand == 1, "+", "undefined"))) %>% 
      filter(strand != "undefined") %>% 
      mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end)) %>% 
      select(junc.id, n.uniq.map)
    
    dataframe_sj.count <- dataframe_sj.count %>% 
      left_join(sj.i, by = "junc.id") %>% 
      mutate(n.uniq.map = ifelse(is.na(n.uniq.map), 0, n.uniq.map))
    
    colnames(dataframe_sj.count)[ncol(dataframe_sj.count)] <- sample_id  # Use clean sample ID
    
    message("[PROGRESS] Counting: ", round((i / length(list_rnaseq.files)) * 100, 2), "%")
  }
}

RUNTIME <- proc.time() - start.time 
print(RUNTIME)

###########################################################################
#  Step 6: Filter for Junctions with Read Count > 10 per Disease Group ----
###########################################################################

# 1. Manipulate the dataframe to have all the counts for each junction+sample in one column
dataframe_sj.count.warp <- dataframe_sj.count %>% 
  gather("sample.id", "count", 2:ncol(dataframe_sj.count)) %>% 
  mutate(judge = ifelse(count < 10, 0, 1)) %>% 
  select(-count)

# 2. Filter for functions with a minimum read count >= 10 in at least MIN_SAMPLES_PCT within each group
# Generalized to "all" group (add more if dataframe_tcga has metadata columns)
case_groups <- list(all = dataframe_tcga$case)

dataframe_sj.psr <- NULL

for (i in seq_along(case_groups)) {
  CASES <- case_groups[[i]]
  
  dataframe_sj.count.case <- dataframe_sj.count.warp %>% 
    filter(sample.id %in% CASES) %>% 
    spread(sample.id, judge)
  
  dataframe_sj.psr.i <- dataframe_sj.count.case %>% 
    mutate(sum_positive.judgement = rowSums(select(., -junc.id), na.rm = TRUE)) %>% 
    mutate(positive.sample.rate = round(sum_positive.judgement / length(CASES), 3)) %>% 
    select(junc.id, positive.sample.rate)
  
  colnames(dataframe_sj.psr.i)[2] <- paste0(names(case_groups)[i], "_n", length(CASES))
  
  if (i == 1) {
    dataframe_sj.psr <- dataframe_sj.psr.i
  } else {
    dataframe_sj.psr <- dataframe_sj.psr %>% 
      full_join(dataframe_sj.psr.i, by = "junc.id")
  }
  
  print(i)
}

# 3. Determine the junctions to retain by whether max PSR >= MIN_SAMPLES_PCT
min_pct <- as.numeric(Sys.getenv("MIN_SAMPLES_PCT", unset = "0.10"))
dataframe_sj.psr <- dataframe_sj.psr %>% 
  mutate(max = apply(select(., -junc.id), 1, max, na.rm = TRUE)) %>% 
  mutate(retain = ifelse(max >= min_pct, 1, 0))

dataframe_sj.psr.retained <- dataframe_sj.psr %>% 
  filter(retain == 1) %>% 
  select(-retain)

table(dataframe_sj.psr$retain) %>% print()

###########################################################################
#  Step 7: Output Data ----------------------------------------------------
###########################################################################

setwd(output_dir)

filename_output.sj.filtered <- paste0("SJ_List_NonAnnotated_Candidates_Protein_Coding_", format(Sys.Date(), "%Y%m%d"), ".tsv")
write_tsv(dataframe_sj.filtered, filename_output.sj.filtered, na = "NA", col_names = TRUE, quote_escape = "double")  # Fixed to write the protein-coding filtered version

filename_output.sj.psr <- paste0("SJ_List_NonAnnotated_ProteinCoding_PSR_", format(Sys.Date(), "%Y%m%d"), ".tsv")
write_tsv(dataframe_sj.psr.retained, filename_output.sj.psr, na = "NA", col_names = TRUE, quote_escape = "double")

message("[INFO] Written filtered SJ to: ", filename_output.sj.filtered)
message("[INFO] Written PSR retained SJ to: ", filename_output.sj.psr)
