#!/usr/bin/env Rscript
# Title: "Step 6: SJ Overlap Table"
# September 24, 2025 | Gaurav Raichand | The Institute of Cancer Research

# Purpose: Generate and prepare an Overlap Table for each candidate splicing junction event
#          This step requires both the annotated (Step 4) and non-annotated (Step 5) SJs in order
#          to identify ALL overlapping junctions

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(tidyverse)
library(readxl)
library(ggsci)
library(data.table)

# Establish Directories ---------------------------------------------------

input_dir <- Sys.getenv("INPUT_DIR")
directory_samples <- Sys.getenv("STAR_SJ_DIR")  # Single directory for Hartwig SJ files

directory_06 <- Sys.getenv("OUTPUT_DIR")  # For outputs

# Load Files --------------------------------------------------------------
# From Step 4: SJ ID's (Annotated) – dynamic filename, load from output_dir based on ls
current_date <- format(Sys.Date(), "%Y%m%d")
filename_sj.annot <- paste0("SJ_List_Filtered_by_GTF_ProteinCoding_ExpressedTranscripts_", current_date, ".tsv")
full_path_annot <- file.path(directory_06, filename_sj.annot)  # Use output_dir
if (!file.exists(full_path_annot)) {
  stop("[ERROR] File does not exist: ", full_path_annot, ". Check output_dir and if Step 4 ran successfully.")
}
dataframe_sj.annot <- fread(full_path_annot, colClasses = list(character = "chr"))

# From Step 5b: SJ (Non-Annotated with Read Count > 10, PSR > 10, and Protein Coding) – dynamic, load from output_dir
filename_sj.nonannot <- paste0("SJ_List_NonAnnotated_Candidates_Protein_Coding_", current_date, ".tsv")
full_path_nonannot <- file.path(directory_06, filename_sj.nonannot)  # Use output_dir
if (!file.exists(full_path_nonannot)) {
  stop("[ERROR] File does not exist: ", full_path_nonannot, ". Check output_dir and if Step 5b ran successfully.")
}
dataframe_sj.nonannot <- fread(full_path_nonannot, colClasses = list(character = "chr"))

###########################################################################
#  Step 1: Generate a Complete List of Junction ID's to Analyze -----------
###########################################################################

# 1. Take the Annotated SJ List and edit it such that it matches the Non-Annotated SJ List
#    with "junc.id", "strand", "int.start", "int.end"

dataframe_sj.nonannot <- dataframe_sj.nonannot %>% 
  select(junc.id) %>% 
  mutate(chr = gsub("chr", "", sapply(strsplit(junc.id, ":"), "[[", 1)),
         strand = sapply(strsplit(junc.id, ":"), "[[", 2),
         int.start = as.numeric(gsub("-.*", "", sapply(strsplit(junc.id, ":"), "[[", 3))) + 1,
         int.end = as.numeric(gsub(".*-", "", sapply(strsplit(junc.id, ":"), "[[", 3))))

# 2. Identify which "annotated" and "non-annotated" junction counts to collect
start.time <- proc.time()

list_sj <- list()
junc.ids <- NULL

for (i in 1:nrow(dataframe_sj.nonannot)) {
  
  # Progress Bar
  message("[PROGRESS] Step 1: ", round((i / nrow(dataframe_sj.nonannot)) * 100, 2), "%")
  
  # Slice out row(i) for every iteration in the for-loop
  sj.i <- dataframe_sj.nonannot %>% slice(i)
  
  # Extract the JUNC.ID, STRAND, START, and END for junction(i)
  CHR <- sj.i %>% pull(chr)
  JUNC.ID <- sj.i %>% pull(junc.id) 
  STRAND <- sj.i %>% pull(strand) 
  START <- sj.i %>% pull(int.start) 
  END <- sj.i %>% pull(int.end) 
  
  # Find all of the overlapping junctions from the "annotated" SJ list
  JUNC.ID.OVERLAP <- dataframe_sj.annot %>% 
    filter(chr == CHR, strand == STRAND & int.start < END & START < int.end) %>% 
    pull(junc.id)
  
  # Create a list of overlapping junctions for each of the "non-annotated" junctions 
  if (length(JUNC.ID.OVERLAP) == 0) {
    list_sj[[i]] <- tibble(
      junc.id = JUNC.ID, 
      junc.id.overlap = NA)
  } else {
    list_sj[[i]] <- tibble(
      junc.id = JUNC.ID, 
      junc.id.overlap = JUNC.ID.OVERLAP)
  }
  
  junc.ids.i <- list_sj[[i]] %>% 
    gather("label", "junc.id", 1:2) %>% 
    filter(!is.na(junc.id)) %>% 
    distinct(junc.id) %>% 
    pull(junc.id)
  
  junc.ids <- unique(c(junc.ids, junc.ids.i))
}

RUNTIME <- proc.time() - start.time 
print(RUNTIME)

head(junc.ids)

###########################################################################
#  Step 2: Analyze the Count Table with the Complete Junction ID List -----
###########################################################################

# 1. Generate an ordered list of junctions to bind to the sj.out.tab files
# 1a. Create the dataframe for SJ Counts with the first column being all of the Junction ID's
#     from the junc.id list made from Step 1
dataframe_sj.count <- tibble(junc.id = junc.ids) %>% 
  mutate(start = as.numeric(gsub("-.*", "", sapply(strsplit(junc.id, ":"), "[[", 3)))) %>% 
  arrange(start) %>% 
  select(junc.id)

dim(dataframe_sj.count) %>% print()

# 2. Bind the ordered list of junctions to the sj.out.tab files
setwd(directory_samples)
list_rnaseq.files <- list.files(pattern = "SJ.out.tab$")  # Assuming files end with SJ.out.tab

start.time <- proc.time()

for (j in 1:length(list_rnaseq.files)) {
  # 2b. Work with each file in the directory iteratively
  filename_j <- list_rnaseq.files[j]
  dataframe_sj.j <- fread(filename_j, col.names = c("chr", "int.start", "int.end", "strand", "int.motif", "annot", "n.uniq.map", "n.mult.map", "max.spl.overhang"),
                          colClasses = list(character = "chr"))
  
  append_i <- dataframe_sj.j %>% 
    mutate(strand = ifelse(strand == 2, "-", ifelse(strand == 1, "+", "undefined"))) %>% 
    filter(strand != "undefined") %>%
    mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end)) %>% 
    select(junc.id, n.uniq.map)
  
  dataframe_sj.count <- left_join(dataframe_sj.count, append_i, by = "junc.id") %>% 
    mutate(n.uniq.map = ifelse(is.na(n.uniq.map), 0, n.uniq.map))
  
  colnames(dataframe_sj.count)[ncol(dataframe_sj.count)] <- gsub("SJ.out.tab", "", filename_j)  # Use full sample name
  
  # Status Bar
  message("[PROGRESS] Step 2: ", round(((ncol(dataframe_sj.count) - 1) / length(list_rnaseq.files)) * 100, 2), "%")
}

RUNTIME <- proc.time() - start.time 
print(RUNTIME)

# 3. Identify the "Most-Dominant" Junctions Among the Multiple Overlaps
start.time <- proc.time()

dataframe_overlap.table <- NULL

for (i in 1:length(list_sj)) {
  if (nrow(list_sj[[i]]) == 1) {
    dataframe_overlap.table <- bind_rows(dataframe_overlap.table, list_sj[[i]])
  } else {
    sj_temp <- list_sj[[i]] %>% 
      select(junc.id.overlap) %>% 
      left_join(dataframe_sj.count, by = c("junc.id.overlap" = "junc.id")) %>% 
      mutate(sum = rowSums(select(., -junc.id.overlap), na.rm = TRUE)) %>% 
      select(junc.id.overlap, sum) %>% 
      arrange(desc(sum)) %>% 
      head(n = 1)
    
    dataframe_overlap.table <- bind_rows(dataframe_overlap.table, 
                                         list_sj[[i]] %>% semi_join(sj_temp, by = "junc.id.overlap"))
    
    # Progress Report
    message("[PROGRESS] Step 3: ", round((i / length(list_sj)) * 100, 2), "%")
  }
}

RUNTIME <- proc.time() - start.time 
print(RUNTIME)

###########################################################################
#  Step 3: Edit the SJ Count Table ----------------------------------------
###########################################################################

# Retain the count data for the overlapping junctions
dataframe_sj.count.retain <- dataframe_sj.count %>% 
  right_join(dataframe_overlap.table %>% 
               gather("label", "junc.id", 1:2) %>% 
               filter(!is.na(junc.id)) %>% 
               distinct(junc.id),
             by = "junc.id")

###########################################################################
#  Step 4: Export Files ---------------------------------------------------
###########################################################################

filename_output.overlap.table <- paste0("SJ_Overlap_Table_", current_date, ".tsv")
full_path_overlap <- file.path(directory_06, filename_output.overlap.table)
write_tsv(dataframe_overlap.table, full_path_overlap, na = "NA", col_names = TRUE, quote_escape = "double")

filename_output.retain <- paste0("SJ_Count_Table_SJ.out.tab_Retained_", current_date, ".tsv")
full_path_retain <- file.path(directory_06, filename_output.retain)
write_tsv(dataframe_sj.count.retain, full_path_retain, na = "NA", col_names = TRUE, quote_escape = "double")
