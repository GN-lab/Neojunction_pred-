#!/usr/bin/env Rscript
# Title: "Step 9: Count, Depth, Frequency, and Judgement Tables for GTEx"
# September 26, 2025 | Gaurav Raichand | The Institute of Cancer Research

# Note: SAME script as Step 7
# Purpose: Generate and prepare the following tables for the GTEx data:
#            1. Count table
#            2. Depth table
#            3. Frequency table
#            4. Judgement table

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(tidyverse)
library(readxl)
library(ggsci)
library(data.table)
current_date <- format(Sys.Date(), "%Y%m%d")

# Establish Directories ---------------------------------------------------

directory_external <- Sys.getenv("GTEX_ANNOTATION_PATH")
directory_07 <- Sys.getenv("STEP07_OUTPUT_DIR")
directory_08 <- Sys.getenv("STEP08_OUTPUT_DIR")
directory_out <- Sys.getenv("OUTPUT_DIR")

# Load Files --------------------------------------------------------------
# From External: GTEx Information Dataframe
setwd(directory_external)
filename_gtex.info <- "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
dataframe_gtex.info <- fread(filename_gtex.info, na = c("", "NA"))

# From Step 7: NJs that passed the PSR > 0.1 filter in Hartwig (dynamic)
setwd(directory_07)
psr_files <- list.files(directory_07,
                        pattern = "^PSR_Table_[0-9]{8}\\.tsv$",
                        full.names = TRUE)
if (length(psr_files) == 0) {
  stop("[ERROR] No PSR_Table_*.tsv found in ", directory_07)
}
latest_psr <- sort(psr_files)[length(psr_files)]
message("[INFO] Loading Hartwig PSR file: ", basename(latest_psr))
dataframe_passed_hartwig_junc <- fread(latest_psr, na = c("", "NA"))

# From Step 8: GTEx Overlap Table (dynamic) - detect latest
setwd(directory_08)
overlap_files <- list.files(directory_08,
                            pattern = "^GTEx_SJ_Overlap_Table_[0-9]{8}\\.tsv$",
                            full.names = TRUE)
if (length(overlap_files) == 0) {
  stop("[ERROR] No GTEx_SJ_Overlap_Table_*.tsv found in ", directory_08)
}
filename_overlap <- sort(overlap_files)[length(overlap_files)]
message("[INFO] Loading latest GTEx Overlap file: ", basename(filename_overlap))
dataframe_overlap.table.1 <- fread(filename_overlap, na = c("", "NA")) %>%
  inner_join(dataframe_passed_hartwig_junc, by = "junc.id") %>%
  select(junc.id, junc.id.overlap)

# From Step 8: GTEx Count Table (dynamic) - detect latest
count_files <- list.files(directory_08,
                          pattern = "^GTEx_SJ_Count_Table_SJ.out.tab_Retained_[0-9]{8}\\.tsv$",
                          full.names = TRUE)
if (length(count_files) == 0) {
  stop("[ERROR] No GTEx_SJ_Count_Table_SJ.out.tab_Retained_*.tsv found in ", directory_08)
}
filename_gtex.count.1 <- sort(count_files)[length(count_files)]
message("[INFO] Loading latest GTEx Count file: ", basename(filename_gtex.count.1))
dataframe_gtex.count.1 <- fread(filename_gtex.count.1, na = c("", "NA"))

# Check if Step 2 output files exist and are non-empty - detect latest set
setwd(directory_out)
count_files_step2 <- list.files(directory_out, pattern = "^Count_ALL_GTEx_[0-9]{8}\\.tsv$", full.names = TRUE)
depth_files_step2 <- list.files(directory_out, pattern = "^Depth_ALL_GTEx_[0-9]{8}\\.tsv$", full.names = TRUE)
freq_files_step2 <- list.files(directory_out, pattern = "^Freq_ALL_GTEx_[0-9]{8}\\.tsv$", full.names = TRUE)
judge_files_step2 <- list.files(directory_out, pattern = "^Judge_ALL_GTEx_[0-9]{8}\\.tsv$", full.names = TRUE)

# Find the latest date with all 4 files
all_dates <- unique(c(sub("^.*_([0-9]{8})\\.tsv$", "\\1", count_files_step2),
                      sub("^.*_([0-9]{8})\\.tsv$", "\\1", depth_files_step2),
                      sub("^.*_([0-9]{8})\\.tsv$", "\\1", freq_files_step2),
                      sub("^.*_([0-9]{8})\\.tsv$", "\\1", judge_files_step2)))
latest_date <- if (length(all_dates) > 0) sort(all_dates, decreasing = TRUE)[1] else NULL

files_exist_and_nonempty <- FALSE
if (!is.null(latest_date)) {
  filename_count <- file.path(directory_out, paste0("Count_ALL_GTEx_", latest_date, ".tsv"))
  filename_depth <- file.path(directory_out, paste0("Depth_ALL_GTEx_", latest_date, ".tsv"))
  filename_freq <- file.path(directory_out, paste0("Freq_ALL_GTEx_", latest_date, ".tsv"))
  filename_judge <- file.path(directory_out, paste0("Judge_ALL_GTEx_", latest_date, ".tsv"))
  
  files_step2 <- c(filename_count, filename_depth, filename_freq, filename_judge)
  files_exist_and_nonempty <- all(sapply(files_step2, function(f) file.exists(f) && file.info(f)$size > 0))
}

if (files_exist_and_nonempty) {
  message("[INFO] Existing Step 2 output files found and non-empty (using date ", latest_date, "), skipping Step 2 processing loop.")
  count <- fread(filename_count, header = TRUE, na = c("", "NA"))
  depth <- fread(filename_depth, header = TRUE, na = c("", "NA"))
  freq <- fread(filename_freq, header = TRUE, na = c("", "NA"))
  judge <- fread(filename_judge, header = TRUE, na = c("", "NA"))
} else {
  message("[INFO] Step 2 output missing or empty, running Step 2 processing loop.")

  ###########################################################################
  #  Step 2: Obtain the Count, Depth, Frequency, and Judgement Values -------
  ###########################################################################

  # Make lists for dataframes
  list_overlap.table <- list(dataframe_overlap.table.1)
  list_sj.count <- list(dataframe_gtex.count.1)

  list_count <- list()
  list_depth <- list()
  list_freq <- list()
  list_judge <- list()

  for (h in 1:length(list_overlap.table)) {
    dataframe_overlap.table.h <- list_overlap.table[[h]]
    dataframe_sj.count.h <- list_sj.count[[h]]
    for (i in 1:nrow(dataframe_overlap.table.h)) {
      # Status bar
      print(h)
      print(i / nrow(dataframe_overlap.table.h) * 100)
      
      # Select out the junction of interest (JUNC.ID) and the overlapping junctions (JUNC.ID.OVERLAP)
      JUNC.ID <- dataframe_overlap.table.h %>% slice(i) %>% pull(junc.id)
      JUNC.ID.OVERLAP <- dataframe_overlap.table.h %>% slice(i) %>% pull(junc.id.overlap)
      
      # If there are NO overlapping junctions
      if (is.na(JUNC.ID.OVERLAP)) {
        res.i <- dataframe_sj.count.h %>% 
          filter(junc.id == JUNC.ID) %>% 
          gather(case, value, 2:ncol(dataframe_sj.count.h)) %>% 
          spread(junc.id, value) %>% 
          rename(count = !!JUNC.ID) %>% 
          mutate(count = ifelse(is.na(count), 0, count)) %>% 
          mutate(depth = count) %>% 
          mutate(freq = round(count / depth, 4)) %>% 
          mutate(freq = ifelse(is.na(freq), 0, freq)) %>% 
          mutate(judge = ifelse(count >= 2 & depth >= 10 & freq >= 0.01, 1, 0)) %>% 
          mutate(junc.id = JUNC.ID) %>% 
          select(junc.id, case, count, depth, freq, judge)
      }
      
      # If there is an overlapping junction
      if (!is.na(JUNC.ID.OVERLAP)) {
        res.i <- bind_rows(dataframe_sj.count.h %>% filter(junc.id == JUNC.ID), 
                           dataframe_sj.count.h %>% filter(junc.id == JUNC.ID.OVERLAP)) %>% 
          gather(case, value, 2:ncol(dataframe_sj.count.h)) %>% 
          spread(junc.id, value) %>% 
          rename(count = !!JUNC.ID, count.o = !!JUNC.ID.OVERLAP) %>% 
          mutate(count = ifelse(is.na(count), 0, count)) %>% 
          mutate(count.o = ifelse(is.na(count.o), 0, count.o)) %>% 
          mutate(depth = count + count.o) %>% 
          mutate(freq = round(count / depth, 4)) %>% 
          mutate(freq = ifelse(is.na(freq), 0, freq)) %>% 
          mutate(judge = ifelse(count >= 2 & depth >= 10 & freq >= 0.01, 1, 0)) %>% 
          mutate(junc.id = JUNC.ID) %>% 
          select(junc.id, case, count, depth, freq, judge)
      }
      
      count.i <- res.i %>% select(junc.id, case, count) %>% spread(case, count)
      depth.i <- res.i %>% select(junc.id, case, depth) %>% spread(case, depth)
      freq.i <- res.i %>% select(junc.id, case, freq) %>% spread(case, freq)
      judge.i <- res.i %>% select(junc.id, case, judge) %>% spread(case, judge)
      
      if (i == 1) {
        count.h <- count.i
        depth.h <- depth.i
        freq.h <- freq.i
        judge.h <- judge.i
      } else {
        count.h <- bind_rows(count.h, count.i)
        depth.h <- bind_rows(depth.h, depth.i)
        freq.h <- bind_rows(freq.h, freq.i)
        judge.h <- bind_rows(judge.h, judge.i)
      }
    }
    list_count[[h]] <- count.h
    list_depth[[h]] <- depth.h
    list_freq[[h]] <- freq.h
    list_judge[[h]] <- judge.h
  }

  setwd(directory_out)
  write_tsv(list_count[[1]], filename_count, na = "NA", col_names = TRUE, escape = "double")
  write_tsv(list_depth[[1]], filename_depth, na = "NA", col_names = TRUE, escape = "double")
  write_tsv(list_freq[[1]], filename_freq, na = "NA", col_names = TRUE, escape = "double")
  write_tsv(list_judge[[1]], filename_judge, na = "NA", col_names = TRUE, escape = "double")
}

###########################################################################
#  Step 3: Prepare 3 Subgroups --------------------------------------------
###########################################################################

# Load judge table if not already in environment
# (Assumes `judge` was loaded above via fread)
available_samples <- colnames(judge)[-1]  # all 50 sample IDs in your judge table

# Define subgroups using only these available samples
CASES <- list(
  GTEx_all        = available_samples,
  GTEx_all_but_ts = available_samples[available_samples %in% 
                       dataframe_gtex.info$SAMPID[dataframe_gtex.info$SMTS != "Testis"]],
  Brain           = available_samples[available_samples %in% 
                       dataframe_gtex.info$SAMPID[dataframe_gtex.info$SMTS == "Brain"]]
)

# Positive Sample Rate (PSR)
list_psr  <- list()
list_pass <- list()
judge.2    <- NULL

for (i in seq_along(CASES)) {
  subgroup_name <- names(CASES)[i]
  samples_i     <- CASES[[i]]
  message("[INFO] Processing subgroup: ", subgroup_name, " (n=", length(samples_i), ")")
  
  # Subset the judge table and remove any duplicate junc.id/sample rows
  judge.i <- judge %>%
    select(junc.id, all_of(samples_i)) %>%
    distinct(junc.id, .keep_all = TRUE)   # ensure one row per junc.id
  
  if (ncol(judge.i) <= 1) {
    warning("No samples found for subgroup ", subgroup_name, "; skipping")
    next
  }
  
  judge.i <- judge.i %>%
    mutate(pos = rowSums(across(-junc.id))) %>%
    mutate(psr = round(pos / length(samples_i), 4)) %>%
    select(junc.id, psr)
  
  colnames(judge.i)[2] <- paste0(subgroup_name, "_n", length(samples_i))
  
  if (is.null(judge.2)) {
    judge.2 <- judge.i
  } else {
    judge.2 <- full_join(
      judge.2, judge.i,
      by           = "junc.id",
      relationship = "many-to-many"
    )
  }
}

# Apply frequency filter
judge.2 <- judge.2 %>%
  mutate(min = pmin(.[[2]], .[[3]]))  # psr columns for GTEx_all and GTEx_all_but_ts

judge.pass <- judge.2 %>%
  filter(min < as.numeric(Sys.getenv("MAX_GTEX_FREQUENCY", unset = "0.01")))

list_psr[[1]]  <- judge.2
list_pass[[1]] <- judge.pass

# Subset final tables
count.pass   <- count  %>% semi_join(judge.pass, by = "junc.id")
depth.pass   <- depth  %>% semi_join(judge.pass, by = "junc.id")
freq.pass    <- freq   %>% semi_join(judge.pass, by = "junc.id")
judge.pass2  <- judge  %>% semi_join(judge.pass, by = "junc.id")
overlap.pass <- dataframe_overlap.table.1 %>%
                  semi_join(judge.pass, by = "junc.id")

# Write outputs (use escape = "double")
write_tsv(count.pass,   paste0("Count_Retained_Passed_GTEx_", current_date, ".tsv"),
          na="NA", col_names=TRUE, escape="double")
write_tsv(depth.pass,   paste0("Depth_Retained_Passed_GTEx_", current_date, ".tsv"),
          na="NA", col_names=TRUE, escape="double")
write_tsv(freq.pass,    paste0("Freq_Retained_Passed_GTEx_", current_date, ".tsv"),
          na="NA", col_names=TRUE, escape="double")
write_tsv(judge.pass2,  paste0("Judge_Retained_Passed_GTEx_", current_date, ".tsv"),
          na="NA", col_names=TRUE, escape="double")
write_tsv(overlap.pass, paste0("OverlapTable_Retained_Passed_GTEx_", current_date, ".tsv"),
          na="NA", col_names=TRUE, escape="double")

# Also write PSR tables
write_tsv(judge.2,      paste0("PSR_ALL_GTEx_", current_date, ".tsv"),
          na="NA", col_names=TRUE, escape="double")
write_tsv(judge.pass,   paste0("PSR_Retained_Passed_GTEx_", current_date, ".tsv"),
          na="NA", col_names=TRUE, escape="double")
