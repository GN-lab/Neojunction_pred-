#!/usr/bin/env Rscript
# Title: "Step 10 Extract Neojunctions"
# September 27, 2025 | Gaurav Raichand | The Institute of Cancer Research

# Purpose: Extract Neojunctions by filtering with GTEx junctions that have a PSR 
#          of < 0.01 (list derived from Step 09)

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(dplyr)
library(tidyverse)
library(readxl)
library(ggsci)
library(data.table)  # For faster reading

# Establish Directories ---------------------------------------------------

directory_external   <- Sys.getenv("META_FILES_PATH")
directory_out.step07 <- Sys.getenv("STEP07_OUTPUT_DIR")
directory_out.step09 <- Sys.getenv("STEP09_OUTPUT_DIR")
directory_out.step10 <- Sys.getenv("OUTPUT_DIR")

###########################################################################
#  Step 1: Load Files with auto-detect for date mismatch -----------------
###########################################################################

# Step 7 PSR Table
setwd(directory_out.step07)
psr_files_hartwig <- list.files(pattern = "^PSR_Table_[0-9]{8}\\.tsv$")
if (length(psr_files_hartwig) == 0) stop("No PSR_Table_*.tsv found in STEP07_OUTPUT_DIR")
latest_psr_hartwig <- sort(psr_files_hartwig)[length(psr_files_hartwig)]
dataframe_psr.hartwig <- fread(latest_psr_hartwig, na = c("", "NA"))

# Step 9 GTEx PSR Table
setwd(directory_out.step09)
psr_files_gtex <- list.files(pattern = "^PSR_Retained_Passed_GTEx_[0-9]{8}\\.tsv$")
if (length(psr_files_gtex) == 0) stop("No PSR_Retained_Passed_GTEx_*.tsv found in STEP09_OUTPUT_DIR")
latest_psr_gtex <- sort(psr_files_gtex)[length(psr_files_gtex)]
dataframe_psr.gtex <- fread(latest_psr_gtex, na = c("", "NA"))

# Step 7 Count/Depth/Freq/Judge Tables
setwd(directory_out.step07)
count_file <- sort(list.files(pattern = "^Count_Table_Retained_and_Passed_Junctions_[0-9]{8}\\.tsv$"))[1]
depth_file <- sort(list.files(pattern = "^Depth_Table_Retained_and_Passed_Junctions_[0-9]{8}\\.tsv$"))[1]
freq_file  <- sort(list.files(pattern = "^Freq_Table_Retained_and_Passed_Junctions_[0-9]{8}\\.tsv$"))[1]
judge_file <- sort(list.files(pattern = "^Judgement_Table_Retained_and_Passed_Junctions_[0-9]{8}\\.tsv$"))[1]

dataframe_count <- fread(count_file, na = c("", "NA"))
dataframe_depth <- fread(depth_file, na = c("", "NA"))
dataframe_freq  <- fread(freq_file,  na = c("", "NA"))
dataframe_judge <- fread(judge_file, na = c("", "NA"))

###########################################################################
#  Step 2: Filter for Neojunctions: PSR -----------------------------------
###########################################################################

dataframe_psr.neo <- right_join(
  dataframe_psr.hartwig %>% rename(max.hartwig = max),
  dataframe_psr.gtex %>% rename(min.gtex = min),
  by = "junc.id"
)

setwd(directory_out.step10)
current_date <- format(Sys.Date(), "%Y%m%d")
filename_psr.p1 <- paste0("PSR_Neojunctions_", current_date, ".tsv")
fwrite(dataframe_psr.neo, filename_psr.p1, sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)

###########################################################################
#  Step 3: Filter for Neojunctions: Count, Depth, Freq, Judge -------------
###########################################################################

list_tables <- list(
  dataframe_count,
  dataframe_depth,
  dataframe_freq,
  dataframe_judge
)
list_filenames <- list(
  paste0("Count_Neojunctions_", current_date, ".tsv"),
  paste0("Depth_Neojunctions_", current_date, ".tsv"),
  paste0("Freq_Neojunctions_", current_date, ".tsv"),
  paste0("Judge_Neojunctions_", current_date, ".tsv")
)

for (i in seq_along(list_tables)) {
  df_i <- list_tables[[i]] %>% semi_join(dataframe_psr.neo, by = "junc.id")
  fwrite(df_i, list_filenames[[i]], sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)
}

print("Neojunction extraction complete.")