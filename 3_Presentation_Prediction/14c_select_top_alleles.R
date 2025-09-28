#!/usr/bin/env Rscript
# Step 14c: Select TOP HLA-allele presentation score for each n-mer (MHCFlurry 2.0)
# September 27, 2025 | Gaurav Raichand | The Institute of Cancer Research

# Purpose: For each of the n-mers, identify the best HLA-allele that it binds to the n-mer.

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(tidyverse)
library(data.table)  # For faster reading and operations
library(parallel)    # For detectCores
library(doParallel)
library(foreach)

#  Load Directories -------------------------------------------------------
directory_14_in <- Sys.getenv("OUTPUT_DIR")
directory_14_out <- Sys.getenv("STEP14_OUTPUT_DIR")

# Hardcode the run date to match input file naming (avoids current date mismatch)
RUN_DATE <- "2023_0812"  # Matches your MHCflurry output files

# Load Files --------------------------------------------------------------
setwd(directory_14_in)
nmer_08 <- fread(paste0("08mers_flank_mhcflurry_", RUN_DATE, ".csv"), na = c("", "NA"))
nmer_09 <- fread(paste0("09mers_flank_mhcflurry_", RUN_DATE, ".csv"), na = c("", "NA"))
nmer_10 <- fread(paste0("10mers_flank_mhcflurry_", RUN_DATE, ".csv"), na = c("", "NA"))
nmer_11 <- fread(paste0("11mers_flank_mhcflurry_", RUN_DATE, ".csv"), na = c("", "NA"))

###########################################################################
#  Step 1: Select highest HLA-allele and associated score -----------------
###########################################################################

# Replace NA with "" (vectorized for efficiency)
nmer_08[is.na(nmer_08)] <- ""
nmer_09[is.na(nmer_09)] <- ""
nmer_10[is.na(nmer_10)] <- ""
nmer_11[is.na(nmer_11)] <- ""

# Set up parallel processing (use SLURM cores if available)
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores() - 1))
num_cores <- max(1, num_cores)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Function to process one n-mer dataset in parallel
process_nmer <- function(dt) {
  # Split into chunks for parallel processing
  chunk_size <- ceiling(nrow(dt) / num_cores)
  indices <- split(1:nrow(dt), ceiling(seq_along(1:nrow(dt)) / chunk_size))
  
  # Parallel foreach over chunks
  result_list <- foreach(idx = indices, .combine = rbind, .packages = c("data.table")) %dopar% {
    dt_chunk <- dt[idx, ]
    dt_chunk[, .SD[which.max(mhcflurry_presentation_score)], by = .(peptide, n_flank, c_flank)]
  }
  
  result_list
}

# Process each n-mer in parallel
nmer_08_edit <- process_nmer(nmer_08)
nmer_09_edit <- process_nmer(nmer_09)
nmer_10_edit <- process_nmer(nmer_10)
nmer_11_edit <- process_nmer(nmer_11)

# Stop the cluster
stopCluster(cl)

# Export Files ------------------------------------------------------------
setwd(directory_14_out)
current_date <- format(Sys.Date(), "%Y%m%d")  # Use today's date for output naming
fwrite(nmer_08_edit, paste0("mhcflurry_08mer_selected_alleles_", current_date, ".tsv"), sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)
fwrite(nmer_09_edit, paste0("mhcflurry_09mer_selected_alleles_", current_date, ".tsv"), sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)
fwrite(nmer_10_edit, paste0("mhcflurry_10mer_selected_alleles_", current_date, ".tsv"), sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)
fwrite(nmer_11_edit, paste0("mhcflurry_11mer_selected_alleles_", current_date, ".tsv"), sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)

print("Step 14c completed: Selected top HLA-alleles for each n-mer and saved output files.")