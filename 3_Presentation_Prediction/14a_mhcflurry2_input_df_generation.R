#!/usr/bin/env Rscript
# Step 14: Generate input CSV file for running MHCFlurry 2.0
# September 27, 2025 | Gaurav Raichand | The Institute of Cancer Research
#     The input CSV file is expected to contain columns 
#     “allele”, “peptide”, and, optionally, “n_flank”, and “c_flank”.

# To run on command line: $ mhcflurry-predict INPUT.csv –out RESULT.csv

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(readxl)
library(tidyverse)
library(data.table)  # For faster reading and operations

# Load Directories -------------------------------------------------------
directory_12 <- Sys.getenv("STEP12_OUTPUT_DIR")
directory_14 <- Sys.getenv("STEP14_OUTPUT_DIR")

# Hardcode the exact input filenames from Step 12 (based on your ls output)
nmers_08_file <- "2023_0812_hlathenalist_msic_08mers.tsv"
nmers_09_file <- "2023_0812_hlathenalist_msic_09mers.tsv"
nmers_10_file <- "2023_0812_hlathenalist_msic_10mers.tsv"
nmers_11_file <- "2023_0812_hlathenalist_msic_11mers.tsv"

# Hardcode the date for output filenames (matching inputs)
run_date <- "2023_0812"

# Load Files -------------------------------------------------------------
setwd(directory_12)
nmers_08 <- fread(nmers_08_file, na = c("", "NA"))
nmers_09 <- fread(nmers_09_file, na = c("", "NA"))
nmers_10 <- fread(nmers_10_file, na = c("", "NA"))
nmers_11 <- fread(nmers_11_file, na = c("", "NA"))

###########################################################################
#  Step 1: Edit dataframes ------------------------------------------------
###########################################################################
# Remove the TPM columns (if present)
if ("TPM" %in% colnames(nmers_08)) nmers_08[, TPM := NULL]
if ("TPM" %in% colnames(nmers_09)) nmers_09[, TPM := NULL]
if ("TPM" %in% colnames(nmers_10)) nmers_10[, TPM := NULL]
if ("TPM" %in% colnames(nmers_11)) nmers_11[, TPM := NULL]

# Change the column names to the suitable ones needed to run MHCFLurry 2.0
setnames(nmers_08, old = colnames(nmers_08), new = c("peptide", "n_flank", "c_flank"))
setnames(nmers_09, old = colnames(nmers_09), new = c("peptide", "n_flank", "c_flank"))
setnames(nmers_10, old = colnames(nmers_10), new = c("peptide", "n_flank", "c_flank"))
setnames(nmers_11, old = colnames(nmers_11), new = c("peptide", "n_flank", "c_flank"))

# Remove all "-" from columns
nmers_08[, n_flank := gsub("-", "", as.character(n_flank))]
nmers_08[, c_flank := gsub("-", "", as.character(c_flank))]
nmers_09[, n_flank := gsub("-", "", as.character(n_flank))]
nmers_09[, c_flank := gsub("-", "", as.character(c_flank))]
nmers_10[, n_flank := gsub("-", "", as.character(n_flank))]
nmers_10[, c_flank := gsub("-", "", as.character(c_flank))]
nmers_11[, n_flank := gsub("-", "", as.character(n_flank))]
nmers_11[, c_flank := gsub("-", "", as.character(c_flank))]

# Define alleles as tibble
allele_a0101 <- tibble(allele = "HLA-A0101")
allele_a0201 <- tibble(allele = "HLA-A0201")
allele_a0301 <- tibble(allele = "HLA-A0301")
allele_a1101 <- tibble(allele = "HLA-A1101")
allele_a2402 <- tibble(allele = "HLA-A2402")

# For each nmer dataframe, replicate for each allele and combine
for (i in 1:4) {
  if (i == 1) {nmer_i <- nmers_08}
  if (i == 2) {nmer_i <- nmers_09}
  if (i == 3) {nmer_i <- nmers_10}
  if (i == 4) {nmer_i <- nmers_11}
  
  nmer_i_a0101 <- cbind(allele_a0101, nmer_i)
  nmer_i_a0201 <- cbind(allele_a0201, nmer_i)
  nmer_i_a0301 <- cbind(allele_a0301, nmer_i)
  nmer_i_a1101 <- cbind(allele_a1101, nmer_i)
  nmer_i_a2402 <- cbind(allele_a2402, nmer_i)
  
  n_mer_i_all <- rbind(nmer_i_a0101, nmer_i_a0201, nmer_i_a0301, nmer_i_a1101, nmer_i_a2402)
  
  if (i == 1) {nmers_08_final <- n_mer_i_all}
  if (i == 2) {nmers_09_final <- n_mer_i_all}
  if (i == 3) {nmers_10_final <- n_mer_i_all}
  if (i == 4) {nmers_11_final <- n_mer_i_all}
}

setwd(directory_14)
fwrite(nmers_08_final, paste0("08mer_mhcflurry_input_", run_date, ".csv"), sep = ",", na = "NA", col.names = TRUE, quote = TRUE)
fwrite(nmers_09_final, paste0("09mer_mhcflurry_input_", run_date, ".csv"), sep = ",", na = "NA", col.names = TRUE, quote = TRUE)
fwrite(nmers_10_final, paste0("10mer_mhcflurry_input_", run_date, ".csv"), sep = ",", na = "NA", col.names = TRUE, quote = TRUE)
fwrite(nmers_11_final, paste0("11mer_mhcflurry_input_", run_date, ".csv"), sep = ",", na = "NA", col.names = TRUE, quote = TRUE)

print("Step 14: MHCflurry input generation complete")
