#!/usr/bin/env Rscript
# Title: "Step 7: Count, Depth, Frequency, and Judgement Tables for Retained Splicing Junctions"
# September 24, 2025 | Gaurav Raichand | The Institute of Cancer Research

# Purpose: Generate and prepare the following tables for the retained SPLICING JUNCTIONS data:
#           1. Count table
#           2. Depth table
#           3. Frequency table
#           4. Judgement table

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(tidyverse)
library(readxl)
library(ggsci)
library(data.table)
library(foreach)
library(doParallel)
library(conflicted)  # Handle conflicts

# Set conflict preferences
conflict_prefer("filter", "dplyr")

# Initialize log file for progress (clears it first)
log_file <- file.path(Sys.getenv("OUTPUT_DIR"), "step7_progress.log")
cat("", file = log_file)  # Clear/Initialize the file

# Establish Directories ---------------------------------------------------

directory_out.step01 <- Sys.getenv("STEP01_OUTPUT_DIR", Sys.getenv("OUTPUT_DIR"))
directory_out.step06 <- Sys.getenv("STEP06_OUTPUT_DIR", Sys.getenv("OUTPUT_DIR"))
directory_out        <- Sys.getenv("OUTPUT_DIR")

# Load Files --------------------------------------------------------------

# From Step 1: Patient/Sample List (Tumor Purity >= MIN_PURITY) - fixed formatting to "0.60"
current_date  <- format(Sys.Date(), "%Y%m%d")
purity_thresh <- as.numeric(Sys.getenv("MIN_PURITY", unset = "0.60"))
filename_tcga <- paste0("Patient_List_Post_TumorPurity_Filter_", format(purity_thresh, nsmall = 2), ".txt")
full_path_tcga <- file.path(directory_out.step01, filename_tcga)
if (!file.exists(full_path_tcga)) {
  stop("[ERROR] File does not exist: ", full_path_tcga, ". Check directory_out.step01 and if Step 1 ran successfully.")
}
dataframe_tcga <- fread(full_path_tcga)  # sample_id + purity expected

# From Step 6: SJ Overlap Table
filename_overlap.table <- paste0("SJ_Overlap_Table_", current_date, ".tsv")
full_path_overlap <- file.path(directory_out.step06, filename_overlap.table)
if (!file.exists(full_path_overlap)) {
  stop("[ERROR] File does not exist: ", full_path_overlap, ". Check directory_out.step06 and if Step 6 ran successfully.")
}
dataframe_overlap.table <- fread(full_path_overlap)  # As data.table

# From Step 6: SJ Count Table
filename_count.table <- paste0("SJ_Count_Table_SJ.out.tab_Retained_", current_date, ".tsv")
full_path_count <- file.path(directory_out.step06, filename_count.table)
if (!file.exists(full_path_count)) {
  stop("[ERROR] File does not exist: ", full_path_count, ". Check directory_out.step06 and if Step 6 ran successfully.")
}
dataframe_count.table <- fread(full_path_count)  # As data.table

###########################################################################
#  Step 1: Obtain the Count | Depth | Freq | Judge Tables -----------------
###########################################################################

start.time <- proc.time()

# Set up parallel backend (adjust n_cores to your SLURM cpus-per-task)
n_cores <- detectCores() - 1  # Leave one core free; adjust as needed
registerDoParallel(cores = n_cores)

# Split the overlap table into chunks for parallel processing
n_rows <- nrow(dataframe_overlap.table)
chunk_size <- ceiling(n_rows / n_cores)
chunks <- split(1:n_rows, rep(1:n_cores, each = chunk_size, length.out = n_rows))

# Parallel loop to compute res.i for each chunk
res_list <- foreach(chunk = chunks, .combine = c, .packages = c("data.table", "tidyr")) %dopar% {
  local_res <- list()
  for (i in chunk) {
    JUNC.ID <- dataframe_overlap.table[i, junc.id]
    JUNC.ID.OVERLAP <- dataframe_overlap.table[i, junc.id.overlap]
    
    if (is.na(JUNC.ID.OVERLAP)) {
      res.i <- dataframe_count.table[junc.id == JUNC.ID] %>%
        pivot_longer(cols = -junc.id, names_to = "case", values_to = "value") %>%
        pivot_wider(names_from = junc.id, values_from = value) %>%
        rename(count = !!JUNC.ID) %>%
        mutate(count = ifelse(is.na(count), 0, count),
               depth = count,
               freq  = round(count / depth, 4),
               freq  = ifelse(is.na(freq), 0, freq),
               judge = ifelse(count >= 10 & depth >= 20 & freq >= 0.01, 1, 0),
               junc.id = JUNC.ID) %>%
        select(junc.id, case, count, depth, freq, judge)
    } else {
      res.i <- rbind(
        dataframe_count.table[junc.id == JUNC.ID],
        dataframe_count.table[junc.id == JUNC.ID.OVERLAP]) %>%
        pivot_longer(cols = -junc.id, names_to = "case", values_to = "value") %>%
        pivot_wider(names_from = junc.id, values_from = value) %>%
        rename(count = !!JUNC.ID, count.o = !!JUNC.ID.OVERLAP) %>%
        mutate(count  = ifelse(is.na(count), 0, count),
               count.o = ifelse(is.na(count.o), 0, count.o),
               depth   = count + count.o,
               freq    = round(count / depth, 4),
               freq    = ifelse(is.na(freq), 0, freq),
               judge   = ifelse(count >= 10 & depth >= 20 & freq >= 0.01, 1, 0),
               junc.id = JUNC.ID) %>%
        select(junc.id, case, count, depth, freq, judge)
    }
    local_res[[length(local_res) + 1]] <- res.i
  }
  
  # Log chunk completion (appends to file from worker)
  cat(paste0("[PARALLEL] ", Sys.time(), " - Completed chunk ending at row ", max(chunk), "\n"), file = log_file, append = TRUE)
  
  local_res
}

# Log after parallel section
cat(paste0("[MAIN] ", Sys.time(), " - Parallel processing complete. Combining results...\n"), file = log_file, append = TRUE)

# Now combine the list of res.i into final tables (batch pivot)
res_all <- rbindlist(res_list)

# **KEY FIX**: Aggregate duplicates before pivoting to ensure unique (junc.id, case) combinations
cat(paste0("[MAIN] ", Sys.time(), " - Aggregating duplicates before pivoting...\n"), file = log_file, append = TRUE)
res_all <- res_all %>%
  group_by(junc.id, case) %>%
  summarise(
    count = sum(count, na.rm = TRUE),
    depth = sum(depth, na.rm = TRUE), 
    freq  = mean(freq, na.rm = TRUE),
    judge = as.integer(mean(judge, na.rm = TRUE) >= 0.5),  # majority vote for binary judge
    .groups = 'drop'
  )

count <- res_all %>% select(junc.id, case, count) %>% pivot_wider(names_from = case, values_from = count)
depth <- res_all %>% select(junc.id, case, depth) %>% pivot_wider(names_from = case, values_from = depth)
freq  <- res_all %>% select(junc.id, case, freq ) %>% pivot_wider(names_from = case, values_from = freq)
judge <- res_all %>% select(junc.id, case, judge) %>% pivot_wider(names_from = case, values_from = judge)

RUNTIME <- proc.time() - start.time
print(RUNTIME)

# Export files without PSR filter
filename_raw_count <- paste0("Raw_Count_Table_", current_date, ".tsv")
full_path_raw_count <- file.path(directory_out, filename_raw_count)
fwrite(count, full_path_raw_count, sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)

filename_raw_depth <- paste0("Raw_Depth_Table_", current_date, ".tsv")
full_path_raw_depth <- file.path(directory_out, filename_raw_depth)
fwrite(depth, full_path_raw_depth, sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)

filename_raw_freq <- paste0("Raw_Freq_Table_", current_date, ".tsv")
full_path_raw_freq <- file.path(directory_out, filename_raw_freq)
fwrite(freq, full_path_raw_freq, sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)

filename_raw_judge <- paste0("Raw_Judge_Table_", current_date, ".tsv")
full_path_raw_judge <- file.path(directory_out, filename_raw_judge)
fwrite(judge, full_path_raw_judge, sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)

###########################################################################
#  Step 2: Calculate Judgement and Filter out Junctions Based on PSR > 0.10
###########################################################################

# --- ID normalization and mapping to judge columns ---
# Judge columns can be altered by name repair (e.g., trailing dots); normalize both sides and map back to actual judge names
judge_cols      <- setdiff(colnames(judge), "junc.id")
judge_cols_norm <- sub("\\.+$", "", trimws(judge_cols))  # strip trailing dots

# Choose patient list column: prefer "case" if present, else "sample_id"
case_col <- if ("case" %in% names(dataframe_tcga)) "case" else if ("sample_id" %in% names(dataframe_tcga)) "sample_id" else stop("[ERROR] Neither case nor sample_id present in patient list.")

patient_ids      <- trimws(dataframe_tcga[[case_col]])
patient_ids_norm <- sub("\\.+$", "", patient_ids)

# Overlap on normalized space, then map back to exact judge column names
overlap_norm <- intersect(unique(patient_ids_norm), unique(judge_cols_norm))
if (length(overlap_norm) == 0) {
  stop("[ERROR] No overlap between patient list IDs and judge columns after normalization; check inputs.")
}
map_norm_to_judge <- setNames(judge_cols, judge_cols_norm)
CASES_JUDGE <- unique(unname(map_norm_to_judge[overlap_norm]))

message("[DEBUG] Using patient column: ", case_col,
        " | Overlap (normalized): ", length(overlap_norm),
        " | CASES used (judge columns): ", length(CASES_JUDGE))
# --- end ID normalization block ---

# 2a. Prepare a list of all disease subtypes (General list - adapted to Hartwig "all")
cases <- list(CASES_JUDGE)
names(cases) <- c("all")

# Depending on which project is being studied - simplified to one project for Hartwig
list_project   <- list(cases)
list_judge.all <- list(judge.all_cases = NULL)

# 2b. Calculate the positive sample rates for all of the disease subtypes
for (h in 1:length(list_project)) {
  project <- list_project[[h]]
  for (i in 1:length(project)) {
    # Status bar
    message("[PROGRESS] Step 2: Group ", h, ", Subgroup ", i)
    
    # Select the disease type for each iteration
    CASES <- project[[i]]
    judge.i <- judge %>%
      select(junc.id, dplyr::all_of(CASES))  # strict selection of actual judge columns
    
    # Debug: show selection
    message("[DEBUG] Number of CASES: ", length(CASES))
    message("[DEBUG] First 5 CASES: ", paste(head(CASES, 5), collapse = ", "))
    message("[DEBUG] Columns in judge.i: ", paste(colnames(judge.i), collapse = ", "))
    
    # Calculate the positive sample rate (PSR) with NA/zero handling
    judge.i <- judge.i %>%
      mutate(total_positive = apply(select(., -junc.id), 1, sum, na.rm = TRUE)) %>%
      mutate(psr = ifelse(length(CASES) == 0 | total_positive == 0, 0,
                          round(total_positive / length(CASES), 4))) %>%
      select(junc.id, psr)
    
    colnames(judge.i)[2] <- paste0(names(cases)[[i]], "_n", length(CASES))
    
    if (i == 1) {
      list_judge.all[[h]] <- judge.i
    } else {
      list_judge.all[[h]] <- list_judge.all[[h]] %>%
        full_join(judge.i, by = "junc.id")
    }
  }
}

judge.all_cases <- list_judge.all[[1]] %>%
  mutate(max = apply(select(., -junc.id), 1, max)) %>%
  filter(max >= as.numeric(Sys.getenv("MIN_SAMPLES_PCT", unset = "0.10")))

# This added variable contains all NJs regardless of PSR threshold
judge.all_cases_allnjs <- list_judge.all[[1]]
judge.all_cases_allnjs <- judge.all_cases_allnjs %>% mutate(max = apply(select(., -junc.id), 1, max))
filename_psr_all <- paste0("PSR_All_NJs_", current_date, ".tsv")
full_path_psr_all <- file.path(directory_out, filename_psr_all)
write_tsv(judge.all_cases_allnjs, full_path_psr_all, na = "NA", col_names = TRUE, quote_escape = "double")

###########################################################################
#  Step 3: Generate the Final Count, Depth, Freq, Judge, and Overlap Tables
###########################################################################

# Generate the Tables (simplified to one project for Hartwig)
dataframe_count.pass  <- count %>%  semi_join(judge.all_cases, by = "junc.id")
dataframe_depth.pass  <- depth %>%  semi_join(judge.all_cases, by = "junc.id")
dataframe_freq.pass   <- freq  %>%  semi_join(judge.all_cases, by = "junc.id")
dataframe_judge.pass  <- judge %>%  semi_join(judge.all_cases, by = "junc.id")
dataframe_overlap.table.pass <- dataframe_overlap.table %>% semi_join(judge.all_cases, by = "junc.id")

# Check for NA values
print(anyNA(dataframe_count.pass))
print(anyNA(dataframe_depth.pass))
print(anyNA(dataframe_freq.pass))
print(anyNA(dataframe_judge.pass))

###########################################################################
#  Step 4: Output Files ---------------------------------------------------
###########################################################################

# Export Tables
filename_count <- paste0("Count_Table_Retained_and_Passed_Junctions_", current_date, ".tsv")
full_path_count <- file.path(directory_out, filename_count)
write_tsv(dataframe_count.pass, full_path_count, na = "NA", col_names = TRUE, quote_escape = "double")

filename_depth <- paste0("Depth_Table_Retained_and_Passed_Junctions_", current_date, ".tsv")
full_path_depth <- file.path(directory_out, filename_depth)
write_tsv(dataframe_depth.pass, full_path_depth, na = "NA", col_names = TRUE, quote_escape = "double")

filename_freq <- paste0("Freq_Table_Retained_and_Passed_Junctions_", current_date, ".tsv")
full_path_freq <- file.path(directory_out, filename_freq)
write_tsv(dataframe_freq.pass, full_path_freq, na = "NA", col_names = TRUE, quote_escape = "double")

filename_judge <- paste0("Judgement_Table_Retained_and_Passed_Junctions_", current_date, ".tsv")
full_path_judge <- file.path(directory_out, filename_judge)
write_tsv(dataframe_judge.pass, full_path_judge, na = "NA", col_names = TRUE, quote_escape = "double")

filename_judge_all <- paste0("PSR_Table_", current_date, ".tsv")
full_path_judge_all <- file.path(directory_out, filename_judge_all)
write_tsv(judge.all_cases, full_path_judge_all, na = "NA", col_names = TRUE, quote_escape = "double")

filename_overlap <- paste0("Overlap_Table_Retained_and_Passed_Junctions_", current_date, ".tsv")
full_path_overlap <- file.path(directory_out, filename_overlap)
write_tsv(dataframe_overlap.table.pass, full_path_overlap, na = "NA", col_names = TRUE, quote_escape = "double")
