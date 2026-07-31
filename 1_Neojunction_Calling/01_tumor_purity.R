#!/usr/bin/env Rscript
# Title: "Step 01: Tumor Purity Filtering ── SSNIP (Hartwig Edition) ───────────────────────────
# October 26, 2025 | Gaurav Raichand | The Insitute of Cancer Research
# Edits: ilters purity >0.60, outputs high-purity list. Loads samples.txt (307 IDs), matches to metadata,

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(stringr)
})

## ------------------------------------------------------------------------
## 0. Setup ---------------------------------------------------------------
## ------------------------------------------------------------------------
base_dir <- getwd()
input_dir   <- Sys.getenv("INPUT_DIR", unset = file.path(base_dir, "0_Input_Files"))
output_dir  <- Sys.getenv("OUTPUT_DIR", unset = file.path(base_dir, "1_Outputs"))
samples_file <- Sys.getenv("SAMPLES_FILE", unset = file.path(input_dir, "samples.txt"))
metadata_file <- Sys.getenv("METADATA_FILE", unset = file.path(input_dir, "mapped_sample_purity.txt"))
thresh <- as.numeric(Sys.getenv("MIN_PURITY", unset = "0.60"))

# File checks
if (!file.exists(samples_file)) stop("[ERROR] samples.txt not found: ", samples_file)
if (!file.exists(metadata_file)) stop("[ERROR] Metadata/purity file not found: ", metadata_file)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
message("[INFO] Samples: ", samples_file)
message("[INFO] Purity file: ", metadata_file)
message("[INFO] Threshold: >= ", thresh)

## ------------------------------------------------------------------------
## 1. Load Data -----------------------------------------------------------
## ------------------------------------------------------------------------
samples_list <- readLines(samples_file)
message("[INFO] Loaded ", length(samples_list), " samples ^_^")

raw_metadata <- fread(metadata_file) |> as_tibble()

# --- Detect file format and normalise to sample_id / tumorPurity ----------
# Format A (KURA):    columns = sample_id, purity
# Format B (Hartwig): columns = sampleId,  tumorPurity  (+ many more columns)
if ("sample_id" %in% colnames(raw_metadata) && "purity" %in% colnames(raw_metadata)) {
  message("[INFO] Detected KURA-style purity file (sample_id / purity columns).")
  metadata <- raw_metadata |>
    rename(sampleId = sample_id, tumorPurity = purity) |>
    mutate(tumorPurity = as.numeric(tumorPurity)) |>
    filter(!is.na(tumorPurity))
} else if ("sampleId" %in% colnames(raw_metadata) && "tumorPurity" %in% colnames(raw_metadata)) {
  message("[INFO] Detected Hartwig-style metadata file (sampleId / tumorPurity columns).")
  metadata <- raw_metadata |>
    mutate(tumorPurity = as.numeric(tumorPurity)) |>
    filter(!is.na(tumorPurity))
} else {
  stop("[ERROR] Purity file must have either (sample_id, purity) or (sampleId, tumorPurity) columns. ",
       "Found: ", paste(colnames(raw_metadata), collapse = ", "))
}

message("[INFO] Purity file: ", nrow(metadata), " rows; purity range: ",
        round(min(metadata$tumorPurity), 2), "-", round(max(metadata$tumorPurity), 2))

## ------------------------------------------------------------------------
## 2. Match samples to purity file ----------------------------------------
## ------------------------------------------------------------------------
# KURA: sample IDs are already identical between samples.txt and purity file
# — direct join. Hartwig: fall back to substring regex matching for safety.
matched_direct <- metadata |>
  filter(sampleId %in% samples_list) |>
  distinct(sampleId, .keep_all = TRUE)

if (nrow(matched_direct) >= length(samples_list) * 0.5) {
  matched_metadata <- matched_direct
  message("[INFO] Direct ID match: ", nrow(matched_metadata), "/", length(samples_list), " samples.")
} else {
  message("[INFO] Direct match low (", nrow(matched_direct), "); falling back to regex substring match...")
  samples_escaped <- str_replace_all(samples_list, "([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1")
  samples_pattern <- paste(samples_escaped, collapse = "|")
  matched_metadata <- metadata |>
    filter(str_detect(sampleId, regex(samples_pattern, ignore_case = TRUE))) |>
    distinct(sampleId, .keep_all = TRUE)
  message("[INFO] Regex match: ", nrow(matched_metadata), "/", length(samples_list), " samples.")
}

n_matched <- nrow(matched_metadata)
if (n_matched < length(samples_list) * 0.5) {
  message("[WARN] Low match rate (<50%). Debug info:")
  unmatched <- setdiff(samples_list, matched_metadata$sampleId)
  cat("First 5 unmatched samples:\n", paste(head(unmatched, 5), collapse = "\n"), "\n")
  cat("First 5 purity file IDs:\n", paste(head(metadata$sampleId, 5), collapse = "\n"), "\n")
}

# Keep only the two columns all downstream steps need; extra columns are optional
matched_metadata <- matched_metadata |>
  select(sampleId, tumorPurity, any_of(c("primaryTumorType_real", "primaryTumorLocation_real", "biopsy_site_simple")))

## ------------------------------------------------------------------------
## 3. Filter High-Purity --------------------------------------------------
## ------------------------------------------------------------------------
high_purity <- matched_metadata |>
  filter(tumorPurity >= thresh)

# Normalization (purity 0.2-1 range looks good; no change needed)
max_p <- max(high_purity$tumorPurity, na.rm = TRUE)
if (max_p > 1) {
  high_purity <- high_purity |> mutate(tumorPurity = tumorPurity / 100)
  message("[INFO] Normalized purity by /100 (max: ", max_p, ")")
} else {
  message("[INFO] Purity OK (max: ", max_p, "); no normalization.")
}

n_filtered <- nrow(high_purity)
message("[INFO] Filtered: ", n_filtered, " high-purity samples (>= ", thresh, "; ~", 
        round((n_filtered / n_matched) * 100, 1), "% of matched)")

if (n_filtered == 0) {
  stop("[ERROR] No high-purity samples found. Lower threshold or check data.")
}

# Stats
cat("[INFO] Purity summary (filtered):\n")
print(summary(high_purity$tumorPurity))

## ------------------------------------------------------------------------
## 4. Outputs (Modified for Pipeline Compatibility) -----------------------
## ------------------------------------------------------------------------
setwd(output_dir)

# Output 1: Pipeline-required format (original naming/structure)
# Two columns: sample_id (extracted short ID if needed) and purity
pipeline_output <- high_purity %>%
  mutate(sample_id = sampleId) %>%  # Use full sampleId; adjust if short needed (see note below)
  select(sample_id, purity = tumorPurity)  # Rename to match pipeline expectations

pipeline_file <- sprintf("Patient_List_Post_TumorPurity_Filter_%.2f.txt", thresh)
write.table(pipeline_output, pipeline_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("[INFO] Pipeline output: ", pipeline_file, " (", nrow(pipeline_output), " rows; 2 cols: sample_id, purity)")

# Output 2: Detailed metadata (optional reference; doesn't break pipeline)
detailed_file <- sprintf("high_purity_metadata_%.2f.tsv", thresh)
write.table(high_purity, detailed_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("[INFO] Detailed metadata: ", detailed_file, " (full columns)")

# Output 3: Sample IDs list (for SLURM/IRFinder loops; one ID per line)
filtered_ids <- high_purity$sampleId
ids_file <- sprintf("high_purity_sample_ids_%.2f.txt", thresh)
writeLines(filtered_ids, ids_file)
message("[INFO] Sample IDs list: ", ids_file, " (", length(filtered_ids), " lines)")

# Peek
cat("[INFO] First 5 pipeline rows (sample_id, purity):\n")
print(head(pipeline_output, 5))

message("[INFO] Script complete! Files ready for pipeline ^_^")
