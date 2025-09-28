#!/usr/bin/env Rscript
# Title: "Step 01: Tumor Purity Filtering ── SSNIP (Hartwig Edition) ───────────────────────────
# September 24, 2025 | Gaurav Raichand | The Insitute of Cancer Research
# Edits: Enhanced normalization for large Hartwig purity values (e.g., divide by 100,000 if >100),
#        flexible column renaming, and improved logging/error handling.

suppressPackageStartupMessages({
  library(data.table)   # fread for efficient reading
  library(tidyverse)    # dplyr / tibble for data manipulation
})

## ------------------------------------------------------------------------
## 0. Load purity table ---------------------------------------------------
## ------------------------------------------------------------------------
input_dir  <- Sys.getenv("INPUT_DIR")          # From config.sh
output_dir <- Sys.getenv("OUTPUT_DIR")
setwd(input_dir)

purity_file <- Sys.getenv("PURITY_FILE")
message("[INFO] Reading purity table: ", basename(purity_file))

# Read file based on extension (TXT or XLSX for Hartwig flexibility)
if (grepl("\\.txt$", purity_file, ignore.case = TRUE)) {
  purity_df <- fread(purity_file) |>
               as_tibble()
} else {
  library(readxl)                        # Only loaded if needed
  purity_df <- read_excel(purity_file) |>
               as_tibble()
}

## Flexible column renaming (tweak for Hartwig datasets) -------------------
# Assuming possible columns like "sample" or "tumor_purity" – adjust if needed
if (!"sample_id" %in% names(purity_df) && "sample" %in% names(purity_df)) {
  purity_df <- purity_df |> rename(sample_id = sample)
}
if (!"purity" %in% names(purity_df) && "tumor_purity" %in% names(purity_df)) {
  purity_df <- purity_df |> rename(purity = tumor_purity)
}

## Expect two columns: sample_id and purity (error if missing) -------------
if (!all(c("sample_id", "purity") %in% names(purity_df))) {
  stop("Purity file must contain columns 'sample_id' and 'purity' (after renaming). Check your Hartwig metadata.")
}

# Convert purity to numeric
purity_df <- purity_df |>
  mutate(purity = as.numeric(purity))

# Enhanced normalization for Hartwig large values (e.g., 60597 -> 0.606 if scaled by 1e5)
max_purity <- max(purity_df$purity, na.rm = TRUE)
if (max_purity > 100) {
  purity_df <- purity_df |> mutate(purity = purity / 100000)  # Adjust factor if not 1e5 (e.g., to get 0-1 scale)
  message("[INFO] Detected large purity values (max=", max_purity, "); normalized by dividing by 100,000.")
} else if (max_purity > 1) {
  purity_df <- purity_df |> mutate(purity = purity / 100)  # Handle percentages
  message("[INFO] Detected percentage-like values (max=", max_purity, "); normalized by dividing by 100.")
} else {
  message("[INFO] Purity values already in 0-1 scale (max=", max_purity, "); no normalization applied.")
}

## ------------------------------------------------------------------------
## 1. Filter by threshold -------------------------------------------------
## ------------------------------------------------------------------------
thresh <- as.numeric(Sys.getenv("MIN_PURITY", unset = "0.60"))
message("[INFO] Applying purity threshold ≥ ", thresh)

filtered_df <- purity_df |>
  filter(!is.na(purity), purity >= thresh)

message("[INFO] Kept ", nrow(filtered_df), " / ", nrow(purity_df), " samples")

## ------------------------------------------------------------------------
## 2. Write output --------------------------------------------------------
## ------------------------------------------------------------------------
setwd(output_dir)
out_file <- sprintf("Patient_List_Post_TumorPurity_Filter_%.2f.txt", thresh)
write.table(filtered_df, out_file,
            sep = "\t", quote = FALSE, row.names = FALSE)

message("[INFO] Written: ", out_file)
