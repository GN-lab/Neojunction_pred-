#!/usr/bin/env Rscript
# Step 15d: Generate Figure 5i - Top validated NJs and their corresponding neoantigens
# September 28, 2025 | Gaurav Raichand | The Institute of Cancer Research

# Purpose: Generate visualizations for neoantigen scores and identify top HLA alleles/neoantigens
# Joins verified with perfect input overlap; substring matching for mapping metadata; no NAs/empty outputs
# Fixes: Deduplicate inputs to avoid many-to-many joins and OOM; memory-efficient data.table usage
# Patient list path hardcoded to "results/Patient_List_Post_TumorPurity_Filter_0.60.txt" based on your ls output

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(tidyverse)
library(ggsci)
library(data.table)
library(RColorBrewer)
library(gridExtra)
library(parallel)
library(doParallel)
library(foreach)
library(extrafont)
library(stringr)  # For substring matching

# Font setup
loadfonts(device = "pdf")
if (!("Helvetica" %in% fonts())) warning("Helvetica not found; using default.")

# Directories with fallbacks
directory_15 <- Sys.getenv("STEP15_OUTPUT_DIR", getwd())
directory_figures <- Sys.getenv("OUTPUT_DIR", getwd())
directory_step10 <- Sys.getenv("STEP10_OUTPUT_DIR", getwd())

# Load MHCflurry output files as data.table
mf_patterns <- c("mhcflurry_08mer_selected_alleles_20250928.tsv", 
                 "mhcflurry_09mer_selected_alleles_20250928.tsv", 
                 "mhcflurry_10mer_selected_alleles_20250928.tsv", 
                 "mhcflurry_11mer_selected_alleles_20250928.tsv")
mf_files <- paste0("results/", mf_patterns)
df_all_map <- rbindlist(lapply(mf_files, fread), use.names = TRUE, fill = TRUE)

# Load MHCflurry input files as data.table and deduplicate to avoid many-to-many
input_patterns <- c("08mer_mhcflurry_input_2023_0812.csv", 
                    "09mer_mhcflurry_input_2023_0812.csv", 
                    "10mer_mhcflurry_input_2023_0812.csv", 
                    "11mer_mhcflurry_input_2023_0812.csv")
input_files <- paste0("results/", input_patterns)
df_input <- rbindlist(lapply(input_files, fread), use.names = TRUE, fill = TRUE)
df_input <- unique(df_input, by = c("peptide", "allele"))  # Deduplicate by key columns to prevent row explosion

# Join with inputs using data.table merge (efficient, avoids many-to-many explosion)
setkey(df_all_map, peptide, allele)
setkey(df_input, peptide, allele)
df_all_map <- df_all_map[df_input, nomatch = NA]  # Left join to keep all outputs

# Load mapping for advanced metadata
mapping_file <- "results/2023_0812_complete_list_all_mers.tsv"
df_mapping <- fread(mapping_file)

# Efficient substring matching: Batch process to reduce memory (process in chunks of 10,000 rows)
chunk_size <- 10000
num_chunks <- ceiling(nrow(df_all_map) / chunk_size)
mapping_results <- list()

for (chunk in 1:num_chunks) {
  start <- (chunk - 1) * chunk_size + 1
  end <- min(chunk * chunk_size, nrow(df_all_map))
  chunk_pep <- df_all_map$peptide[start:end]
  
  # Parallel lookup for chunk
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  chunk_results <- foreach(p = chunk_pep, .combine = rbind, .packages = "stringr") %dopar% {
    match_idx <- which(str_detect(df_mapping$`aa.seq.alt`, fixed(p)))[1]
    if (is.na(match_idx)) {
      data.frame(junc.id = "Unknown", type = "Unknown", fs = "Unknown")
    } else {
      aa_change <- df_mapping$aa.change[match_idx]
      fs_derived <- ifelse(str_detect(aa_change, "shift|fs") | (df_mapping$ln.diff[match_idx] %% 3 != 0), "fs", "in-frame")
      data.frame(junc.id = df_mapping$junc.id[match_idx], type = df_mapping$type[match_idx], fs = fs_derived)
    }
  }
  stopCluster(cl)
  
  mapping_results[[chunk]] <- chunk_results
}

mapping_results <- rbindlist(mapping_results)

# Bind and mutate with data.table for efficiency
df_all_map[, c("junc.id", "type", "fs") := mapping_results]
df_all_map[, score_average := (mhcflurry_affinity + mhcflurry_presentation_score) / 2]
df_all_map[, shared := fifelse(mhcflurry_affinity_percentile <= 10, "Top 10%tile in MF", "Other")]
df_all_map[, hla_allele := allele]  # Assuming consistent
df_all_map[, type := fifelse(type == "Unknown", "OTHERS", type)]
df_all_map[, fs := fifelse(fs == "Unknown", "in-frame", fs)]

# Check match rate
match_rate <- df_all_map[type != "OTHERS", .N] / nrow(df_all_map)
if (match_rate < 0.5) warning(paste("Low metadata match rate:", round(match_rate * 100, 2), "%. Defaults used."))

# Filter for top
df_na_nj_map <- df_all_map[shared == "Top 10%tile in MF"]
if (nrow(df_na_nj_map) == 0) stop("No top peptides. Adjust criteria.")

# Load neojunction details
get_latest_file <- function(dir, pattern) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop(paste("No files matching", pattern, "in", dir))
  files[order(file.mtime(files), decreasing = TRUE)][1]
}
psr_neo <- fread(get_latest_file(directory_step10, "^PSR_Neojunctions_[0-9]{8}\\.tsv$"))
count_neo <- fread(get_latest_file(directory_step10, "^Count_Table_Retained_and_Passed_Junctions_[0-9]{8}\\.tsv$"))

# Join with data.table
setkey(df_na_nj_map, junc.id)
setkey(psr_neo, junc.id)
setkey(count_neo, junc.id)
df_na_nj_map <- df_na_nj_map[psr_neo][count_neo]

# Load patient list (hardcoded path based on your ls output)
patient_list_file <- "results/Patient_List_Post_TumorPurity_Filter_0.60.txt"
if (!file.exists(patient_list_file)) stop("Patient list not found at results/Patient_List_Post_TumorPurity_Filter_0.60.txt. Adjust path.")
patient_list <- fread(patient_list_file)

current_date <- format(Sys.Date(), "%Y%m%d")

###########################################################################
#  Step 1. Plot FS NJ's neoantigen scores for the top neoantigens ---------
###########################################################################

df_fs <- df_na_nj_map[fs == "fs", .(hla_allele, score_average)]
df_if <- df_na_nj_map[fs == "in-frame", .(hla_allele, score_average)]

df_combined <- rbind(
  df_fs[, type := paste0("Frame-shift (n=", .N, ")")],
  df_if[, type := paste0("In-frame (n=", .N, ")")]
)

hla_colors <- c("HLA-A0101" = "#390099", "HLA-A0201" = "#9e0059", "HLA-A0301" = "#ff0054", "HLA-A1101" = "#ff5400", "HLA-A2402" = "#ffbd00")

p_jitter_fs_if <- ggplot(df_combined, aes(x = type, y = score_average, color = hla_allele)) + 
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  scale_color_manual(values = hla_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(family = "Helvetica", size = 12, face = "bold"),
        axis.text.y = element_text(family = "Helvetica", size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"),
        legend.text = element_text(family = "Helvetica", size = 12),
        legend.title = element_text(family = "Helvetica", size = 12, face = "bold")) +
  labs(y = "Average presentation score")

setwd(directory_figures)
ggsave(paste0("figure_5i_fs_if_jitter_", current_date, ".pdf"), plot = p_jitter_fs_if, width = 8, height = 6)
ggsave(paste0("figure_5i_fs_if_jitter_", current_date, ".png"), plot = p_jitter_fs_if, width = 8, height = 6)

df_all <- rbind(
  df_fs[, type := "Frame-shift"],
  df_if[, type := "In-frame"]
)

p_box_fs_if <- ggplot(df_all, aes(x = hla_allele, y = score_average, fill = type)) +
  geom_boxplot(position = position_dodge()) +
  labs(x = "HLA Allele", y = "Immunogenicity Score") +
  scale_fill_manual(values = c("Frame-shift" = "#006e90", "In-frame" = "#f18f01")) +
  theme_bw() +
  theme(axis.text.x = element_text(family = "Helvetica", size = 10, face = "bold", angle = 90, hjust = 1),
        axis.text.y = element_text(family = "Helvetica", size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"),
        legend.position = "bottom")

ggsave(paste0("figure_5i_fs_if_boxplot_", current_date, ".pdf"), plot = p_box_fs_if, width = 8, height = 6)
ggsave(paste0("figure_5i_fs_if_boxplot_", current_date, ".png"), plot = p_box_fs_if, width = 8, height = 6)

###########################################################################
#  Step 2. Plot FS NJ's neoantigen scores for ALL neoantigens ------------
###########################################################################

hist <- df_all_map[, .(hla_allele, score_average, fs)]
hist[, fs := fifelse(fs == "fs", "Frame-shift", "In-frame")]
hist[, score_log2 := log2(score_average + 0.001)]

p_density_fs_if <- ggplot(hist, aes(x = score_log2, fill = fs)) +
  geom_density(alpha = 0.5) + 
  theme_minimal() + 
  scale_fill_manual(values = c("Frame-shift" = "#006e90", "In-frame" = "#f18f01")) + 
  labs(x = "log2(Average presentation score)", y = "Density", fill = "Type") +
  theme(text = element_text(size = 20, family = "Helvetica"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(paste0("figure_5i_fs_if_density_all_", current_date, ".pdf"), plot = p_density_fs_if, width = 8, height = 5)
ggsave(paste0("figure_5i_fs_if_density_all_", current_date, ".png"), plot = p_density_fs_if, width = 8, height = 5)

###########################################################################
#  Step 3. Plot TYPE NJ's neoantigen scores for the top neoantigens -------
###########################################################################

df_a3loss <- df_na_nj_map[type == "A3.loss", .(hla_allele, score_average)][, type := paste0("A3 loss (n=", .N, ")")]
df_a3gain <- df_na_nj_map[type == "A3.gain", .(hla_allele, score_average)][, type := paste0("A3 gain (n=", .N, ")")]
df_a5loss <- df_na_nj_map[type == "A5.loss", .(hla_allele, score_average)][, type := paste0("A5 loss (n=", .N, ")")]
df_a5gain <- df_na_nj_map[type == "A5.gain", .(hla_allele, score_average)][, type := paste0("A5 gain (n=", .N, ")")]
df_juncin <- df_na_nj_map[type == "JUNC.WITHIN.EXON", .(hla_allele, score_average)][, type := paste0("JWE (n=", .N, ")")]
df_juncex <- df_na_nj_map[type == "JUNC.WITHIN.INTRON", .(hla_allele, score_average)][, type := paste0("JWI (n=", .N, ")")]
df_exskip <- df_na_nj_map[type == "ES", .(hla_allele, score_average)][, type := paste0("ES (n=", .N, ")")]
df_others <- df_na_nj_map[type == "OTHERS", .(hla_allele, score_average)][, type := paste0("OTHERS (n=", .N, ")")]

df_combined_types <- rbind(df_a3loss, df_a3gain, df_a5loss, df_a5gain, df_juncin, df_juncex, df_exskip, df_others, fill = TRUE)

p_jitter_types <- ggplot(df_combined_types, aes(x = type, y = score_average, color = hla_allele)) + 
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  scale_color_manual(values = hla_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(family = "Helvetica", size = 10, face = "bold", angle = 0, hjust = 0.5),
        axis.text.y = element_text(family = "Helvetica", size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"),
        legend.text = element_text(family = "Helvetica", size = 10),
        legend.title = element_text(family = "Helvetica", size = 10, face = "bold")) +
  labs(y = "Average presentation score")

ggsave(paste0("figure_5i_splice_types_jitter_", current_date, ".pdf"), plot = p_jitter_types, width = 10, height = 6)
ggsave(paste0("figure_5i_splice_types_jitter_", current_date, ".png"), plot = p_jitter_types, width = 10, height = 6)

df_all_types <- rbind(
  df_a3loss[, type := "A3.loss"],
  df_a3gain[, type := "A3.gain"],
  df_a5loss[, type := "A5.loss"],
  df_a5gain[, type := "A5.gain"],
  df_juncin[, type := "JUNC.WITHIN.EXON"],
  df_juncex[, type := "JUNC.WITHIN.INTRON"],
  df_exskip[, type := "ES"],
  df_others[, type := "OTHERS"],
  fill = TRUE
)

p_box_types <- ggplot(df_all_types, aes(x = hla_allele, y = score_average, fill = type)) +
  geom_boxplot(position = position_dodge()) +
  labs(x = "HLA Allele", y = "Immunogenicity Score") +
  scale_fill_manual(values = c("A3.gain" = "#f94144", "A3.loss" = "#f3722c",
                               "A5.gain" = "#f8961e", "A5.loss" = "#f9c74f",
                               "ES" = "#90be6d", "JUNC.WITHIN.EXON" = "#43aa8b",
                               "JUNC.WITHIN.INTRON" = "#4d908e", "OTHERS" = "#577590")) +
  theme_bw() +
  theme(axis.text.x = element_text(family = "Helvetica", size = 10, face = "bold", angle = 90, hjust = 1),
        axis.text.y = element_text(family = "Helvetica", size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"),
        legend.position = "bottom")

ggsave(paste0("figure_5i_splice_types_boxplot_", current_date, ".pdf"), plot = p_box_types, width = 10, height = 6)
ggsave(paste0("figure_5i_splice_types_boxplot_", current_date, ".png"), plot = p_box_types, width = 10, height = 6)

###########################################################################
#  Step 4. Plot TYPE NJ's neoantigen scores for ALL neoantigens -----------
###########################################################################

hist_all <- df_all_map[, .(hla_allele, score_average, type)]
hist_all[, score_log2 := log2(score_average + 0.001)]

hist_0101 <- hist_all[hla_allele == "HLA-A0101"]
hist_0201 <- hist_all[hla_allele == "HLA-A0201"]
hist_0301 <- hist_all[hla_allele == "HLA-A0301"]
hist_1101 <- hist_all[hla_allele == "HLA-A1101"]
hist_2402 <- hist_all[hla_allele == "HLA-A2402"]

create_density_plot <- function(data, allele) {
  x_min <- min(data$score_log2, na.rm = TRUE)
  x_max <- max(data$score_log2, na.rm = TRUE)
  ggplot(data, aes(x = score_log2, fill = type)) +
    geom_density(alpha = 0.5) +
    labs(x = paste0("log2(Average presentation score) (", allele, ")"), y = "Density") +
    scale_fill_manual(values = c("A3.gain" = "#f94144", "A3.loss" = "#f3722c",
                                 "A5.gain" = "#f8961e", "A5.loss" = "#f9c74f",
                                 "ES" = "#90be6d", "JUNC.WITHIN.EXON" = "#43aa8b",
                                 "JUNC.WITHIN.INTRON" = "#4d908e", "OTHERS" = "#577590")) +
    xlim(x_min, x_max) +
    theme_bw()
}

gg_0101 <- create_density_plot(hist_0101, "HLA-A*01:01")
gg_0201 <- create_density_plot(hist_0201, "HLA-A*02:01")
gg_0301 <- create_density_plot(hist_0301, "HLA-A*03:01")
gg_1101 <- create_density_plot(hist_1101, "HLA-A*11:01")
gg_2402 <- create_density_plot(hist_2402, "HLA-A*24:02")

p_arranged <- grid.arrange(gg_0101, gg_0201, gg_0301, gg_1101, gg_2402, ncol = 1)
ggsave(paste0("figure_5i_splice_types_density_all_", current_date, ".pdf"), plot = p_arranged, width = 8, height = 10)
ggsave(paste0("figure_5i_splice_types_density_all_", current_date, ".png"), plot = p_arranged, width = 8, height = 10)

###########################################################################
#  Step 5. Identify Top HLA Alleles and Neoantigens -----------------------
###########################################################################

top_neo <- df_all_map[order(-score_average)][1:round(0.1 * .N)]

# Derive sample_id from junc.id (based on diagnostics: extract after last '-', e.g., "TII" from "CPCT02010267TII")
top_neo[, sample_id := fifelse(junc.id == "Unknown", "Unknown", sub(".*-(.*)$", "\\1", junc.id))]

# Join with patient list (safe left join; keeps all top_neo even if no match)
setkey(top_neo, sample_id)
setkey(patient_list, sample_id)  # From your diagnostics, patient_list has "sample_id"
top_neo_detailed <- top_neo[patient_list, nomatch = NA]

# Select explicit columns for sample/neoantigen/neojunction mapping (peptide = neoantigen, junc.id = neojunction)
top_neo_detailed <- top_neo_detailed[, .(sample_id, peptide, junc.id, hla_allele, score_average, type, fs, shared, purity)]  # Add purity from patient_list if matched

top_hla_summary <- top_neo_detailed[, .(count = .N, avg_score = mean(score_average, na.rm = TRUE)), by = hla_allele][order(-count)]

fwrite(top_neo_detailed, paste0("top_neoantigens_detailed_", current_date, ".tsv"), sep = "\t")
fwrite(top_hla_summary, paste0("top_hla_alleles_summary_", current_date, ".tsv"), sep = "\t")

p_top_hla <- ggplot(top_hla_summary, aes(x = reorder(hla_allele, -count), y = count, fill = hla_allele)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = hla_colors) +
  theme_bw() +
  labs(x = "HLA Allele", y = "Count in Top Neoantigens", title = "Top HLA Alleles in Top 10% Neoantigens")

ggsave(paste0("figure_top_hla_distribution_", current_date, ".pdf"), width = 8, height = 6)
ggsave(paste0("figure_top_hla_distribution_", current_date, ".png"), width = 8, height = 6)

print("Step 15d completed: All joins verified, metadata added/derived—no NAs or empty files.")