#!/usr/bin/env Rscript
# Step 14d: Make histogram 
# September 24, 2025 | Gaurav Raichand | The Institute of Cancer Research

# Purpose: Generate a histogram illustrating the MHCFlurry 2.0 presentation scores

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(readxl)
library(tidyverse)
library(ggplot2)
library(extrafont)
library(data.table)  # For faster reading if needed

# Font Import (run once, then comment out)
# font_import()
loadfonts(device = "pdf")

#  Load Directories -------------------------------------------------------
directory_14 <- Sys.getenv("STEP14_OUTPUT_DIR")
directory_figures <- Sys.getenv("OUTPUT_DIR")

# Load Files --------------------------------------------------------------
setwd(directory_14)
current_date <- format(Sys.Date(), "%Y%m%d")
nmer_08 <- fread(paste0("mhcflurry_08mer_selected_alleles_", current_date, ".tsv"), na = c("", "NA"))
nmer_09 <- fread(paste0("mhcflurry_09mer_selected_alleles_", current_date, ".tsv"), na = c("", "NA"))
nmer_10 <- fread(paste0("mhcflurry_10mer_selected_alleles_", current_date, ".tsv"), na = c("", "NA"))
nmer_11 <- fread(paste0("mhcflurry_11mer_selected_alleles_", current_date, ".tsv"), na = c("", "NA"))

###########################################################################
#  Step 1: Remove duplicate n-mers ----------------------------------------
###########################################################################

nmer_08_unique <- distinct(nmer_08, allele, peptide, n_flank, c_flank, .keep_all = TRUE)
nmer_09_unique <- distinct(nmer_09, allele, peptide, n_flank, c_flank, .keep_all = TRUE)
nmer_10_unique <- distinct(nmer_10, allele, peptide, n_flank, c_flank, .keep_all = TRUE)
nmer_11_unique <- distinct(nmer_11, allele, peptide, n_flank, c_flank, .keep_all = TRUE)

nmer_all_unique <- bind_rows(nmer_08_unique, nmer_09_unique, nmer_10_unique, nmer_11_unique)

###########################################################################
#  Step 2: Make stacked histogram with all n-mers and associated scores ---
###########################################################################

setwd(directory_figures)

for (i in 1:5) {
  if (i == 1) {plot_i <- nmer_08_unique; n_count <- nrow(plot_i); title_i <- paste0("8-mers (n=", n_count, ")"); filename_i <- paste0("histogram_mhcflurry_all_08mer_n", n_count, "_", current_date, ".pdf")}
  if (i == 2) {plot_i <- nmer_09_unique; n_count <- nrow(plot_i); title_i <- paste0("9-mers (n=", n_count, ")"); filename_i <- paste0("histogram_mhcflurry_all_09mer_n", n_count, "_", current_date, ".pdf")}
  if (i == 3) {plot_i <- nmer_10_unique; n_count <- nrow(plot_i); title_i <- paste0("10-mers (n=", n_count, ")"); filename_i <- paste0("histogram_mhcflurry_all_10mer_n", n_count, "_", current_date, ".pdf")}
  if (i == 4) {plot_i <- nmer_11_unique; n_count <- nrow(plot_i); title_i <- paste0("11-mers (n=", n_count, ")"); filename_i <- paste0("histogram_mhcflurry_all_11mer_n", n_count, "_", current_date, ".pdf")}
  if (i == 5) {plot_i <- nmer_all_unique; n_count <- nrow(plot_i); title_i <- paste0("All n-mers (n=", n_count, ")"); filename_i <- paste0("histogram_mhcflurry_all_all_nmers_n", n_count, "_", current_date, ".pdf")}
  
  ggplot(plot_i, aes(x = mhcflurry_presentation_score, fill = allele)) +
    geom_histogram(bins = 50) + 
    theme_minimal() + 
    scale_fill_manual(values=c("#2c9061", "#147ab0", "#cf0e41", "#3d387e", "#dc7320")) + 
    scale_color_manual(values=c("#2c9061", "#147ab0", "#cf0e41", "#3d387e", "#dc7320")) + 
    xlab("MHCFlurry 2.0 Presentation Score") +
    ylab("Count") +
    ggtitle(title_i) +
    theme(plot.title = element_text(size = 20,  family="Helvetica", face = "bold", hjust = 0.5),  # Change the text of the title and adjust to center with hjust
          text = element_text(size = 20,  family="Helvetica"),                                    # Change the text of the legend and the axis  
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +               # Remove the background grids
    labs(fill = "HLA-A allele")
  
  ggsave(filename_i, limitsize = FALSE, width = 8, height = 5)
}

###########################################################################
#  Step 3: Top 10 percentile ----------------------------------------------
###########################################################################
# 1. Generate a list of the top 10 percentile for each n_mer
nmer_08_ordered <- nmer_08_unique[order(nmer_08_unique$mhcflurry_presentation_score, decreasing = TRUE), ]
nmer_09_ordered <- nmer_09_unique[order(nmer_09_unique$mhcflurry_presentation_score, decreasing = TRUE), ]
nmer_10_ordered <- nmer_10_unique[order(nmer_10_unique$mhcflurry_presentation_score, decreasing = TRUE), ]
nmer_11_ordered <- nmer_11_unique[order(nmer_11_unique$mhcflurry_presentation_score, decreasing = TRUE), ]
nmer_all_ordered <- nmer_all_unique[order(nmer_all_unique$mhcflurry_presentation_score, decreasing = TRUE), ]

for (i in 1:5) {
  if (i == 1) {nmer_i <- nmer_08_ordered}
  if (i == 2) {nmer_i <- nmer_09_ordered}
  if (i == 3) {nmer_i <- nmer_10_ordered}
  if (i == 4) {nmer_i <- nmer_11_ordered}
  if (i == 5) {nmer_i <- nmer_all_ordered}
  
  percentile_10 <- as.integer(nrow(nmer_i) / 10)
  
  if (i == 1) {nmer_10percentile_08 <- nmer_i[1:percentile_10, ]}
  if (i == 2) {nmer_10percentile_09 <- nmer_i[1:percentile_10, ]}
  if (i == 3) {nmer_10percentile_10 <- nmer_i[1:percentile_10, ]}
  if (i == 4) {nmer_10percentile_11 <- nmer_i[1:percentile_10, ]}
  if (i == 5) {nmer_10percentile_all <- nmer_i[1:percentile_10, ]}
}

# 2. Generate a pie chart for the number of n_mers representing each allele in the top 10 percentile
pie_08 <- as.data.frame(t(table(nmer_10percentile_08$allele)))[, c(2, 3)]
pie_09 <- as.data.frame(t(table(nmer_10percentile_09$allele)))[, c(2, 3)]
pie_10 <- as.data.frame(t(table(nmer_10percentile_10$allele)))[, c(2, 3)]
pie_11 <- as.data.frame(t(table(nmer_10percentile_11$allele)))[, c(2, 3)]

# Generate a custom table for the all n-mers
nmer_10percentile_all_edit <- NULL
for (i in 1:nrow(nmer_10percentile_all)) {
  row_i <- nmer_10percentile_all %>% slice(i)
  peptide_i <- row_i %>% pull(peptide)
  LEN <- nchar(peptide_i)
  row_i <- row_i %>% mutate(len = LEN)
  nmer_10percentile_all_edit <- rbind(nmer_10percentile_all_edit, row_i)
}
pie_all_nmer <- as.data.frame(t(table(nmer_10percentile_all_edit$len)))[, c(2, 3)]
pie_all_hla <- as.data.frame(t(table(nmer_10percentile_all_edit$allele)))[, c(2, 3)]

for (i in 1:6) {
  if (i == 1) {pie_i <- pie_08; title_i <- "8-mers (Top 10 Percentile)"; filename_i <- paste0("pie_chart_mhcflurry_08mer_top10percentile_", current_date, ".pdf")}
  if (i == 2) {pie_i <- pie_09; title_i <- "9-mers (Top 10 Percentile)"; filename_i <- paste0("pie_chart_mhcflurry_09mer_top10percentile_", current_date, ".pdf")}
  if (i == 3) {pie_i <- pie_10; title_i <- "10-mers (Top 10 Percentile)"; filename_i <- paste0("pie_chart_mhcflurry_10mer_top10percentile_", current_date, ".pdf")}
  if (i == 4) {pie_i <- pie_11; title_i <- "11-mers (Top 10 Percentile)"; filename_i <- paste0("pie_chart_mhcflurry_11mer_top10percentile_", current_date, ".pdf")}
  if (i == 5) {pie_i <- pie_all_hla; title_i <- "All n-mers (Top 10 Percentile)"; filename_i <- paste0("pie_chart_mhcflurry_all_nmers_top10percentile_hla_distribution_", current_date, ".pdf")}
  if (i == 6) {pie_i <- pie_all_nmer; title_i <- "All n-mers (Top 10 Percentile)"; filename_i <- paste0("pie_chart_mhcflurry_all_nmers_top10percentile_len_distribution_", current_date, ".pdf")}
  
  if (i < 6) {
    ggplot(pie_i, aes(x = "", y = Freq, fill = Var2)) +
      scale_fill_manual(values=c("#2c9061", "#147ab0", "#cf0e41", "#3d387e", "#dc7320")) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0) +
      theme_void() + # remove background, grid, numeric labels
      ggtitle(title_i) +
      theme(plot.title = element_text(size = 20,  family="Helvetica", face = "bold", hjust = 0.5),
            text = element_text(size = 20,  family="Helvetica")) + 
      labs(fill = "HLA-A allele")
    
    ggsave(filename_i, limitsize = FALSE, width = 8, height = 5)
  }
  
  if (i == 6) {
    ggplot(pie_i, aes(x = "", y = Freq, fill = Var2)) +
      scale_fill_manual(values=c("#e6c229", "#f17105", "#d11149", "#6610f2")) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0) +
      theme_void() + # remove background, grid, numeric labels
      ggtitle(title_i) +
      theme(plot.title = element_text(size = 20,  family="Helvetica", face = "bold", hjust = 0.5),
            text = element_text(size = 20,  family="Helvetica")) + 
      labs(fill = "n-mer Length")
    
    ggsave(filename_i, limitsize = FALSE, width = 8, height = 5)
  }
}

nmer_10percentile_all_edit$len <- as.character(nmer_10percentile_all_edit$len)
nmer_10percentile_all_edit_len <- nmer_10percentile_all_edit
nmer_10percentile_all_edit_len$len[nmer_10percentile_all_edit_len$len == "8"] <- "08"
nmer_10percentile_all_edit_len$len[nmer_10percentile_all_edit_len$len == "9"] <- "09"
nmer_10percentile_all_edit_len <- nmer_10percentile_all_edit_len[order(nmer_10percentile_all_edit_len$len), ]

# 3. Generate a histogram chart for the number of n_mers representing each allele in the top 10 percentile
for (i in 1:6) {
  if (i == 1) {plot_i <- nmer_10percentile_08; title_i <- "8-mers (Top 10 Percentile)"; filename_i <- paste0("histogram_mhcflurry_top10percentile_08mer_", current_date, ".pdf")}
  if (i == 2) {plot_i <- nmer_10percentile_09; title_i <- "9-mers (Top 10 Percentile)"; filename_i <- paste0("histogram_mhcflurry_top10percentile_09mer_", current_date, ".pdf")}
  if (i == 3) {plot_i <- nmer_10percentile_10; title_i <- "10-mers (Top 10 Percentile)"; filename_i <- paste0("histogram_mhcflurry_top10percentile_10mer_", current_date, ".pdf")}
  if (i == 4) {plot_i <- nmer_10percentile_11; title_i <- "11-mers (Top 10 Percentile)"; filename_i <- paste0("histogram_mhcflurry_top10percentile_11mer_", current_date, ".pdf")}
  if (i == 5) {plot_i <- nmer_10percentile_all_edit; title_i <- "All n-mers (Top 10 Percentile)"; filename_i <- paste0("histogram_mhcflurry_top10percentile_all_nmers_hla_distribution_", current_date, ".pdf")}
  if (i == 6) {plot_i <- nmer_10percentile_all_edit_len; title_i <- "All n-mers (Top 10 Percentile)"; filename_i <- paste0("histogram_mhcflurry_top10percentile_all_nmers_len_distribution_", current_date, ".pdf")}
  
  if (i < 6) {
    ggplot(plot_i, aes(x = mhcflurry_presentation_score, fill = allele)) +
      geom_histogram(bins = 30) + 
      theme_minimal() + 
      scale_fill_manual(values=c("#2c9061", "#147ab0", "#cf0e41", "#3d387e", "#dc7320")) + 
      scale_color_manual(values=c("#2c9061", "#147ab0", "#cf0e41", "#3d387e", "#dc7320")) + 
      xlab("MHCFlurry 2.0 Presentation Score") +
      ylab("Count") +
      ggtitle(title_i) +
      theme(plot.title = element_text(size = 20,  family="Helvetica", face = "bold", hjust = 0.5),
            text = element_text(size = 20,  family="Helvetica"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      labs(fill = "HLA-A allele")
    
    ggsave(filename_i, limitsize = FALSE, width = 8, height = 5)
  }
  
  if (i == 6) {
    ggplot(plot_i, aes(x = mhcflurry_presentation_score, fill = len)) +
      geom_histogram(bins = 30) + 
      theme_minimal() + 
      scale_fill_manual(values=c("#e6c229", "#f17105", "#d11149", "#6610f2")) + 
      scale_color_manual(values=c("#e6c229", "#f17105", "#d11149", "#6610f2")) + 
      xlab("MHCFlurry 2.0 Presentation Score") +
      ylab("Count") +
      ggtitle(title_i) +
      theme(plot.title = element_text(size = 20,  family="Helvetica", face = "bold", hjust = 0.5),
            text = element_text(size = 20,  family="Helvetica"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      labs(fill = "n-mer Length")
    
    ggsave(filename_i, limitsize = FALSE, width = 8, height = 5)
  }
}

# Export files
setwd(directory_14)

print("finished")