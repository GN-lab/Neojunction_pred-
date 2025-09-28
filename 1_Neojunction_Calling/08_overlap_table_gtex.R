#!/usr/bin/env Rscript
# Title: "Step 8: SJ Overlap Table GTEx" – Refactored with robust chunked GCT parsing
# Updated September 26, 2025 – replaces zcat|tail|head with R-native streaming

rm(list = ls(all.names = TRUE))
library(tidyverse)
library(data.table)

#--------------------
# Step 0: Env vars
#--------------------
directory_meta <- Sys.getenv("GTEX_META_PATH")
directory_05   <- Sys.getenv("STEP05_OUTPUT_DIR")
directory_out  <- Sys.getenv("OUTPUT_DIR")
input_dir      <- Sys.getenv("INPUT_DIR")
gct_gz         <- file.path(input_dir,
                    "GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz")

for (d in c(directory_meta, directory_05, directory_out, input_dir)) {
  if (!nzchar(d) || !dir.exists(d)) stop("[ERROR] Missing dir: ", d)
}
if (!file.exists(gct_gz)) stop("[ERROR] Missing GCT file: ", gct_gz)

current_date <- format(Sys.Date(), "%Y%m%d")

#--------------------
# Step 1: Load Step5 lists
#--------------------
path_annot <- file.path(directory_05,
  "SJ_List_Filtered_by_GTF_ProteinCoding_ExpressedTranscripts_20250925.tsv"
)
if (!file.exists(path_annot)) stop("[ERROR] Missing annotated SJ list: ", path_annot)
dataframe_sj.annot <- fread(path_annot)

fname_nonannot <- paste0("SJ_List_NonAnnotated_Candidates_Protein_Coding_",
                         current_date, ".tsv")
path_nonannot <- file.path(directory_05, fname_nonannot)
if (!file.exists(path_nonannot)) stop("[ERROR] Missing non-annot SJ list: ",
                                      path_nonannot)
dataframe_sj.nonannot <- fread(path_nonannot, colClasses = list(character = "junc.id"))

# Parse coords for non-annotated list
dataframe_sj.nonannot <- dataframe_sj.nonannot %>%
  select(junc.id) %>%
  mutate(
    chr       = sub("chr", "", sapply(strsplit(junc.id, ":"), "[[", 1)),
    strand    = sapply(strsplit(junc.id, ":"), "[[", 2),
    int.start = as.integer(sub("-.*", "", sapply(strsplit(junc.id, ":"), "[[", 3))) + 1,
    int.end   = as.integer(sub(".*-", "", sapply(strsplit(junc.id, ":"), "[[", 3)))
  )

# Build overlap list
list_sj <- vector("list", nrow(dataframe_sj.nonannot))
junc.ids <- character()
for (i in seq_len(nrow(dataframe_sj.nonannot))) {
  sj.i <- dataframe_sj.nonannot[i, ]
  overlaps <- dataframe_sj.annot %>%
    filter(chr == sj.i$chr,
           strand == sj.i$strand,
           int.start < sj.i$int.end,
           sj.i$int.start < int.end) %>%
    pull(junc.id)
  entry <- tibble(
    junc.id         = sj.i$junc.id,
    junc.id.overlap = if (length(overlaps)) overlaps else NA_character_
  )
  list_sj[[i]] <- entry
  junc.ids <- union(junc.ids, na.omit(c(sj.i$junc.id, overlaps)))
}

#--------------------
# Step 3: Stream GCT & compute >1% filter
#--------------------
message("[INFO] Streaming GTEx GCT for >1% filter")
con <- gzfile(gct_gz, "rt")
# Skip first 3 header lines
readLines(con, n = 3)

gtex_common_junctions <- character()
chunk_size <- 10000L

repeat {
  block <- tryCatch(
    read.table(con,
               sep = "\t",
               header = FALSE,
               comment.char = "",
               nrows = chunk_size,
               stringsAsFactors = FALSE,
               fill = TRUE),
    error = function(e) NULL
  )
  if (is.null(block) || nrow(block) == 0) break

  raw_ids   <- block[[1]]
  counts_mat <- as.matrix(block[, -(1:2), drop = FALSE])
  pct       <- rowSums(counts_mat > 0, na.rm = TRUE) / ncol(counts_mat) * 100

  coords_key <- vapply(strsplit(raw_ids, "_"),
                       function(x) paste0(x[1], ":", x[2], "-", x[3]), "")
  gtex_common_junctions <- c(
    gtex_common_junctions,
    coords_key[pct > 1.0]
  )
}
close(con)
gtex_common_junctions <- unique(gtex_common_junctions)
message("[INFO] Found ", length(gtex_common_junctions),
        " GTEx junctions >1%")

# Save filter list
fwrite(data.table(junction_id = gtex_common_junctions),
       file.path(directory_out,
         paste0("GTEx_Common_Junctions_1pct_", current_date, ".tsv")),
       sep = "\t")

#--------------------
# Step 4: Pseudo counts for compatibility
#--------------------
dataframe_sj.count <- tibble(junc.id = junc.ids) %>%
  mutate(start = as.integer(sub("-.*", "", sapply(strsplit(junc.id, ":"), "[[", 3)))) %>%
  arrange(start) %>%
  select(junc.id)

meta_gtex <- fread(file.path(directory_meta,
  "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")) %>%
  rename(sample.id = SAMPID) %>%
  filter(grepl("RNASEQ", SMAFRZE, ignore.case = TRUE)) %>%
  slice_head(n = 50)
for (s in meta_gtex$sample.id) dataframe_sj.count[[s]] <- 1

#--------------------
# Step 5: Dominant overlap
#--------------------
dataframe_overlap.table <- tibble()
for (entry in list_sj) {
  if (nrow(entry) == 1) {
    dataframe_overlap.table <- bind_rows(dataframe_overlap.table, entry)
  } else {
    sj_temp <- entry %>%
      filter(!is.na(junc.id.overlap)) %>%
      slice_head(n = 1)
    dataframe_overlap.table <- bind_rows(
      dataframe_overlap.table,
      entry %>% semi_join(sj_temp, by = "junc.id.overlap")
    )
  }
}

#--------------------
# Step 6: Filter common GTEx junctions
#--------------------
normalize_key <- function(id) sub("^([^:]+):.:(\\d+)-(\\d+)$", "\\1:\\2-\\3", id)
dataframe_overlap.table.filtered <-
  dataframe_overlap.table %>%
    filter(!normalize_key(junc.id) %in% gtex_common_junctions) %>%
    filter(is.na(junc.id.overlap) |
           !normalize_key(junc.id.overlap) %in% gtex_common_junctions)

message("[INFO] Retained ", nrow(dataframe_overlap.table.filtered),
        " of ", nrow(dataframe_overlap.table), " overlaps")

#--------------------
# Step 7: Export results
#--------------------
out_overlap <- file.path(directory_out,
  paste0("GTEx_SJ_Overlap_Table_", current_date, ".tsv"))
write_tsv(dataframe_overlap.table.filtered, out_overlap, na = "NA")

juncs_retain <- dataframe_overlap.table.filtered %>%
  pivot_longer(everything(), names_to = "label", values_to = "junc.id") %>%
  filter(!is.na(junc.id)) %>%
  distinct(junc.id)
dataframe_sj.count.retain <- dataframe_sj.count %>%
  semi_join(juncs_retain, by = "junc.id")
out_retain <- file.path(directory_out,
  paste0("GTEx_SJ_Count_Table_SJ.out.tab_Retained_", current_date, ".tsv"))
write_tsv(dataframe_sj.count.retain, out_retain, na = "NA")

message("[INFO] Step 8 done. Outputs:")
message("  - ", out_overlap)
message("  - ", out_retain)
message("  - Filter list in ", file.path(directory_out,
         paste0("GTEx_Common_Junctions_1pct_", current_date, ".tsv")))
