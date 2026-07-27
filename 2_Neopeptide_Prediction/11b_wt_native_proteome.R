#!/usr/bin/env Rscript
# Title: "Step 11b: Patient-Specific WT Native Proteome Peptide Generation"
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# Purpose: Generate 8-11mer WT peptides from the UNION of each patient's own
#          expressed, protein-coding transcriptome (TPM >= MIN_TPM in that
#          patient), rather than the entire canonical UniProt proteome.
#
#          Rationale: the previous version of this script tiled the whole
#          UniProt canonical FASTA -- not patient-specific, no TPM filter,
#          and it recreates peptides for genes no patient in this cohort
#          ever expresses. This version instead:
#            1. Reads the per-sample TPM matrix from Step 03
#               (TPM_Filter10_ProteinCodingTx_*.tsv -- one column per patient)
#            2. For each transcript, keeps it if ANY patient's OWN TPM value
#               is >= MIN_TPM (i.e. "expressed in at least one patient")
#            3. Pulls the real protein sequence for each kept transcript from
#               EnsDb (already translated -- no genomic reconstruction needed,
#               since a transcript's protein sequence doesn't vary by patient,
#               only whether a given patient expresses it does)
#            4. Tiles into 8-11mers, deduplicates
#
#          This IS the "union of patient-specific WT" you wanted: each
#          transcript enters the pool only because some patient's own
#          expression justified it, and the final n-mer set is the union
#          across the whole 307-sample cohort.
#
#          Output files are written with EXACTLY the same filenames that
#          Steps 13a and 14a expect -- so nothing downstream needs to change:
#
#            2023_0812_hlathenalist_msic_08mers_wt.tsv   (-> 13a NetMHCpan)
#            2023_0812_hlathenalist_msic_09mers_wt.tsv
#            2023_0812_hlathenalist_msic_10mers_wt.tsv
#            2023_0812_hlathenalist_msic_11mers_wt.tsv
#
#          This script runs AFTER Step 03 and BEFORE Step 13a.
#
# Config variables used (all set in config.sh):
#   OUTPUT_DIR              -- where Step 03's TPM file lives, and where
#                               13a/14a expect the _wt.tsv files
#   MIN_TPM                 -- expression threshold (default 10, same as
#                               used in Step 03's cohort-median filter)
#   SLURM_CPUS_PER_TASK     -- parallelism (set by SLURM automatically)

###########################################################################
# 0. Packages
###########################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(dplyr)
  library(tibble)
  library(doParallel)
  library(foreach)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)
})

start_time <- proc.time()

###########################################################################
# 1. Paths / config
###########################################################################

output_dir <- Sys.getenv("OUTPUT_DIR")
if (output_dir == "") stop("[ERROR] OUTPUT_DIR is not set. Source config.sh first.")

min_tpm <- as.numeric(Sys.getenv("MIN_TPM", unset = "10"))

num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
num_cores <- max(1L, num_cores - 1L)   # leave one core free for OS

cat("[CONFIG] OUTPUT_DIR:   ", output_dir, "\n")
cat("[CONFIG] MIN_TPM:      ", min_tpm, "\n")
cat("[CONFIG] Cores to use: ", num_cores, "\n\n")

###########################################################################
# 2. Load per-sample TPM matrix from Step 03 and derive the union of
#    patient-expressed transcripts (ANY patient's own TPM >= min_tpm)
###########################################################################

setwd(output_dir)
tpm_pattern <- "^TPM_Filter[0-9]+_ProteinCodingTx_[0-9]{8}\\.tsv$"
tpm_files   <- list.files(pattern = tpm_pattern)
if (length(tpm_files) == 0) {
  stop("[ERROR] No file matching '", tpm_pattern, "' found in ", output_dir,
       " -- run Step 03 first.")
}
tpm_file <- sort(tpm_files, decreasing = TRUE)[1]
cat("[INFO] Loading per-sample TPM matrix:", tpm_file, "\n")

tpm_dt <- fread(tpm_file)
if (!"enst" %in% colnames(tpm_dt)) stop("[ERROR] Expected an 'enst' column in ", tpm_file)

sample_cols <- setdiff(colnames(tpm_dt), "enst")
cat("[INFO] Transcripts in file:", nrow(tpm_dt), "| Patients:", length(sample_cols), "\n")

# Per-transcript max TPM across all patients -- a transcript is kept if it
# clears the threshold in AT LEAST ONE patient's own value. This is what
# turns "per-patient expression" into a cohort union: each transcript's
# inclusion is justified by some specific patient's real expression, not
# a cohort-wide average.
tpm_dt[, max_tpm := do.call(pmax, c(.SD, list(na.rm = TRUE))), .SDcols = sample_cols]
tpm_dt[, n_patients_expressing := rowSums(.SD >= min_tpm, na.rm = TRUE), .SDcols = sample_cols]

expressed <- tpm_dt[max_tpm >= min_tpm, .(enst, max_tpm, n_patients_expressing)]
cat("[INFO] Transcripts expressed (TPM >=", min_tpm, ") in >=1 patient:",
    nrow(expressed), "out of", nrow(tpm_dt), "\n\n")

if (nrow(expressed) == 0) {
  stop("[ERROR] No transcripts passed the per-patient TPM >= ", min_tpm, " filter.")
}

enst_list <- unique(expressed$enst)

###########################################################################
# 3. Pull protein sequences for these transcripts from EnsDb
#    (already translated -- no genomic reconstruction needed here, since
#    Step 11 already handles that for the ALT/junction-specific case; this
#    is the normal/canonical protein for each expressed transcript)
###########################################################################

cat("[INFO] Querying EnsDb for", length(enst_list), "transcript protein sequences...\n")

# Chunk the tx_id filter to keep each ensembldb query reasonably sized
chunk_size <- 2000
chunks     <- split(enst_list, ceiling(seq_along(enst_list) / chunk_size))

prot_list <- vector("list", length(chunks))
for (i in seq_along(chunks)) {
  prot_list[[i]] <- tryCatch({
    proteins(
      EnsDb.Hsapiens.v86,
      columns = c("tx_id", "protein_sequence"),
      filter  = TxIdFilter(chunks[[i]])
    ) %>% as_tibble()
  }, error = function(e) {
    cat("[WARN] Chunk", i, "failed:", conditionMessage(e), "\n")
    tibble(tx_id = character(0), protein_sequence = character(0))
  })
  if (i %% 10 == 0) cat("[INFO]   ", i, "/", length(chunks), "chunks done\n")
}

prot_df <- bind_rows(prot_list) %>%
  distinct(tx_id, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(protein_sequence), nchar(protein_sequence) > 0)

n_mapped <- nrow(prot_df)
n_missed <- length(enst_list) - n_mapped
cat("[INFO] Mapped", n_mapped, "/", length(enst_list),
    "expressed transcripts to a protein sequence (", n_missed,
    "had no coding sequence in EnsDb -- likely non-coding biotype or a",
    "transcript version EnsDb doesn't carry, and are dropped).\n\n")

if (n_mapped == 0) stop("[ERROR] Zero transcripts mapped to protein sequences -- check tx_id format.")

fasta_seqs <- prot_df$protein_sequence

# Quality filter: remove sequences with ambiguous/stop characters mid-sequence
before <- length(fasta_seqs)
fasta_seqs <- fasta_seqs[!grepl("[*BJOUXZ]", fasta_seqs)]
cat("[INFO] After quality filter:", length(fasta_seqs),
    "(removed", before - length(fasta_seqs), ")\n\n")

###########################################################################
# 4. Tile all sequences into 8-11mers in parallel
#    (identical tiling logic to the previous version of this script --
#    only the SOURCE of fasta_seqs has changed, from UniProt canonical to
#    patient-expressed EnsDb transcripts)
###########################################################################

cat("[INFO] Starting peptide tiling (", num_cores, "cores)...\n")

cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(cl, varlist = "fasta_seqs", envir = environment())

all_mers_list <- foreach(
  h         = 8:11,
  .combine  = "c",
  .packages = c("data.table", "stringr")
) %dopar% {

  tile_one <- function(seq, h) {
    len <- nchar(seq)
    if (is.na(len) || len < h) return(NULL)
    starts <- seq_len(len - h + 1L)
    data.table(
      n_mer   = substring(seq, starts, starts + h - 1L),
      n_flank = substring(seq, pmax(1L, starts - 30L), starts - 1L),
      c_flank = substring(seq, starts + h, pmin(len, starts + h + 29L))
    )
  }

  chunks_t <- lapply(fasta_seqs, tile_one, h = h)
  dt       <- rbindlist(chunks_t, use.names = FALSE, fill = FALSE)

  dt <- dt[!duplicated(dt$n_mer)]
  dt <- dt[!grepl("[^ACDEFGHIKLMNPQRSTVWY]", n_mer)]

  list(dt)
}

stopCluster(cl)

cat("[INFO] Tiling complete.\n\n")

###########################################################################
# 5. Write output files in exactly the format Steps 13a and 14a expect
###########################################################################

mer_lengths <- 8:11

for (idx in seq_along(mer_lengths)) {

  h      <- mer_lengths[idx]
  lenpad <- sprintf("%02d", h)
  dt_h   <- all_mers_list[[idx]]

  n_raw <- nrow(dt_h)
  cat(sprintf("[INFO] %smer: %d unique peptides before final filters\n", lenpad, n_raw))

  dt_h[, ctex_up := str_pad(n_flank, width = 30L, side = "left",  pad = "-")]
  dt_h[, ctex_dn := str_pad(c_flank, width = 30L, side = "right", pad = "-")]
  dt_h[, TPM     := ""]

  out_dt   <- dt_h[, .(n_mer, ctex_up, ctex_dn, TPM)]
  out_file <- file.path(output_dir,
                        paste0("2023_0812_hlathenalist_msic_", lenpad, "mers_wt.tsv"))

  fwrite(out_dt, out_file, sep = "\t", na = "NA",
         col.names = TRUE, quote = FALSE)

  cat(sprintf("[INFO] Written: %s (%d rows)\n", out_file, nrow(out_dt)))
}

###########################################################################
# 6. Done
###########################################################################

runtime <- proc.time() - start_time
cat("\n[DONE] Step 11b (patient-specific) complete.\n")
cat(sprintf("[DONE] Total runtime: %.1f seconds (%.1f minutes)\n",
            runtime[3], runtime[3] / 60))
cat("[DONE] Source: ", n_mapped, " transcripts expressed (TPM >= ", min_tpm,
    ") in at least one of ", length(sample_cols), " patients.\n", sep = "")
cat("[DONE] The four _wt.tsv files in", output_dir,
    "now reflect the patient-specific expressed transcriptome union.\n")
cat("[DONE] Steps 13a and 14a will pick them up automatically -- no changes needed.\n")
