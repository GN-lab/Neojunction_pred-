#!/usr/bin/env Rscript
# Step 15a: Cross-analysis -- NetMHCpan 4.2 x MHCflurry 2.0 Concordance
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# TWO-STEP LOGIC:
#
# STEP 1 -- Clean ALT immunopeptidome:
#   Load ALL ALT predictions (NetMHCpan + MHCflurry, full universe).
#   Remove ANY peptide sequence found in WT predictions from EITHER tool
#   at ANY binding level (NB, WB, SB -- doesn't matter).
#   What remains is the clean tumour-specific peptide universe.
#
# STEP 2 -- Concordance tiering (applied separately to ALT and WT):
#   Tier 1 (High confidence):   NMP EL rank < 0.5% AND MHCflurry < 500nM
#   Tier 2 (Medium confidence): NMP EL rank < 2.0% AND MHCflurry < 500nM
#   Tier 3 (Discordant):        One tool binder, other NB -- flagged
#   Excluded:                   Both NB -- dropped
#
# Output files:
#   ALT (tumour-specific, WT-filtered):
#     alt_concordance_tier1_YYYYMMDD.tsv   -- both tools strong
#     alt_concordance_tier2_YYYYMMDD.tsv   -- moderate agreement
#     alt_concordance_tier3_YYYYMMDD.tsv   -- discordant
#     alt_concordance_all_YYYYMMDD.tsv     -- all tiers combined
#   WT (native immunopeptidome reference):
#     wt_concordance_tier1_YYYYMMDD.tsv
#     wt_concordance_tier2_YYYYMMDD.tsv
#     wt_concordance_tier3_YYYYMMDD.tsv
#     wt_concordance_all_YYYYMMDD.tsv
#   Shared:
#     cross_alg_all_nmers_YYYYMMDD.tsv    -- clean ALT full join (for 15b/c/d)
#     wt_exclusion_peptides_YYYYMMDD.txt  -- list of excluded WT peptides

###########################################################################
#  Step 0: Packages and config
###########################################################################

rm(list = ls(all.names = TRUE))
library(data.table)

output_dir   <- Sys.getenv("OUTPUT_DIR")
directory_13 <- Sys.getenv("STEP13_OUTPUT_DIR")
directory_14 <- Sys.getenv("STEP14_OUTPUT_DIR")
directory_15 <- Sys.getenv("STEP15_OUTPUT_DIR")

if (nchar(output_dir)   == 0) stop("OUTPUT_DIR not set -- source config.sh first")
if (nchar(directory_13) == 0) directory_13 <- output_dir
if (nchar(directory_14) == 0) directory_14 <- output_dir
if (nchar(directory_15) == 0) directory_15 <- output_dir

dir.create(directory_15, showWarnings = FALSE, recursive = TRUE)

# Binding thresholds
NMP_SB_RANK <- 0.5    # NetMHCpan EL rank strictly < this = Strong Binder
NMP_WB_RANK <- 2.0    # NetMHCpan EL rank strictly < this = Weak Binder
MHC_AFF_NM  <- 500    # MHCflurry affinity strictly < this nM = binder

###########################################################################
#  Step 1: Detect file dates from existing files
###########################################################################

# 13a raw date: netmhcpan_09mer_YYYY_MMDD.tsv
nmp_raw_scan <- list.files(directory_13,
                            pattern = "^netmhcpan_09mer_[0-9]{4}_[0-9]{4}\\.tsv$")
if (length(nmp_raw_scan) == 0)
  stop("[ERROR] No netmhcpan_09mer_YYYY_MMDD.tsv in: ", directory_13)
nmp_run_date <- sub("^netmhcpan_09mer_([0-9]{4}_[0-9]{4})\\.tsv$", "\\1",
                    nmp_raw_scan[which.max(
                      file.info(file.path(directory_13, nmp_raw_scan))$mtime)])

# 14b ALT date: 09mers_flank_mhcflurry_YYYY_MMDD.csv
mhc_scan <- list.files(directory_14,
                        pattern = "^09mers_flank_mhcflurry_[0-9]{4}_[0-9]{4}\\.csv$")
if (length(mhc_scan) == 0)
  stop("[ERROR] No 09mers_flank_mhcflurry_YYYY_MMDD.csv in: ", directory_14)
mhc_date <- sub("^09mers_flank_mhcflurry_([0-9]{4}_[0-9]{4})\\.csv$", "\\1",
                mhc_scan[which.max(
                  file.info(file.path(directory_14, mhc_scan))$mtime)])

# 14b WT date: 09mers_flank_mhcflurry_wt_YYYY_MMDD.csv
mhc_wt_scan <- list.files(directory_14,
                           pattern = "^09mers_flank_mhcflurry_wt_[0-9]{4}_[0-9]{4}\\.csv$")
if (length(mhc_wt_scan) == 0)
  stop("[ERROR] No 09mers_flank_mhcflurry_wt_YYYY_MMDD.csv in: ", directory_14)
mhc_wt_date <- sub("^09mers_flank_mhcflurry_wt_([0-9]{4}_[0-9]{4})\\.csv$", "\\1",
                   mhc_wt_scan[which.max(
                     file.info(file.path(directory_14, mhc_wt_scan))$mtime)])

current_date <- format(Sys.Date(), "%Y%m%d")

cat("=== File dates detected ===\n")
cat("[INFO] NetMHCpan raw (13a):      ", nmp_run_date, "\n")
cat("[INFO] MHCFlurry ALT (14b):      ", mhc_date,     "\n")
cat("[INFO] MHCFlurry WT  (14b):      ", mhc_wt_date,  "\n")
cat("[INFO] Output date:              ", current_date,  "\n")
cat("[INFO] Tier 1: NMP EL rank <",  NMP_SB_RANK, "% AND MHCflurry <", MHC_AFF_NM, "nM\n")
cat("[INFO] Tier 2: NMP EL rank <",  NMP_WB_RANK, "% AND MHCflurry <", MHC_AFF_NM, "nM\n")

###########################################################################
#  Helper: load NetMHCpan raw TSV (all binding levels)
###########################################################################

load_nmp_raw <- function(len, type = "alt") {
  if (type == "wt") {
    f <- file.path(directory_13,
                   paste0("netmhcpan_", len, "mer_wt_", nmp_run_date, ".tsv"))
  } else {
    f <- file.path(directory_13,
                   paste0("netmhcpan_", len, "mer_", nmp_run_date, ".tsv"))
  }
  if (!file.exists(f)) { warning("[WARN] Missing: ", f); return(data.table()) }
  cat(sprintf("[INFO]   Loading NMP %s %smer...\n", type, len))
  dt <- fread(f, na.strings = c("", "NA"),
              select = c("allele", "peptide", "netmhcpan_EL_score",
                         "netmhcpan_EL_rank", "netmhcpan_BA_score",
                         "netmhcpan_BA_rank", "binder"))
  dt[, netmhcpan_EL_rank  := as.numeric(netmhcpan_EL_rank)]
  dt[, netmhcpan_EL_score := as.numeric(netmhcpan_EL_score)]
  dt
}

###########################################################################
#  Helper: load MHCflurry raw flank CSV (full prediction universe)
###########################################################################

load_mhc_raw <- function(len, type = "alt") {
  if (type == "wt") {
    f <- file.path(directory_14,
                   paste0(len, "mers_flank_mhcflurry_wt_", mhc_wt_date, ".csv"))
  } else {
    f <- file.path(directory_14,
                   paste0(len, "mers_flank_mhcflurry_", mhc_date, ".csv"))
  }
  if (!file.exists(f)) { warning("[WARN] Missing: ", f); return(data.table()) }
  cat(sprintf("[INFO]   Loading MHC %s %smer...\n", type, len))
  dt <- fread(f, na.strings = c("", "NA"),
              select = c("allele", "peptide", "n_flank", "c_flank",
                         "mhcflurry_affinity", "mhcflurry_presentation_score",
                         "mhcflurry_presentation_percentile"))
  dt[, mhcflurry_affinity           := as.numeric(mhcflurry_affinity)]
  dt[, mhcflurry_presentation_score := as.numeric(mhcflurry_presentation_score)]
  dt
}

###########################################################################
#  Helper: assign concordance tiers to a joined data.table
###########################################################################

assign_tiers <- function(dt) {
  dt[, nmp_EL_rank := as.numeric(netmhcpan_EL_rank)]
  dt[, mhc_aff     := as.numeric(mhcflurry_affinity)]

  # NMP call -- strict < thresholds
  dt[, nmp_call := fcase(
    nmp_EL_rank < NMP_SB_RANK, "SB",
    nmp_EL_rank < NMP_WB_RANK, "WB",
    default = "NB"
  )]

  # MHCflurry call -- strict < 500nM
  dt[, mhc_call := fcase(
    mhc_aff < MHC_AFF_NM, "binder",
    default = "NB"
  )]

  # Concordance tier
  dt[, concordance_tier := fcase(
    nmp_call == "SB" & mhc_call == "binder",              "Tier1_HighConfidence",
    nmp_call == "WB" & mhc_call == "binder",              "Tier2_MediumConfidence",
    nmp_call == "SB" & mhc_call == "NB",                  "Tier3_Discordant_NMPstrong",
    nmp_call %in% c("WB","NB") & mhc_call == "binder",   "Tier3_Discordant_MHCstrong",
    default = "Excluded"
  )]

  # Combined score: mean of NMP EL score + MHCflurry presentation score (both 0-1)
  dt[, combined_score := rowMeans(
    cbind(as.numeric(netmhcpan_EL_score),
          as.numeric(mhcflurry_presentation_score)),
    na.rm = TRUE
  )]
  dt
}

###########################################################################
#  Helper: join NMP + MHCflurry on allele + peptide, keep best row each
###########################################################################

join_tools <- function(nmp_dt, mhc_dt) {
  # Best row per allele+peptide in each tool
  # setorder() + unique() avoids the very large .SD allocation made by .SD[1].
  setorder(nmp_dt, netmhcpan_EL_rank, na.last = TRUE)
  setorder(mhc_dt, mhcflurry_affinity, na.last = TRUE)
  nmp_best <- unique(nmp_dt, by = c("allele", "peptide"))
  mhc_best <- unique(mhc_dt, by = c("allele", "peptide"))
  merge(nmp_best, mhc_best, by = c("allele", "peptide"), all = TRUE)
}

###########################################################################
#  Helper: split tiered table into files and write
###########################################################################

write_tiers <- function(dt, prefix) {
  setwd(directory_15)

  tier1 <- dt[concordance_tier == "Tier1_HighConfidence"]
  tier2 <- dt[concordance_tier == "Tier2_MediumConfidence"]
  tier3 <- dt[grepl("Tier3", concordance_tier)]
  all_t <- dt[concordance_tier != "Excluded"]

  setorder(tier1, nmp_EL_rank)
  setorder(tier2, nmp_EL_rank)
  setorder(tier3, concordance_tier, nmp_EL_rank)
  setorder(all_t, concordance_tier, nmp_EL_rank)

  fwrite(tier1, paste0(prefix, "_tier1_", current_date, ".tsv"), sep="\t", na="NA")
  fwrite(tier2, paste0(prefix, "_tier2_", current_date, ".tsv"), sep="\t", na="NA")
  fwrite(tier3, paste0(prefix, "_tier3_", current_date, ".tsv"), sep="\t", na="NA")
  fwrite(all_t, paste0(prefix, "_all_",   current_date, ".tsv"), sep="\t", na="NA")

  cat(sprintf("\n  [%s] Tier 1 (High confidence):   %d rows\n", prefix, nrow(tier1)))
  cat(sprintf("  [%s] Tier 2 (Medium confidence): %d rows\n", prefix, nrow(tier2)))
  cat(sprintf("  [%s] Tier 3 (Discordant):        %d rows\n", prefix, nrow(tier3)))
  cat(sprintf("  [%s] Excluded (both NB):         %d rows\n", prefix,
              nrow(dt[concordance_tier == "Excluded"])))

  list(tier1=tier1, tier2=tier2, tier3=tier3, all=all_t)
}

###########################################################################
#  ======================== STEP 1: BUILD WT EXCLUSION SET ===============
#  Load ALL WT predictions from BOTH tools at ALL binding levels.
#  Every peptide sequence seen in WT is added to the exclusion set.
#  This is done BEFORE any ALT processing.
###########################################################################

cat("\n========== STEP 1: Building WT exclusion set (by length) ==========\n")

# Never retain the 345M-row WT prediction tables.  Only the unique peptide
# vectors are kept here; the full data are re-read one length at a time in 2b.
wt_exclusion_parts <- vector("list", 4L)
lengths <- c("08", "09", "10", "11")
for (i in seq_along(lengths)) {
  len <- lengths[i]
  cat(sprintf("[INFO] Collecting unique WT peptides for %smer...\n", len))
  nmp_file <- file.path(directory_13,
                        paste0("netmhcpan_", len, "mer_wt_", nmp_run_date, ".tsv"))
  mhc_file <- file.path(directory_14,
                        paste0(len, "mers_flank_mhcflurry_wt_", mhc_wt_date, ".csv"))
  # Step 1 needs no scores: reading one column instead of fourteen materially
  # reduces both resident memory and I/O.
  nmp_peptides <- fread(nmp_file, select = "peptide", na.strings = c("", "NA"))
  mhc_peptides <- fread(mhc_file, select = "peptide", na.strings = c("", "NA"))
  wt_exclusion_parts[[i]] <- unique(c(nmp_peptides$peptide, mhc_peptides$peptide))
  cat(sprintf("[INFO]   %smer rows: NMP=%d, MHC=%d; unique peptides=%d\n",
              len, nrow(nmp_peptides), nrow(mhc_peptides),
              length(wt_exclusion_parts[[i]])))
  rm(nmp_peptides, mhc_peptides)
  gc(verbose = FALSE)
}
wt_exclusion <- unique(unlist(wt_exclusion_parts, use.names = FALSE))
rm(wt_exclusion_parts)
gc(verbose = FALSE)
cat(sprintf("\n[INFO] Total unique WT peptides to exclude: %d\n", length(wt_exclusion)))

# Write exclusion list for traceability
setwd(directory_15)
writeLines(wt_exclusion,
           paste0("wt_exclusion_peptides_", current_date, ".txt"))
cat(sprintf("[INFO] WT exclusion list written: wt_exclusion_peptides_%s.txt\n",
            current_date))

###########################################################################
#  ======================== STEP 2a: CLEAN ALT IMMUNOPEPTIDOME ===========
#  Load ALL ALT predictions, remove WT peptides, then tier what remains.
###########################################################################

cat("\n========== STEP 2a: Clean ALT immunopeptidome (by length) ==========\n")

# Process and write one peptide length at a time.  Excluded rows are counted
# but never accumulated: downstream scripts only consume retained tiers.
process_in_chunks <- function(type, apply_wt_exclusion = FALSE) {
  prefix <- paste0(type, "_concordance")
  paths <- c(
    tier1 = file.path(directory_15, paste0(prefix, "_tier1_", current_date, ".tsv")),
    tier2 = file.path(directory_15, paste0(prefix, "_tier2_", current_date, ".tsv")),
    tier3 = file.path(directory_15, paste0(prefix, "_tier3_", current_date, ".tsv")),
    all   = file.path(directory_15, paste0(prefix, "_all_",   current_date, ".tsv"))
  )
  unlink(paths)
  counts <- setNames(integer(4), c("Tier1_HighConfidence",
                                   "Tier2_MediumConfidence",
                                   "Tier3_Discordant_NMPstrong",
                                   "Tier3_Discordant_MHCstrong"))
  excluded <- 0L
  for (len in lengths) {
    cat(sprintf("\n[INFO] Processing %s %smer chunk...\n", toupper(type), len))
    nmp <- load_nmp_raw(len, type)
    mhc <- load_mhc_raw(len, type)
    if (apply_wt_exclusion) {
      nmp <- nmp[!peptide %chin% wt_exclusion]
      mhc <- mhc[!peptide %chin% wt_exclusion]
    }
    tiered <- assign_tiers(join_tools(nmp, mhc))
    tab <- tiered[, .N, by = concordance_tier]
    excluded <- excluded + tab[concordance_tier == "Excluded", sum(N)]
    for (nm in names(counts)) counts[nm] <- counts[nm] + tab[concordance_tier == nm, sum(N)]

    retained <- tiered[concordance_tier != "Excluded"]
    pieces <- list(
      tier1 = retained[concordance_tier == "Tier1_HighConfidence"],
      tier2 = retained[concordance_tier == "Tier2_MediumConfidence"],
      tier3 = retained[grepl("Tier3", concordance_tier)],
      all = retained
    )
    for (nm in names(paths)) {
      if (nrow(pieces[[nm]]) > 0L)
        fwrite(pieces[[nm]], paths[[nm]], sep = "\t", na = "NA",
               append = file.exists(paths[[nm]]), col.names = !file.exists(paths[[nm]]))
    }
    cat(sprintf("[INFO]   retained=%d, excluded=%d\n", nrow(retained),
                tab[concordance_tier == "Excluded", sum(N)]))
    rm(nmp, mhc, tiered, retained, pieces, tab)
    gc(verbose = FALSE)
  }
  list(counts = counts, excluded = excluded, paths = paths)
}

alt_results <- process_in_chunks("alt", apply_wt_exclusion = TRUE)

# Exact compatibility alias used by downstream 15b/c/d.
file.copy(alt_results$paths[["all"]],
          file.path(directory_15, paste0("cross_alg_all_nmers_", current_date, ".tsv")),
          overwrite = TRUE)

###########################################################################
#  ======================== STEP 2b: WT NATIVE IMMUNOPEPTIDOME ===========
#  Tier the WT predictions independently as a reference set.
###########################################################################

cat("\n========== STEP 2b: WT native immunopeptidome (by length) ==========\n")
wt_results <- process_in_chunks("wt", apply_wt_exclusion = FALSE)

###########################################################################
#  Summary
###########################################################################

cat("\n=== 15a Summary ===\n")
cat("\nALT (tumour-specific, WT-filtered):\n")
cat(sprintf("  Tier 1 -- High confidence (NMP SB + MHC <500nM):  %d\n",
            alt_results$counts[["Tier1_HighConfidence"]]))
cat(sprintf("  Tier 2 -- Medium confidence (NMP WB + MHC <500nM): %d\n",
            alt_results$counts[["Tier2_MediumConfidence"]]))
cat(sprintf("  Tier 3 -- Discordant (tools disagree):             %d\n",
            sum(alt_results$counts[grepl("Tier3", names(alt_results$counts))])))

cat("\nWT native immunopeptidome (reference):\n")
cat(sprintf("  Tier 1 -- High confidence:   %d\n", wt_results$counts[["Tier1_HighConfidence"]]))
cat(sprintf("  Tier 2 -- Medium confidence: %d\n", wt_results$counts[["Tier2_MediumConfidence"]]))
cat(sprintf("  Tier 3 -- Discordant:        %d\n",
            sum(wt_results$counts[grepl("Tier3", names(wt_results$counts))])))

cat(sprintf("\n  WT peptides excluded from ALT: %d\n", length(wt_exclusion)))
cat(sprintf("  All files written to: %s\n", directory_15))
cat("\n[DONE] Step 15a complete.\n")
