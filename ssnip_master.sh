#!/bin/bash
#SBATCH --job-name=ssnip_driver
#SBATCH --partition=compute
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=8042
#SBATCH --time=48:00:00
#SBATCH --output=logs/master_%j.log

set -euo pipefail
set -E
trap 'rc=$?; echo "**ERROR** rc=$rc at line $LINENO cmd: $BASH_COMMAND" >&2; exit $rc' ERR

# Source config directly
WORKDIR="$(pwd)"
source "${WORKDIR}/config.sh"

# Generate samples.txt if missing
if [[ ! -s "${WORKDIR}/samples.txt" ]]; then
  echo "[INFO] Generating samples.txt from STAR SJ.out.tab files..."
  find "${STAR_SJ_DIR}" -type f -name "*.SJ.out.tab" | \
    awk -F'/' '{fname=$NF; sub(/\.SJ\.out\.tab$/, "", fname); print fname}' | \
    sort -u > "${WORKDIR}/samples.txt"
  echo "[INFO] Created samples.txt with $(wc -l < "${WORKDIR}/samples.txt") samples."
fi

export SAMPLES_TXT="${WORKDIR}/samples.txt"

THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_NUM_THREADS="$THREADS"

# Load modules and activate environments
module load R  # Load R for Rscript (adjust version if needed)
module load bedtools/2.29.2
source "../venv/bin/activate"  # Update if venv path differs
source /data/scratch/DMP/UCEC/EVOLIMMU/csalas/miniconda3/etc/profile.d/conda.sh

echo "=== Submitting SSNIP Pipeline (dependency-driven, individual scripts) ==="

# Pre-submission check: Ensure key inputs exist
for f in "${GTF_FILE}" "${PURITY_FILE}" "${TPM_FILE}" "${STAR_SJ_DIR}" "${GTEX_FILE}"; do
  [[ -e "$f" ]] || { echo "[ERROR] Missing required file/directory: $f" >&2; exit 1; }
done

# Helper: Check if step is done via checkpoint
step_done() { [[ -f "${WORKDIR}/.checkpoints/$1.done" ]]; }

# Rolling dependency: the last job ID
jobPrev=""

# Function to submit a job and get ID
submit_job() {
  local name="$1"
  local command="$2"
  local dep="${3:-}"
  local jobid=""
  if [[ -n "$dep" ]]; then
    jobid=$(sbatch --parsable --job-name="$name" --partition=compute --cpus-per-task=4 --mem-per-cpu=8042 --time=12:00:00 --output=logs/${name}_%j.log --dependency=afterok:${dep} --wrap="$command; touch ${WORKDIR}/.checkpoints/${name}.done")
  else
    jobid=$(sbatch --parsable --job-name="$name" --partition=compute --cpus-per-task=4 --mem-per-cpu=8042 --time=12:00:00 --output=logs/${name}_%j.log --wrap="$command; touch ${WORKDIR}/.checkpoints/${name}.done")
  fi
  echo "Submitted $name as $jobid"
  return $jobid
}

# ---- Neojunction Calling ----
if ! step_done "01_tumor_purity"; then
  jobPrev=$(submit_job "01_tumor_purity" "Rscript ${NEOJUNCTION_DIR}/01_tumor_purity.R")
else echo "01_tumor_purity already done."; fi

if ! step_done "02_protein_coding_genes"; then
  jobPrev=$(submit_job "02_protein_coding_genes" "Rscript ${NEOJUNCTION_DIR}/02_protein_coding_genes.R" "$jobPrev")
else echo "02_protein_coding_genes already done."; fi

if ! step_done "03_tpm_filter_10"; then
  jobPrev=$(submit_job "03_tpm_filter_10" "Rscript ${NEOJUNCTION_DIR}/03_tpm_filter_10.R" "$jobPrev")
else echo "03_tpm_filter_10 already done."; fi

if ! step_done "04_extract_annot_sj_to_analyze"; then
  jobPrev=$(submit_job "04_extract_annot_sj_to_analyze" "Rscript ${NEOJUNCTION_DIR}/04_extract_annot_sj_to_analyze.R" "$jobPrev")
else echo "04_extract_annot_sj_to_analyze already done."; fi

if ! step_done "05a_extract_sj_nonannot_tcga.count.min10"; then
  jobPrev=$(submit_job "05a_extract_sj_nonannot_tcga.count.min10" "Rscript ${NEOJUNCTION_DIR}/05a_extract_sj_nonannot_tcga.count.min10.R" "$jobPrev")
else echo "05a_extract_sj_nonannot_tcga.count.min10 already done."; fi

if ! step_done "05b_extract_sj_nonannot_tcga.count.min10_protein.coding_psr10"; then
  jobPrev=$(submit_job "05b_extract_sj_nonannot_tcga.count.min10_protein.coding_psr10" "Rscript ${NEOJUNCTION_DIR}/05b_extract_sj_nonannot_tcga.count.min10_protein.coding_psr10.R" "$jobPrev")
else echo "05b_extract_sj_nonannot_tcga.count.min10_protein.coding_psr10 already done."; fi

if ! step_done "06_prep_overlap_table"; then
  jobPrev=$(submit_job "06_prep_overlap_table" "Rscript ${NEOJUNCTION_DIR}/06_prep_overlap_table.R" "$jobPrev")
else echo "06_prep_overlap_table already done."; fi

if ! step_done "07_count_depth_freq_judge"; then
  jobPrev=$(submit_job "07_count_depth_freq_judge" "Rscript ${NEOJUNCTION_DIR}/07_count_depth_freq_judge.R" "$jobPrev")
else echo "07_count_depth_freq_judge already done."; fi

if ! step_done "08_overlap_table_gtex"; then
  jobPrev=$(submit_job "08_overlap_table_gtex" "Rscript ${NEOJUNCTION_DIR}/08_overlap_table_gtex.R" "$jobPrev")
else echo "08_overlap_table_gtex already done."; fi

if ! step_done "09_count_depth_freq_judge_psr_gtex"; then
  jobPrev=$(submit_job "09_count_depth_freq_judge_psr_gtex" "Rscript ${NEOJUNCTION_DIR}/09_count_depth_freq_judge_psr_gtex.R" "$jobPrev")
else echo "09_count_depth_freq_judge_psr_gtex already done."; fi

if ! step_done "10_extract_neojunctions"; then
  jobPrev=$(submit_job "10_extract_neojunctions" "Rscript ${NEOJUNCTION_DIR}/10_extract_neojunctions.R" "$jobPrev")
else echo "10_extract_neojunctions already done."; fi

# ---- Neopeptide Prediction ----
if ! step_done "11_aaseq_prediction"; then
  jobPrev=$(submit_job "11_aaseq_prediction" "Rscript ${NEOPEPTIDE_DIR}/11_aaseq_prediction.R" "$jobPrev")
else echo "11_aaseq_prediction already done."; fi

if ! step_done "12_nmer_generation"; then
  jobPrev=$(submit_job "12_nmer_generation" "Rscript ${NEOPEPTIDE_DIR}/12_nmer_generation.R" "$jobPrev")
else echo "12_nmer_generation already done."; fi

# ---- Presentation Prediction ----

if ! step_done "14a_mhcflurry2_input_df_generation"; then
  jobPrev=$(submit_job "14a_mhcflurry2_input_df_generation" "Rscript ${PRESENTATION_DIR}/14a_mhcflurry2_input_df_generation.R" "$jobPrev")
else echo "14a_mhcflurry2_input_df_generation already done."; fi

if ! step_done "14b_mhcflurry2_analysis_with_flank"; then
  jobPrev=$(submit_job "14b_mhcflurry2_analysis_with_flank" "bash ${PRESENTATION_DIR}/14b_mhcflurry2_analysis_with_flank.sh" "$jobPrev")
else echo "14b_mhcflurry2_analysis_with_flank already done."; fi

if ! step_done "14c_select_top_alleles"; then
  jobPrev=$(submit_job "14c_select_top_alleles" "Rscript ${PRESENTATION_DIR}/14c_select_top_alleles.R" "$jobPrev")
else echo "14c_select_top_alleles already done."; fi

if ! step_done "14d_generate_figures"; then
  jobPrev=$(submit_job "14d_generate_figures" "Rscript ${PRESENTATION_DIR}/14d_generate_figures.R" "$jobPrev")
else echo "14d_generate_figures already done."; fi

if ! step_done "15d_scores_for_nj_types"; then
  jobPrev=$(submit_job "15d_scores_for_nj_types" "Rscript ${PRESENTATION_DIR}/15d_scores_for_nj_types.R" "$jobPrev")
else echo "15d_scores_for_nj_types already done."; fi

echo "=== Submission complete. Jobs chained with afterok dependencies (where applicable). ==="
