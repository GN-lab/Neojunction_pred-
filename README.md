# Neojunction_prediction

Neojunction prediction from 307 transcriptomic dataset acquired from Hartwig (Breast sample)

## 1. Neojunction calling
The characterization of neojunctions is first performed by evaluating splice sites in tumor samples derived from TCGA and filtering against splice sites with putative read expression in a normal tissue repository, GTEx (n=9166).

### Step 01. Tumor purity
From meta data files, remove samples with a tumor purity of lower than 60%. Various TCGA publications will include meta data quantifying tumor purity. For our studies, we have used Kahles et al. 2018 (Cancer Cell).

### Step 02. Protein-coding genes
Identify and filter for protein-coding genes with a corresponding GTF file (Homo_sapiens.GRCh37.87.chr.gtf in our study). This will generate a dataframe of protein-encoding genes and their respectie coordinates.

### Step 03. TPM filter > 10
From tumor samples with a tumor purity > 0.60, identify transcripts with TPM values > 10.

### Step 04. Characterize annotated splicing junctions
Next characterize all splicing junctions that are identified at alignment. From STAR aligner's SJ.out.tab files, the splice junction site count data can be extracted. Based on Step 02 and Step 03, select for splice sites that are found in protein-coding genes with a TPM filter > 10.

### Step 05. Extract putative non-annotated splice junctions
Filter out annotated splicing junctions in the corresponding GTF file (sjdbList.fromGTF.out.tab) to identify non-annotated splicing junctions. Note that in this step, different GTF versions will affect the characterization of neojunctions that are identified. At the time of this code being used, our GRCh37.87 GTF sj.out.tab file version was GENCODE v33. 

### Step 06. Prepare splicing junction overlap table


**Goal** (Summary)

- Start with STAR sj.out.tab files.
- Remove junctions annotated in GENCODE(GRch38).
- Keep only junctions overlapping non-mitochondrial, protein-coding genes.
- Remove junctions with <10 spliced reads (per sample) or <20 total reads (cohort).
- Compute spliced frequency:
- Frequency = total target spliced reads / (target + canonical spliced reads).
- Retain junctions with frequency >1%.
- Define “public” junctions:
- Expressed in ≥10% of cohort, with above criteria.
- Remove junctions found in >1% of GTEx normal samples.
- Output: high-confidence, cancer-specific, non-annotated junctions.

**STAR aligner** 

STAR produces SJ.out.tab: Contains every splice junction STAR detected, with per-junction read counts (unique + multi-mapping reads) and annotation flags.

- srun nextflow run nf-core/rnaseq \
 
  --aligner star_salmon \
  --skip_quantification \
  --skip_qc \
   --star_align_args "--outSAMtype BAM SortedByCoordinate \
                     --quantMode TranscriptomeSAM \
                     --outSAMunmapped Within \
                     --outSJfilterOverhangMin 15 20 20 20 \
                     --alignSJoverhangMin 8 \
                     --alignIntronMin 20 \
                     --alignIntronMax 1000000" \
  -profile singularity \
  -c /data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/config_files/icr_alma.config \
  -r 3.18.0 \
  -resume

--outSJfilterOverhangMin 15 20 20 20
Controls the minimum overhang (length of mapped sequence flanking a splice junction) for different splice junction types to be reported. This increases stringency on junction detection:
  -  15 for canonical GT/AG junctions
  -  20 for other types (e.g., GC/AG, AT/AC, non-canonical)
This helps reduce false positives and ensures only well-supported junctions are output.

--alignSJoverhangMin 8
Sets the minimum overhang length required on each side of a splice junction to consider it a valid splice junction during alignment. A value of 8 bases means STAR must see at least 8 matching bases flanking the junction.

--alignIntronMin 20
Specifies the minimum allowed intron length (20 bases) for splice junctions. Very short introns are usually sequencing or alignment artifacts, so this filters those out.

--alignIntronMax 1000000
Sets the maximum allowed intron length (1,000,000 bases). This limits extremely long introns that could be biologically implausible or mapping artifacts.

<img width="2400" height="1800" alt="figure_5i_fs_if_boxplot_20250928" src="https://github.com/user-attachments/assets/afada8f3-e0b4-4ad9-8843-6276ad1510eb" />

<img width="3000" height="1800" alt="figure_5i_splice_types_jitter_20250928" src="https://github.com/user-attachments/assets/cdb3c880-395a-4753-9251-12ddf68b5f1f" />

<img width="2400" height="1500" alt="figure_5i_fs_if_density_all_20250928" src="https://github.com/user-attachments/assets/1b15fbf6-babc-4fb4-bcd7-caaca6ac60d3" />

<img width="2400" height="3000" alt="figure_5i_splice_types_density_all_20250928" src="https://github.com/user-attachments/assets/73be694c-fd42-4bdf-89d6-a7c3983c3f9f" />

<img width="2400" height="1800" alt="figure_5i_fs_if_jitter_20250928" src="https://github.com/user-attachments/assets/40a758b0-f2d4-464c-b5f4-51e530a5b98b" />


