#!/bin/bash

# =============================================================================
# RNA-seq Analysis Pipeline Script for HPC
# =============================================================================
# This script performs a complete RNA-seq analysis workflow:
# 1. Prepares reference indices for STAR and RSEM.
# 2. Iterates through all paired-end FASTQ files in a given directory.
# 3. Aligns reads to the genome using STAR.
# 4. Performs post-alignment QC using Samtools, QualiMap, and RSeQC.
# 5. Quantifies gene expression using RSEM.
#
# Designed to be run in the background using nohup.
# =============================================================================

# --- Script Behavior ---
# Exit immediately if a command exits with a non-zero status.
set -e
# Treat unset variables as an error when substituting.
set -u
# The return value of a pipeline is the status of the last command
# to exit with a non-zero status.
set -o pipefail

# =============================================================================
# START OF NEW CODE BLOCK TO ADD
# =============================================================================
# --- Activate Conda Environment ---
# IMPORTANT: Source the conda initialization script first.
# REPLACE the path below with the one you found in Step 1.
source /home/haoxm_pkuhpc/haoxm_cls/miniconda3/etc/profile.d/conda.sh

# Now, activate the specific environment for this pipeline
conda activate rna_seq_env
# --- End of Activation ---
# =============================================================================
# END OF NEW CODE BLOCK
# =============================================================================

# --- Configuration Section: PLEASE EDIT PATHS AND PARAMETERS HERE ---

# Number of threads to use for multithreaded tasks
THREADS=24

# Read length of your FASTQ files (e.g., 150 for PE150)
READ_LENGTH=150

# Main working directory (assuming you run the script from here)
# It's better to use absolute paths for robustness in HPC environments
MAIN_DIR="/home/haoxm_pkuhpc/haoxm_cls/gpfs1/zhw" # <-- IMPORTANT: REPLACE WITH THE ABSOLUTE PATH to your zhw directory

# Input/Output Directories based on your structure
FASTQ_DIR="${MAIN_DIR}/fastp_results"
OUTPUT_DIR="${MAIN_DIR}/genome_alignment"
REF_DIR="${OUTPUT_DIR}/reference"

# Reference Genome and Annotation Files
GENOME_FA_GZ="${REF_DIR}/GRCh38.primary_assembly.genome.fa.gz"
ANNOTATION_GTF_GZ="${REF_DIR}/gencode.v49.primary_assembly.basic.annotation.gtf.gz"

# --- End of Configuration ---

# --- Derived File Paths ---
# Remove .gz suffix for decompressed filenames
GENOME_FA="${GENOME_FA_GZ%.gz}"
ANNOTATION_GTF="${ANNOTATION_GTF_GZ%.gz}"
ANNOTATION_BED="${REF_DIR}/gencode.v49.primary_assembly.basic.annotation.bed"

# Derived Index Directories
STAR_INDEX_DIR="${REF_DIR}/STAR_index_basic"
RSEM_REF_DIR="${REF_DIR}/RSEM_ref_basic"


# =============================================================================
# STEP 0: PREPARE REFERENCE INDICES (Done only once)
# =============================================================================
echo "================================================="
echo "STEP 0: PREPARING REFERENCE INDICES"
echo "================================================="

# Decompress reference files if they are still gzipped
if [ -f "$GENOME_FA_GZ" ]; then
    echo "Decompressing genome FASTA..."
    gunzip "$GENOME_FA_GZ"
fi
if [ -f "$ANNOTATION_GTF_GZ" ]; then
    echo "Decompressing annotation GTF..."
    gunzip "$ANNOTATION_GTF_GZ"
fi

# 0.1: Build STAR Index
# --genomeSAindexNbases指构建的索引长度，默认14，建议取10-15。该值越大会消耗越多的内存，但是检索的更快。
# 但是对于小基因组来说，不能太大，如果索引太长就会造成索引总数少的问题，可以通过(log2(GenomeLength)/2 - 1)计算得到
if [ ! -d "$STAR_INDEX_DIR" ]; then
    echo "STAR index not found. Building..."
    mkdir -p "$STAR_INDEX_DIR"
    STAR --runMode genomeGenerate \
         --runThreadN ${THREADS} \
         --genomeDir ${STAR_INDEX_DIR} \
         --genomeFastaFiles ${GENOME_FA} \
         --sjdbGTFfile ${ANNOTATION_GTF} \
         --sjdbOverhang $((READ_LENGTH - 1))
else
    echo "STAR index found. Skipping build."
fi

# 0.2: Prepare RSEM Reference
if [ ! -d "$RSEM_REF_DIR" ]; then
    echo "RSEM reference not found. Preparing..."
    # Manually create the output directory for RSEM, as it does not create it automatically.
    mkdir -p "$RSEM_REF_DIR"
    rsem-prepare-reference --gtf ${ANNOTATION_GTF} \
                           --star \
                           -p ${THREADS} \
                           ${GENOME_FA} \
                           ${RSEM_REF_DIR}/human_ref
else
    echo "RSEM reference found. Skipping preparation."
fi

# --- 0.3: Create a clean BED file for RSeQC using the official UCSC toolchain ---
# This is the definitive, robust method to convert GTF to a RSeQC-compatible BED file.

# Define intermediate filename
GENEPRED_FILE="${REF_DIR}/gencode.v49.primary_assembly.basic.annotation.genePred"

if [ ! -f "$ANNOTATION_BED" ]; then
    echo "RSeQC-compatible BED file not found. Creating from GTF..."
    
    # Step 1: Convert GTF to genePred format using UCSC's validated tool.
    echo "    -> Step 1/2: Converting GTF to genePred..."
    gtfToGenePred -genePredExt -ignoreGroupsWithoutExons ${ANNOTATION_GTF} ${GENEPRED_FILE}

    # Step 2: Convert the genePred file to a perfectly formatted BED12 file.
    echo "    -> Step 2/2: Converting genePred to BED..."
    genePredToBed ${GENEPRED_FILE} ${ANNOTATION_BED}
    
    # Clean up the intermediate genePred file
    rm -f ${GENEPRED_FILE}

    echo "--- BED file for RSeQC created successfully ---"
else
    echo "RSeQC-compatible BED file found. Skipping creation."
fi

echo "--- Reference preparation complete ---"


# =============================================================================
# MAIN LOOP: PROCESS EACH SAMPLE
# =============================================================================
echo "================================================="
echo "STARTING MAIN ANALYSIS LOOP"
echo "================================================="

for R1_FILE in ${FASTQ_DIR}/*_clean_1.fq.gz; do
    
    # --- A. Setup file names and directories for the current sample ---
    BASENAME=$(basename "$R1_FILE")
    SAMPLE_PREFIX=${BASENAME%_clean_1.fq.gz}
    R2_FILE="${FASTQ_DIR}/${SAMPLE_PREFIX}_clean_2.fq.gz"

    echo "-------------------------------------------------"
    echo "Processing Sample: ${SAMPLE_PREFIX}"
    echo "-------------------------------------------------"

    # Define and create output directories for this sample
    STAR_OUT_DIR="${OUTPUT_DIR}/01_star_alignment/${SAMPLE_PREFIX}"
    QC_OUT_DIR="${OUTPUT_DIR}/02_post_alignment_qc/${SAMPLE_PREFIX}"
    RSEM_OUT_DIR="${OUTPUT_DIR}/03_rsem_quantification/${SAMPLE_PREFIX}"
    
    mkdir -p ${STAR_OUT_DIR} ${QC_OUT_DIR} ${RSEM_OUT_DIR}

    # Define major output filenames
    BAM_SORTED_BY_COORD="${STAR_OUT_DIR}/${SAMPLE_PREFIX}_Aligned.sortedByCoord.out.bam"
    BAM_TO_TRANSCRIPTOME="${STAR_OUT_DIR}/${SAMPLE_PREFIX}_Aligned.toTranscriptome.out.bam"

    # --- B. STEP 1: STAR Alignment ---
    echo "--> Step 1: Running STAR alignment..."
    STAR --genomeDir ${STAR_INDEX_DIR} \
         --runThreadN ${THREADS} \
         --readFilesIn ${R1_FILE} ${R2_FILE} \
         --readFilesCommand zcat \
         --outFileNamePrefix "${STAR_OUT_DIR}/${SAMPLE_PREFIX}_" \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --quantMode TranscriptomeSAM \
         --twopassMode Basic

    # --- C. STEP 2: Post-Alignment QC ---
    echo "--> Step 2: Running Post-Alignment QC..."
    
    # 2.1 Samtools index and flagstat
    echo "    -> Samtools index & flagstat..."
    samtools index ${BAM_SORTED_BY_COORD}
    samtools flagstat ${BAM_SORTED_BY_COORD} > "${QC_OUT_DIR}/${SAMPLE_PREFIX}_flagstat.txt"
    
    # 2.2 QualiMap
    echo "    -> QualiMap rnaseq..."
    qualimap rnaseq -bam ${BAM_SORTED_BY_COORD} \
                     -gtf ${ANNOTATION_GTF} \
                     -outdir "${QC_OUT_DIR}/qualimap" \
                     -pe \
                     --java-mem-size=8G

    # 2.3 RSeQC
    echo "    -> RSeQC read_distribution..."
    read_distribution.py -i ${BAM_SORTED_BY_COORD} -r ${ANNOTATION_BED} > "${QC_OUT_DIR}/${SAMPLE_PREFIX}_read_distribution.txt"
    
    echo "    -> RSeQC infer_experiment (checking strandedness)..."
    infer_experiment.py -i ${BAM_SORTED_BY_COORD} -r ${ANNOTATION_BED} > "${QC_OUT_DIR}/${SAMPLE_PREFIX}_strand_specificity.txt"

    # --- D. STEP 3: RSEM Quantification ---
    echo "--> Step 3: Running RSEM quantification..."
    rsem-calculate-expression --paired-end \
                              --bam \
                              --strandedness none \
                              -p ${THREADS} \
                              ${BAM_TO_TRANSCRIPTOME} \
                              ${RSEM_REF_DIR}/human_ref \
                              "${RSEM_OUT_DIR}/${SAMPLE_PREFIX}"

    echo "--- Finished processing sample ${SAMPLE_PREFIX} ---"

done

echo "================================================="
echo "ALL SAMPLES PROCESSED. PIPELINE COMPLETE."
echo "================================================="
