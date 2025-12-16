#!/bin/bash
#SBATCH --job-name=Plant_aDNA_pipeline
#SBATCH --output=Plant_aDNA_pipeline_%j.out
#SBATCH --error=Plant_aDNA_pipeline_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

################################################################################
########   This is a demonstration of the pipeline described in:        ########
########                                                                ########
########  Latorre S.M., Lang, P.L.M., Burbano, H.A., Gutaker, R.M. 2020 ########
########   Isolation, Library Preparation, and Bioinformatic Analysis   ########
########               of Historical and Ancient Plant DNA              ########
########           d.o.i.: https://doi.org/10.1002/cppb.20121           ########
########					version 4.0 - 14/11/2025					########
################################################################################

# Main aDNA pipeline - processing samples
# Make sure to update PROJECT_NAME to match your prep script

PROJECT_NAME="ENTER_DIRECTORY_NAME"

# Function to check if a command was successful
check_command() {
    if [ $? -ne 0 ]; then
        echo "Error: $1 failed"
        exit 1
    fi
}
# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}
log "Starting aDNA pipeline processing..."

# Activate the environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate Plant_aDNA_pipeline
check_command "Activating conda environment"

# Create project directory and subdirectories
PROJECT_DIR="$(pwd)/${PROJECT_NAME}"
INITIAL_DATA_DIR="${PROJECT_DIR}/1_initial_data"
REFERENCE_GENOME_DIR="${PROJECT_DIR}/1_initial_data/reference_genome"

# Find the reference genome file (should be the .fna file from prep script)
REFERENCE_GENOME=$(ls ${REFERENCE_GENOME_DIR}/*.fna)
if [[ -z "$REFERENCE_GENOME" ]]; then
    log "Error: No reference genome (.fna file) found in ${REFERENCE_GENOME_DIR}"
    exit 1
fi
log "Using reference genome: $(basename $REFERENCE_GENOME)"

# Check if project directory exists
if [[ ! -d "$PROJECT_DIR" ]]; then
    log "Error: Project directory $PROJECT_DIR not found. Run the preparation script first."
    exit 1
fi

# aDNA_screening main pipeline:
cd "${PROJECT_DIR}"

# Find all fastq.gz files in the initial data directory
sample_count=0
for fastq_file in ${INITIAL_DATA_DIR}/*_1.fastq.gz; do
    # Skip if no files match
    [ -e "$fastq_file" ] || continue

    # Extract sample_id from the filename
    sample_id=$(basename "$fastq_file" _1.fastq.gz)
    
    # Check if paired file exists
    if [[ ! -f "${INITIAL_DATA_DIR}/${sample_id}_2.fastq.gz" ]]; then
        log "Skipping $sample_id: Paired file not found"
        continue
    fi

    ((sample_count++))
    log "Processing sample $sample_count: $sample_id"

    # Trim adapters and merge paired reads
    log "  Running AdapterRemoval..."
    AdapterRemoval --file1 "${INITIAL_DATA_DIR}/${sample_id}_1.fastq.gz" \
                   --file2 "${INITIAL_DATA_DIR}/${sample_id}_2.fastq.gz" \
                   --basename "2_trimmed_merged/${sample_id}" \
                   --collapse --gzip --threads 16
    check_command "AdapterRemoval for $sample_id"

    # Remove the "collapsed" tag - useful for MultiQC downstream
    mv "2_trimmed_merged/${sample_id}.collapsed.gz" "2_trimmed_merged/${sample_id}.gz"
    check_command "Renaming collapsed file for $sample_id"

    # FastQC
    log "  Running FastQC..."
    fastqc -t 16 -o 3_quality_control/ \
           "2_trimmed_merged/${sample_id}.gz" \
           "2_trimmed_merged/${sample_id}.pair1.truncated.gz" \
           "2_trimmed_merged/${sample_id}.pair2.truncated.gz"
    check_command "FastQC for $sample_id"

    # Map reads to indexed reference 
    log "  Running BWA alignment..."
    bwa aln -t 16 -l 1024 -f "4_mapping/${sample_id}.collapsed.sai" \
        "${REFERENCE_GENOME}" \
        "2_trimmed_merged/${sample_id}.gz"
    check_command "BWA alignment for $sample_id"

    # Convert to SAM
    log "  Converting to SAM..."
    bwa samse -r "@RG\\tID:${sample_id}\\tSM:${sample_id}" \
              -f "4_mapping/${sample_id}.sam" \
              "${REFERENCE_GENOME}" \
              "4_mapping/${sample_id}.collapsed.sai" \
              "2_trimmed_merged/${sample_id}.gz"
    check_command "SAM conversion for $sample_id"

    # Calculate proportion of mapped reads
    log "  Calculating mapping statistics..."
    samtools flagstat "4_mapping/${sample_id}.sam" > "4_mapping/${sample_id}.flagstat.log"
    check_command "Samtools flagstat for $sample_id"

    # Create compressed BAM from the mapped reads
    log "  Creating BAM file..."
    samtools view -@ 16 -F 4 -Sbh -o "4_mapping/${sample_id}.mapped.bam" "4_mapping/${sample_id}.sam"
    check_command "BAM creation for $sample_id"

    # Sort BAM by chromosome 
    log "  Sorting BAM file..."
    samtools sort -@ 16 -o "4_mapping/${sample_id}.mapped.sorted.bam" "4_mapping/${sample_id}.mapped.bam"
    check_command "BAM sorting for $sample_id"

    # Remove SAM and unsorted BAM
    rm "4_mapping/${sample_id}.sam" "4_mapping/${sample_id}.mapped.bam"
    check_command "Removing intermediate files for $sample_id"

    # Filter PCR duplicates
    log "  Removing PCR duplicates..."
    dedup -i "4_mapping/${sample_id}.mapped.sorted.bam" -m -o 4_mapping/
    check_command "Deduplication for $sample_id"

    # Assess DNA damage
    log "  Assessing DNA damage..."
    mapDamage -i "4_mapping/${sample_id}.mapped.sorted_rmdup.bam" \
              -r "${REFERENCE_GENOME}" \
              --no-stats -y 0.05 -d "5_aDNA_characteristics/${sample_id}"
    check_command "mapDamage for $sample_id"

    # Run AMBER to assess mapping biases
    log "  Assessing mapping biases with AMBER..."
    
    # Copy and sort BAM for AMBER (AMBER needs properly sorted BAMs)
    samtools sort -@ 16 -o "6_AMBER/${sample_id}.mapped.sorted_rmdup_sorted.bam" "4_mapping/${sample_id}.mapped.sorted_rmdup.bam"
    check_command "Sorting BAM for AMBER for $sample_id"
    
    # Index the BAM file (AMBER may require this)
    samtools index "6_AMBER/${sample_id}.mapped.sorted_rmdup_sorted.bam"
    check_command "Indexing BAM for AMBER for $sample_id"
    
    # Download AMBER if not already present
    if [[ ! -f "${PROJECT_DIR}/6_AMBER/AMBER" ]]; then
        cd "${PROJECT_DIR}/6_AMBER"
        curl -L -o AMBER https://github.com/tvandervalk/AMBER/raw/main/AMBER
        chmod +x ./AMBER
        check_command "Downloading AMBER"
    fi
    
    # Create BAMs.txt file for this sample
    cd "${PROJECT_DIR}/6_AMBER"
    echo -e "${sample_id}\t./${sample_id}.mapped.sorted_rmdup_sorted.bam" > "${sample_id}_BAMs.txt"
    # Process with AMBER - but don't fail the entire pipeline if AMBER fails
    if ./AMBER --bamfiles "${sample_id}_BAMs.txt" --errorbars --counts --output "${sample_id}_AMBER"; then
        log "  AMBER analysis completed successfully for $sample_id"
    else
        log "  Warning: AMBER analysis failed for $sample_id - continuing with other samples"
        log "  AMBER error likely due to insufficient data, possibly screening blanks?"
    fi
    # Clean up temporary files
    rm "${sample_id}_BAMs.txt"
    cd "${PROJECT_DIR}"
    
    # Run preseq analyses
    # preseq c_curve
    log "  Running preseq c_curve complexity curve..."
    if preseq c_curve -v -s 50 -o "7_preseq/${sample_id}_complexity.txt" -B "4_mapping/${sample_id}.mapped.sorted.bam"; then
        log "  Preseq complexity curve completed successfully for $sample_id"
    else
        log "  Warning: Preseq c_curve complexity curve failed for $sample_id - continuing with other analyses"
        log "  This is often due to insufficient unique reads in the sample, try reduce -s value "
    fi

    # preseq lc_extrap
    log "  Running preseq library complexity extrapolation..."
    if preseq lc_extrap -v -s 500 -D -o "7_preseq/${sample_id}_future_yields.txt" -B "4_mapping/${sample_id}.mapped.sorted.bam"; then
        log "  Preseq lc_extrap library complexity extrapolation completed successfully for $sample_id"
    else
        log "  Warning: Preseq lc_extrap library complexity extrapolation failed for $sample_id - continuing with other analyses"
        log "  This is often due to insufficient duplicate reads for extrapolation, try reduce -s value"
    fi

    # preseq gc_extrap
    log "  Running preseq genomic coverage extrapolation..."
    # Convert BAM to mapped reads format
    if to-mr -o "7_preseq/${sample_id}.mapped.mr" -L 10000 "4_mapping/${sample_id}.mapped.sorted.bam"; then
        # Sort the mapped reads file
        if sort -k 1,1 -k 2,2n -k 3,3n "7_preseq/${sample_id}.mapped.mr" > "7_preseq/${sample_id}.mapped.sorted.mr"; then
            # Run genomic coverage extrapolation
            cd "${PROJECT_DIR}/7_preseq"
            if preseq gc_extrap -v -defects -s 500 -o "${sample_id}_future_coverage.txt" "${sample_id}.mapped.sorted.mr"; then
                log "  Preseq genomic coverage extrapolation completed successfully for $sample_id"
            else
                log "  Warning: Preseq genomic coverage extrapolation failed for $sample_id - continuing with other analyses"
                log "  This is often due to insufficient coverage data for extrapolation"
            fi
            cd "${PROJECT_DIR}"
        else
            log "  Warning: Failed to sort mapped reads file for $sample_id - skipping genomic coverage extrapolation"
        fi
    else
        log "  Warning: Failed to convert BAM to mapped reads format for $sample_id - skipping genomic coverage extrapolation"
    fi
    
    log "Finished processing sample: $sample_id"
done

if [[ $sample_count -eq 0 ]]; then
    log "No paired FastQ files found in ${INITIAL_DATA_DIR}"
    log "Expected files: *_1.fastq.gz and *_2.fastq.gz"
    exit 1
fi

log "All $sample_count samples processed successfully."

# Run MultiQC to collate data
log "Running MultiQC to generate summary report..."
cd "${PROJECT_DIR}" 
multiqc 1_initial_data/ 2_trimmed_merged/ 3_quality_control/ 4_mapping/ 5_aDNA_characteristics/ 6_AMBER/ 7_preseq/ \
        --cl-config "read_count_multiplier: 0.001" --cl-config "read_count_prefix: K" --cl-config "read_count_desc: thousands"
check_command "MultiQC report generation"

log "Pipeline complete! Check the MultiQC report for results summary."
