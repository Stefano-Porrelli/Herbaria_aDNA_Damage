#!/bin/bash
#SBATCH --job-name=Plant_aDNA_pipeline_prep
#SBATCH --output=Plant_aDNA_pipeline_prep_%j.out
#SBATCH --error=Plant_aDNA_pipeline_prep_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

# Preparation: create conda environment, install programs and Dependencies,
# Retrieve reference genome and index it
# Change variables below to project name and FTP link to reference genome assembly

PROJECT_NAME="Oryza_Sediment_KALIMANTAN"
REFERENCE_GENOME_FTP="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.fna.gz"

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
log "Starting aDNA pipeline preparation..."

# Check if the conda environment already exists, if not create it
if ! conda env list | grep -q "Plant_aDNA_pipeline"; then
    log "Plant_aDNA_pipeline environment not found. Creating it now."
    # Create YAML file for conda environment
    cat << 'EOF' > Plant_aDNA_pipeline_env.yml
name: Plant_aDNA_pipeline
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - adapterremoval
  - bwa
  - fastqc
  - samtools
  - dedup
  - mapdamage2
  - multiqc
  - preseq
  - r-base
  - r-essentials
  - r-gam
  - r-rcppgsl
  - r-rcpp
  - r-inline
  - r-ggplot2
  - numpy
  - matplotlib
  - pysam
EOF
    check_command "Creating YAML file"
    
    # Create conda environment from YAML file
    log "Creating conda environment (this may take several minutes)..."
    conda env create -f Plant_aDNA_pipeline_env.yml
    check_command "Creating conda environment"
    
    log "Conda environment created successfully."
else
    log "Plant_aDNA_pipeline environment already exists. Skipping creation."
fi

# Activate the environment
log "Activating conda environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate Plant_aDNA_pipeline
check_command "Activating conda environment"

# Create project directory and subdirectories
log "Setting up project directory structure..."
PROJECT_DIR="$(pwd)/${PROJECT_NAME}"
mkdir -p "$PROJECT_DIR"/{1_initial_data/reference_genome,2_trimmed_merged,3_quality_control,4_mapping,5_aDNA_characteristics,6_AMBER,7_preseq}
check_command "Creating project directories"

# Download and prepare reference genome
cd "${PROJECT_DIR}/1_initial_data/reference_genome"
log "Downloading reference genome..."
curl -O "${REFERENCE_GENOME_FTP}"
check_command "Downloading reference genome"

# Unzip
log "Unzipping reference genome..."
gunzip *.gz
check_command "Unzipping reference genome"

# Set reference genome variable to the unzipped file
REFERENCE_GENOME=$(ls *.fna)
log "Reference genome file: ${REFERENCE_GENOME}"

# Index reference genome
log "Indexing reference genome with BWA..."
bwa index "${REFERENCE_GENOME}"
check_command "Indexing reference genome with BWA"

log "Indexing reference genome with Samtools..."
samtools faidx "${REFERENCE_GENOME}"
check_command "Indexing reference genome with Samtools"

log "Setup complete! Conda environment and directory structure created, reference genome indexed."
log "Reference genome file: ${PROJECT_DIR}/1_initial_data/reference_genome/${REFERENCE_GENOME}"
log "Project directory: ${PROJECT_DIR}"
