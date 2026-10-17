#!/bin/bash
#SBATCH -J fetchngs
#SBATCH -A r01340
#SBATCH -p general
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

source /etc/profile
module load java/21.0.8 apptainer/1.3.6
export NXF_HOME=/N/project/BacInteraction/schapag_cowrumen/nextflow_home
export NXF_SINGULARITY_CACHEDIR=/N/scratch/schapag/apptainer_cache/nxf
export APPTAINER_CACHEDIR=/N/scratch/schapag/apptainer_cache
export NXF_OPTS="-Xms1g -Xmx4g"
export NXF_VER=26.04.6

cd /N/scratch/schapag/projectA_vmb_ptb

/N/u/schapag/Quartz/.local/bin/nextflow run nf-core/fetchngs -r 1.13.0 \
    -profile singularity \
    -c conf/fetchngs_local.config \
    -w work/fetchngs \
    --input assets/ids.csv \
    --outdir raw/fetchngs \
    --download_method ftp \
    --nf_core_pipeline ampliseq \
    -resume
