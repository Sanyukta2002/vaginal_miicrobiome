#!/bin/bash
#SBATCH -J fetchngs
#SBATCH -A r01340
#SBATCH -p general
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
# Download all SRP115697 FASTQs with nf-core/fetchngs. Submit from 01_pipeline/.

cd "$SLURM_SUBMIT_DIR"
source conf/env_quartz.sh

nextflow run nf-core/fetchngs -r 1.13.0 \
    -profile singularity \
    -c conf/fetchngs_local.config \
    -w work/fetchngs \
    --input assets/ids.csv \
    --outdir raw/fetchngs \
    --download_method ftp \
    --nf_core_pipeline ampliseq \
    -resume
