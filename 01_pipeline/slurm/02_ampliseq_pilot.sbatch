#!/bin/bash
#SBATCH -J ampliseq_pilot
#SBATCH -A r01340
#SBATCH -p general
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
# Pilot run (166 samples) at the study's DADA2 settings. Submit from 01_pipeline/.

cd "$SLURM_SUBMIT_DIR"
source conf/env_quartz.sh

nextflow run nf-core/ampliseq -r 2.18.0 \
    -profile singularity \
    -c conf/ampliseq_local.config \
    -w work/ampliseq_pilot \
    --input assets/samplesheet_pilot.csv \
    --outdir results/pilot \
    --skip_cutadapt \
    --trunclenf 245 \
    --trunclenr 235 \
    --max_ee 2 \
    --skip_taxonomy \
    --skip_qiime \
    --ignore_empty_input_files \
    --ignore_failed_trimming \
    --ignore_failed_filtering \
    -resume
