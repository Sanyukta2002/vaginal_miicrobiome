#!/bin/bash
#SBATCH -J ampliseq_test
#SBATCH -A r01340
#SBATCH -p general
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
# nf-core/ampliseq test profile on the cluster. Submit from 01_pipeline/.

cd "$SLURM_SUBMIT_DIR"
source conf/env_quartz.sh

nextflow run nf-core/ampliseq -r 2.18.0 \
    -profile test,singularity \
    -c conf/quartz.config \
    -w work/ampliseq_test \
    --outdir results/ampliseq_test
