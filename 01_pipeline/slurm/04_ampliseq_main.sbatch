#!/bin/bash
#SBATCH -J ampliseq_main
#SBATCH -A r01340
#SBATCH -p general
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
# Full run, all 2,367 samples. Head job only: each task is its own SLURM job (conf/quartz_main.config).
# Submit from 01_pipeline/.

cd "$SLURM_SUBMIT_DIR"
source conf/env_quartz.sh
export NXF_OPTS="-Xms2g -Xmx12g"

nextflow run nf-core/ampliseq -r 2.18.0 \
    -profile singularity \
    -c conf/quartz_main.config \
    -params-file params/main.yaml \
    -w work/main \
    -resume ampliseq_main
