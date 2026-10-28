#!/bin/bash
#SBATCH -A r01340
#SBATCH -p general
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
# Pilot at another truncation setting. Usage (from 01_pipeline/):
#   sbatch -J amp_<TF>_<TR> slurm/03_ampliseq_truncation.sbatch <trunclenf> <trunclenr>

TF=$1
TR=$2
if [ -z "$TF" ] || [ -z "$TR" ]; then
  echo "ERROR: need two args: trunclenf trunclenr" >&2; exit 1
fi
if [ $((TF + TR)) -lt 265 ]; then
  echo "ERROR: TF+TR=$((TF+TR)) < 265; insufficient overlap for a 253 bp insert" >&2; exit 1
fi

cd "$SLURM_SUBMIT_DIR"
source conf/env_quartz.sh

nextflow run nf-core/ampliseq -r 2.18.0 \
    -profile singularity \
    -c conf/ampliseq_local.config \
    -w work/ampliseq_t${TF}_${TR} \
    --input assets/samplesheet_pilot.csv \
    --outdir results/pilot_t${TF}_${TR} \
    --skip_cutadapt \
    --trunclenf ${TF} \
    --trunclenr ${TR} \
    --max_ee 2 \
    --skip_taxonomy \
    --skip_qiime \
    --ignore_empty_input_files \
    --ignore_failed_trimming \
    --ignore_failed_filtering \
    -name sweep_t${TF}_${TR}_$(date +%s)
