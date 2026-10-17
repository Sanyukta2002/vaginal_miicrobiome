#!/bin/bash
#SBATCH -J amplipilot
#SBATCH -A r01340
#SBATCH -p general
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
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

/N/u/schapag/Quartz/.local/bin/nextflow run nf-core/ampliseq -r 2.18.0 \
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
