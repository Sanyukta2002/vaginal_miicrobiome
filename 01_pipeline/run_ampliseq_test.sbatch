#!/bin/bash
#SBATCH -J ampliseq_test
#SBATCH -A r01340
#SBATCH -p general
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
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
    -profile test,singularity \
    -c conf/quartz.config \
    -w work/ampliseq_test \
    --outdir results/ampliseq_test
