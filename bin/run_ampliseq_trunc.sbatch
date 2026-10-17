#!/bin/bash
# Usage: sbatch -J amp_<TF>_<TR> bin/run_ampliseq_trunc.sbatch <trunclenf> <trunclenr>
#SBATCH -A r01340
#SBATCH -p general
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

TF=$1
TR=$2
if [ -z "$TF" ] || [ -z "$TR" ]; then
  echo "ERROR: need two args: trunclenf trunclenr" >&2; exit 1
fi
if [ $((TF + TR)) -lt 265 ]; then
  echo "ERROR: TF+TR=$((TF+TR)) < 265; insufficient overlap for a 253 bp insert" >&2; exit 1
fi

source /etc/profile
module load java/21.0.8 apptainer/1.3.6
export NXF_HOME=/N/project/BacInteraction/schapag_cowrumen/nextflow_home
export NXF_SINGULARITY_CACHEDIR=/N/scratch/schapag/apptainer_cache/nxf
export APPTAINER_CACHEDIR=/N/scratch/schapag/apptainer_cache
export NXF_OPTS="-Xms1g -Xmx4g"
export NXF_VER=26.04.6

cd /N/scratch/schapag/projectA_vmb_ptb

echo "truncation sweep: trunclenf=$TF trunclenr=$TR  (overlap $((TF + TR - 253)) bp)"

/N/u/schapag/Quartz/.local/bin/nextflow run nf-core/ampliseq -r 2.18.0 \
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
