#!/bin/bash
#SBATCH -J ampliseq_main
#SBATCH -A r01340
#SBATCH -p general
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

# Head job only: it schedules, it does not compute. Each pipeline task is
# submitted as its own slurm job (see conf/quartz_main.config).

source /etc/profile
module load java/21.0.8 apptainer/1.3.6
export NXF_HOME=/N/project/BacInteraction/schapag_cowrumen/nextflow_home
export NXF_SINGULARITY_CACHEDIR=/N/scratch/schapag/apptainer_cache/nxf
export APPTAINER_CACHEDIR=/N/scratch/schapag/apptainer_cache
# raised from -Xmx4g: the interrupted fetchngs run died with a native JVM OOM
# tracking 2,367 accessions (hs_err_pid541634.log). ~7,500 tasks here.
export NXF_OPTS="-Xms2g -Xmx12g"
export NXF_VER=26.04.6

cd /N/scratch/schapag/projectA_vmb_ptb

/N/u/schapag/Quartz/.local/bin/nextflow run nf-core/ampliseq -r 2.18.0 \
    -profile singularity \
    -c conf/quartz_main.config \
    -params-file params/main.yaml \
    -w work/main \
    -resume ampliseq_main
