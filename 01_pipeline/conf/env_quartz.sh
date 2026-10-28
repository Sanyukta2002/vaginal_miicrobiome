# Environment for IU Quartz (sourced by every slurm/ script)
source /etc/profile
module load java/21.0.8 apptainer/1.3.6
export NXF_HOME=/N/project/BacInteraction/schapag_cowrumen/nextflow_home
export NXF_SINGULARITY_CACHEDIR=/N/scratch/schapag/apptainer_cache/nxf
export APPTAINER_CACHEDIR=/N/scratch/schapag/apptainer_cache
export NXF_OPTS="-Xms1g -Xmx4g"
export NXF_VER=26.04.6
