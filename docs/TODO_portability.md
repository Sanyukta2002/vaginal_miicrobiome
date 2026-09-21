# Portability cleanup (do at Step 13)
- Replace hard-coded project path (cd /N/scratch/schapag/...) with $SLURM_SUBMIT_DIR / relative paths
- Move SLURM account (r01340) and partition into conf/quartz.config only
- Move NXF_HOME / cache exports into conf/env.sh (user edits once)
- Call `nextflow` from PATH; keep NXF_VER pinned
- Move run parameters into params/*.yaml
- Keep: pinned pipeline versions (-r), NXF_VER, container versions
- Final check: fresh clone -> run from SRP115697 -> compare outputs
