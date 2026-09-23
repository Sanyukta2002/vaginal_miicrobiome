# Tool versions (Step 0)
- Nextflow: 26.04.6 (build 12646)
- nf-core/tools: 4.1.0
- Apptainer: 1.3.6
- Java: 21.0.8
- nf-core/ampliseq: 2.18.0
- nf-core/fetchngs: 1.13.0
- Cluster: IU Quartz, SLURM account r01340
- Project root: /N/scratch/schapag/projectA_vmb_ptb
- Step 0 checkpoint: ampliseq 2.18.0 -profile test,singularity passed on Quartz (2026-09-21), 39 min, 125 succeeded.
  2 tasks (QIIME2_EXPORT_RELTAX, QIIME2_FEATURETABLE_GROUP) exited 137 at 1 GB (test-profile request), succeeded on automatic retry. Not caused by quartz.config caps.

## Step 1 - metadata and samplesheet (2026-09-23)
- Clinical/run metadata from authors' Stanford Digital Repository files (purl.stanford.edu/yb681vm1809):
  input/processed.rda (df: 2,179 swabs x 61 fields) and input/ps_allsams.rds (6,088 samples: SampleID, BodySite, Run)
- Join key: authors' SampleID == ENA sample_alias. 2,179/2,179 matched, 0 unmatched.
- SRP115697 composition: 2,179 vaginal swabs + 101 EC (extraction controls) + 61 NTC + 26 mock communities = 2,367
- Cohort counts reproduce the paper exactly: Stanford 39 women / 897 swabs (157 PTB, 740 term);
  UAB 96 women / 1,282 swabs (392 PTB, 890 term)
- Sequencing runs: IL02, IL03, IL05-IL09 (IL04 had no vaginal swabs, not deposited)
- NOTE: cohort is perfectly confounded with run (Stanford = IL02/03/05/07/08; UAB = IL06/09).
  All association tests are within-cohort; between-cohort contrasts are descriptive only.
- Forward primer (from df): 515F-Y GTGYCAGCMGCCGCGGTAA. Reverse (806R) still to confirm from paper SI.
- Outputs: assets/samplesheet_all.csv (2,367), assets/samplesheet_swabs.csv (2,179),
  assets/metadata_ampliseq.tsv. Sample IDs prefixed "S" (ampliseq requires non-numeric start).
