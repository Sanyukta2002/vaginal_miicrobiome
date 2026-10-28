# Decisions log — Vaginal microbiome & preterm birth: CSTs vs continuous frequencies

Analysis of public raw 16S reads (SRA SRP115697; study: Callahan et al. 2017, PNAS,
doi:10.1073/pnas.1705899114) via nf-core/ampliseq, comparing VALENCIA CSTs with continuous
taxon frequencies. Published results are used only to validate the pipeline. Every parameter choice and its evidence is recorded here.

Project root: `/N/scratch/schapag/projectA_vmb_ptb`

---

## Environment

| Component | Version |
|---|---|
| Nextflow | 26.04.6 |
| nf-core tools | 4.1.0 |
| nf-core/ampliseq | 2.18.0 |
| nf-core/fetchngs | 1.13.0 |
| Apptainer | 1.3.6 |
| Java | 21.0.8 |
| Python | 3.12 |
| Platform | IU Quartz, SLURM, account r01340 |

---

## Decisions

| Step | Decision | Options considered | Evidence / reason | Date |
|---|---|---|---|---|
| 0 | Apptainer for containers | Docker vs Apptainer | No admin privileges on IU HPC; Apptainer is the supported option | Sep 2026 |
| 1 | `nf-core/fetchngs` as downloader | fetchngs vs sra-tools `prefetch`/`fasterq-dump` | Nextflow-native, resumable, emits an ampliseq samplesheet | Sep 2026 |
| 1 | Per-run error models via `run` column | Pooled vs per-run | Matches the paper's method ("error rates estimated for each sequencing run to reduce batch effects", SI); required to resolve 1–2 nt *Gardnerella* variants | Sep 2026 |
| 1 | Sample types resolved from Stanford `processedrds` | SRA metadata vs authors' deposit | SRA lacks body site and run; authors' object gave both, closing the 188-run gap exactly | Sep 2026 |
| 2 | QC via a pilot run of the pipeline, not standalone FastQC | standalone FastQC+MultiQC vs pilot through ampliseq | Standalone QC creates an unversioned in-house stage outside the workflow provenance; and truncation cannot be judged without post-denoising merge rates, which only the pipeline produces | 2026-09-24 |
| 2 | 96 IL02 runs have R1/R2 labels transposed — corrected in samplesheet | rename FASTQs vs fix samplesheet vs drop IL02 | Archive files left untouched; samplesheet points `forwardReads` at `_2`. Verified by merge rate (see Findings 2) | 2026-09-24 |
| 2 | No `--ignore_binned_quality` | set vs unset | 31 distinct Phred characters spanning Q2 and Q10–Q40; binned data shows 4–8. Encoding is continuous, DADA2 error model unaffected | 2026-09-24 |
| 3 | `--skip_cutadapt`; primers NOT trimmed | pass `--FW_primer`/`--RV_primer` vs `--retain_untrimmed` vs `--skip_cutadapt` | Primers absent from reads (0/40,000 matches; reads begin at the first biological base). EMP fusion-primer design sequences from the gene-specific primer onward | 2026-09-24 |
| 3 | `--trunclenf 245`, `--trunclenr 235` (explicit, not auto) | auto via `--trunc_qmin` vs explicit | Values are the authors' own (SI Bioinformatics). Auto-detection would choose per-batch values, producing different ASV lengths per run and fracturing `mergeSequenceTables()` | 2026-09-24 |
| 3 | `--max_ee 2`; `maxN=0` (DADA2 default) | — | SI: "no ambiguous bases allowed, and each read required to have less than two expected errors" | 2026-09-24 |
| 3 | `--ignore_empty_input_files`, `--ignore_failed_trimming`, `--ignore_failed_filtering` all set | omit vs set | Required, not optional: NTCs carry 6–117 reads and would fail post-filter gates, halting the whole run. Controls must stay in for `decontam` (Step 7) | 2026-09-24 |
| 4 | Chimera removal by consensus (ampliseq default) | pooled vs consensus vs per-sample | SI: "ASVs were removed if identified as chimeric in a sufficient fraction of the samples in which they were present" | 2026-09-24 |
| 8 | Wilcoxon p-values with the pre-R 4.6 rule (`03_analysis/functions.R`) | R 4.6 default (exact conditional with ties) vs exact only without ties | R 4.6.0 made `wilcox.test()` exact with ties, which changed every tied Stanford test. The pre-4.6 rule reproduces the published counts exactly (43 variants, 44 genera, Mycoplasma p = 0.042) and makes results independent of R version | 2026-10-28 |
| 5 | `--skip_dada_addspecies` NOT used in the full run | skip (cheaper) vs keep | Species resolution is the paper's central result — low *L. crispatus* associates with PTB in both cohorts, *L. iners* does not; at genus level the signal vanishes. Tune `--dada_assign_chunksize` if memory is tight rather than dropping resolution | 2026-09-24 |
| 5 | SILVA 138 | SILVA 128/132/138, GTDB, sbdi-gtdb (ampliseq default) | Paper used SILVA v123, not offered by ampliseq 2.18.0; 138 is the nearest supported release. Expect species-level differences vs the published table, notably from the 2020 *Lactobacillus* reclassification and the *Gardnerella* genus split — relevant to Step 6 validation | 2026-09-24 |
| pilot | `--skip_taxonomy`, `--skip_qiime` | full pipeline vs trimmed | Taxonomy (SILVA download + `assignTaxonomy`) is the expensive tail and irrelevant to parameter selection. Restored for the full run | 2026-09-24 |
| 7 | rarefy or not | McMurdie & Holmes 2014 vs Schloss 2024 | *pending* | |

---

## Findings

### 1. Amplicon architecture — why no primer trimming

Primers from the SI fusion-primer constructs (gene-specific portions only):
FW 515F-Y GTGYCAGCMGCCGCGGTAA 19 nt
RV 806RB GGACTACNVGGGTWTCTAAT 20 nt

This is the Parada/Apprill updated EMP V4 pair, **not** the older Caporaso 515F/806R —
note `Y` at forward position 4 and `NV` in the reverse. Using the original sequences
would mis-trim.

Neither primer appears in the reads. A 40,000-read scan returned 0 matches, and reads
begin `TACGTAGGTGGCAAGCGTT…` — the first biological base after 515F. The EMP design's
sequencing primer terminates at the gene-specific primer, so the primer is never read.

Arithmetic, which the ASV length distribution later confirmed:

- E. coli 515→806 spans 292 nt **including** both primers
- primer-free insert = 292 − 19 − 20 = **253 bp**
- reads are a uniform 250 nt → neither read spans the insert → no adapter read-through
- at 245 + 235 = 480, overlap = 480 − 253 = **227 bp**, against `mergePairs()`'s minimum of 12
- binding constraint is therefore `truncLenF + truncLenR ≥ 265`; the authors used 480

Truncation here is driven purely by quality, not by overlap.

### 2. IL02 R1/R2 transposition (96 runs)

**What.** 96 of IL02's 182 runs were deposited with R1/R2 labels swapped: the `_1` file
holds reverse-oriented reads, the `_2` file forward-oriented. No other run is affected.
The 182 runs are 182 distinct sample aliases — not duplicate submissions.

**How detected.** Forward V4 reads begin `TAC…`; reverse reads begin `CCT…`
(revcomp of the conserved region inside the 806R site). Scanning the first read of all
2,367 R1 files split the dataset cleanly into two groups, all 96 outliers in IL02.
Confirmed by deep check: flagged files are ~99% reverse in R1 and ~99% forward in R2,
the exact mirror of a normal file.

**Two false starts, both classifier artifacts rather than data problems** — worth
recording because the first two answers were wrong:

1. Five samples appeared to form a third category. The motif pattern required `G` at
   position 4 (`^[TN]ACG`), but forward reads legitimately start `TACA`, `TACT`, `TGCA`.
   All five were normal.
2. Twenty-two samples appeared possibly mixed-orientation — both orientations inside one
   file — which would have been serious and unfixable by a samplesheet edit. Those samples
   are dominated by a taxon reading `CCTATT` rather than `CCTGTT`, one base off the
   pattern. Relaxing the discriminator to `TAC` vs `CCT` resolved all of them.

Lesson: anchor on the conserved part of a motif, not a literal string.

**Fix.** FASTQ files untouched. `bin/build_samplesheet.py` reads
`assets/reverse_oriented_runs.txt` and writes `forwardReads` = `_2`, `reverseReads` = `_1`
for those accessions, with guards asserting that every listed accession matched a row,
that the applied count equals the list length, and that all swapped samples are IL02.
`assets/metadata_ampliseq.tsv` carries a `MatesSwapped` column so the correction travels
with the data.

**Verification (pilot, 2026-09-24).** 10 corrected + 10 normal IL02 swabs were forced
into the pilot:

| Group | n | merged % (min / median / max) |
|---|---|---|
| IL02 R1/R2-corrected | 10 | 81.7 / **86.9** / 90.1 |
| IL02 normal | 10 | 77.8 / **84.6** / 90.1 |
| All other runs | 107 | 68.6 / **86.4** / 96.1 |

Indistinguishable. Had the correction been wrong, merge rates would have collapsed to
near zero. Left uncorrected, IL02 would have silently failed to merge and been easy to
misdiagnose as a truncation problem.

### 3. Dataset composition and depth

2,367 runs, 79.5 GB, reconciling exactly against the paper:

| Type | n |
|---|---|
| VaginalSwab | 2,179 |
| EC (extraction control) | 101 |
| NTC (no-template control) | 61 |
| Mock community | 26 |

Cohort and sequencing run are perfectly confounded (Stanford: IL02, IL03, IL05, IL07,
IL08; UAB: IL06, IL09). This is tolerable because the paper never tests Stanford vs UAB —
every Wilcoxon runs *within* a cohort, so the batch effect is constant across the
PTB/term comparison and cannot manufacture an association. The cohorts meet only at the
Fisher combination, where a signal must survive independently in two batches, two
populations and two risk profiles.

Depth (ENA `read_count`): 119 runs <1k, 39 at 1–5k, 17 at 5–10k, 2,192 ≥10k. The
shallowest runs are all EC/NTC — which is what healthy negative controls should look
like, and good news for `decontam`.

**Outlier flagged for Step 7:** `SRR5970816`, an extraction control on IL03, carries
**90,157 reads** — three orders of magnitude above every other EC. Candidate explanations:
genuine heavy contamination in that extraction batch, index hopping / well-to-well
spillover, or a sample-label swap. Force-included in the pilot via `FORCE_SRR` so it is
characterised early. Directly relevant to the paper's own Fig. 4 contamination analysis.

Mock communities exist in IL02, IL05, IL06, IL07, IL08 (26 total). This upgrades Step 6:
measured false-positive ASV and chimera rates against known ground truth are available,
not just sequence matching to the authors' table.

### 4. Pilot validation (2026-09-24, job 10604867)

166 samples (127 swabs stratified across all 7 runs with a floor of 10, plus 15 EC,
14 NTC, 10 mocks). 535 tasks, all succeeded. 47 min wall, 9.5 CPU-hours, peak RSS
12.5 GB of 64 GB requested.

- merge rate, swabs: median **86.4%**, range 68.6–96.1%; **zero** swabs below 50%
- chimeric: median **3.5%**, max 17.9% (normal band is 5–25%)
- ASV lengths: 1,543 at 253 bp, 513 at 252, 88 at 254 — **98.3% within 252–254**.
  Off-target negligible (23 at 273, 3 at 315). Independently confirms the 253 bp
  primer-free arithmetic.
- 2,182 ASVs total
- per-run `DADA2_ERR` / `DENOISING` / `RMCHIMERA` / `STATS` all report 7 of 7,
  confirming the `run` column partitions correctly

**Per-run merge rates reveal a real compromise:**

| Run | n | median merged % |
|---|---|---|
| IL02 | 20 | 86.4 |
| IL03 | 15 | 80.1 |
| IL05 | 19 | 93.3 |
| IL06 | 26 | 83.7 |
| IL07 | 13 | **73.3** |
| IL08 | 10 | 84.7 |
| IL09 | 24 | 93.5 |

IL07 is 20 points below IL05/IL09, most likely from earlier reverse-read quality decay.
73% is acceptable, but it is direct evidence that one global truncation is a compromise
across runs — and the question the sensitivity analysis should answer is whether a
shorter reverse truncation rescues IL07 without perturbing community profiles elsewhere.
A single truncation remains mandatory regardless, since per-run values would fracture the
merged ASV table.

Controls merge at 71–93%, so their reads are real community sequence rather than noise —
further material for `decontam`.

---

### 5. Nextflow head-JVM memory for the full run

The interrupted fetchngs run left `hs_err_pid541634.log`: `Out of Memory Error
(arena.cpp:168)` — a native JVM allocation failure in the head process, not a task
failure, running under `NXF_OPTS="-Xms1g -Xmx4g"` while tracking 2,367 accessions.

The pilot (166 samples, 535 tasks) ran fine at the same setting. The full run is
~7,500–8,000 tasks, roughly 15x the head-process bookkeeping. Raise to
`NXF_OPTS="-Xms2g -Xmx12g"` for the full run; the compute node has 64 GB and the head
JVM competes with little else. Cheap insurance against losing a multi-hour run.

Related: the Step 0 test profile saw two QIIME2 tasks exit 137 at a 1 GB request,
recovered on automatic retry — a test-profile resource request, not a `quartz.config` cap.

---

## Open items

- [ ] Truncation sensitivity analysis: rerun the pilot at ~240/220 and ~220/180; compare
      reads retained, merge rate, and Bray–Curtis concordance (Mantel) against 245/235.
      Not industry-standard practice, but this project's stated purpose is documented
      judgment on exactly this parameter — and IL07 gives it real content.
- [ ] Vaginal swabs under 10k reads: count and decide exclusions
- [ ] Characterise `SRR5970816` (90k-read extraction control)
- [x] Versions recorded in `docs/versions.md` (authoritative; Environment table above is a summary)
- [ ] R environment on HPC (phyloseq, vegan, picante, decontam) — needed from Step 7/9
- [ ] Rarefy or not (Step 7)

### 6. Concurrent Nextflow runs and `-resume`

Job 10610317 (240/220 sweep) failed 14 s in: `Unable to acquire lock on session`.
Bare `-resume` resolves to the most recent session in the launch directory, so two
jobs submitted together resolved to the *same* session ID and contended on one
LevelDB cache lock.

Sweep runs therefore drop `-resume`: with different parameters no cached task can be
reused, so it offers nothing and actively prevents concurrency. `-resume` is retained
for the full run, where interruption recovery matters — but two resumable runs must
never be launched concurrently from the same directory.

### 7. Truncation sensitivity analysis (2026-09-24)

The same 166-sample pilot was rerun at two alternative truncations and compared to
the authors' 245/235. ASV sequences differ between settings, so tables are not
directly comparable; instead a Bray-Curtis distance matrix was built independently
per setting on the shared samples and compared to the reference by Mantel test
(999 permutations). This asks whether the parameter changes the *relationships
between samples* — the thing the biology depends on.

| setting | ASVs | median merged% | min merged% | on-target reads | Mantel r | p |
|---|---|---|---|---|---|---|
| 245/235 (paper) | 2,182 | 86.4 | 68.6 | 99.96% | — | — |
| 240/220 | 2,267 | 92.2 | 78.9 | 99.96% | 0.9999 | 0.0010 |
| 220/180 | 2,282 | 94.9 | 88.8 | 99.95% | 0.9997 | 0.0010 |

p = 0.0010 is the floor for 999 permutations, i.e. as significant as the test allows.

**Community structure is invariant to truncation.** r ≈ 1.0 in both comparisons, and
the on-target read fraction (ASVs of 252–254 bp) does not move. The ~100 extra ASVs
at shorter truncation are short off-target products carrying ~0.01% of reads — most
likely short amplicons whose 250 nt reads run into adapter, which a shorter truncation
trims away so they merge at their true length. ASV *count* is a misleading summary here;
reads are not moving.

**Read retention is not invariant.** Median merge rises monotonically, and every run
improves:

| run | 245/235 | 240/220 | 220/180 |
|---|---|---|---|
| IL02 | 86.4 | 91.7 | 94.6 |
| IL03 | 80.1 | 90.3 | 94.3 |
| IL05 | 93.3 | 94.7 | 96.4 |
| IL06 | 83.7 | 91.1 | 94.2 |
| IL07 | **73.3** | 80.9 | **93.0** |
| IL08 | 84.7 | 91.9 | 94.8 |
| IL09 | 93.5 | 95.1 | 96.9 |

IL07 gains 20 points. Note the loss at 245/235 is run-structured, and Stanford spans
five runs (IL02/03/05/07/08) — so within that cohort an IL07 sample loses far more
reads than an IL05 sample. If PTB and term women are not balanced across runs within a
cohort, that is a technical gradient sitting inside the comparison the paper makes.
220/180 largely flattens it (93.0–96.9% vs 73.3–93.5%).

**Caveat, stated rather than glossed:** part of the merge gain is mechanical, not pure
improvement. Shorter reads accumulate fewer expected errors so fewer are rejected by
`maxEE 2`, and a shorter overlap has fewer opportunities to violate `mergePairs`'s
zero-mismatch default. More merging is partly less stringency. That ASV count rose only
4% and on-target read fraction did not move indicates the leniency is not importing
noise, but the mechanism should be acknowledged.

**DECISION: retain 245/235 for the primary analysis.** Using the original study's parameters
means the pipeline can be validated against its published table: any divergence is
attributable to software versions rather than to our choices, which is what Step 6
checks. The sensitivity analysis is reported as a finding:
the parameter everyone picks arbitrarily does not move the biology here (r = 0.9997),
and we can now quantify what the authors' choice cost in read retention. A supplementary
220/180 run is justified if depth becomes limiting.

Script: `bin/compare_truncations.py`. Output: `results/truncation_sensitivity.txt`.
