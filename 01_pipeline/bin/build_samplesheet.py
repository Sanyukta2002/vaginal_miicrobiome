#!/usr/bin/env python3
"""Build the ampliseq samplesheets and metadata from the ENA run table and the study's sample tables.

Runs listed in reverse_oriented_runs.txt (96, all IL02) were deposited with R1/R2 transposed;
the samplesheet swaps them back.
"""
import csv, glob, os, re, sys, collections

PROJ = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
ENA  = os.path.join(PROJ, "assets/SRP115697_ena_runs.tsv")
PS   = os.path.join(PROJ, "assets/callahan_allsamples_run_bodysite.tsv")
META = os.path.join(PROJ, "assets/callahan_sample_metadata.tsv")
SWAP = os.path.join(PROJ, "assets/reverse_oriented_runs.txt")
FQ   = os.path.join(PROJ, "raw/fetchngs/fastq")

def read_tsv(p):
    with open(p, newline="") as f:
        rows = []
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append({(k.strip() if k else k): (v.strip() if isinstance(v, str) else v)
                         for k, v in r.items()})
        return rows

if not os.path.exists(SWAP):
    sys.exit(f"ERROR: swap list not found: {SWAP}")
with open(SWAP) as f:
    swap_ids = {ln.strip() for ln in f if ln.strip()}
print(f"swap list: {len(swap_ids)} accessions")

# FASTQ pairs by run accession
files = {}
for mate in (1, 2):
    for path in glob.glob(f"{FQ}/*_{mate}.fastq.gz"):
        m = re.search(r"(SRR\d+)", os.path.basename(path))
        if m:
            files.setdefault(m.group(1), {})[mate] = path

# Run and body site for every sample; clinical metadata for the 2,179 swabs
ps = {r["SampleID"]: r for r in read_tsv(PS)}
meta = {r["SampleID"]: r for r in read_tsv(META)}

rows, missing, applied = [], [], set()
for r in read_tsv(ENA):
    srr, alias = r["run_accession"], r["sample_alias"]
    if srr not in files or 1 not in files[srr] or 2 not in files[srr]:
        missing.append(srr); continue

    fwd, rev = files[srr][1], files[srr][2]
    swapped = srr in swap_ids
    if swapped:
        fwd, rev = rev, fwd
        applied.add(srr)

    info  = ps.get(alias, {})
    body  = info.get("BodySite", "NA")
    stype = "Mock" if body.startswith("mock") else body
    rows.append({
        "sampleID":     "S" + alias,
        "forwardReads": fwd,
        "reverseReads": rev,
        "run":          info.get("Run", "NA"),
        "_alias": alias, "_srr": srr, "_type": stype,
        "_swapped": "yes" if swapped else "no",
    })

# Stop on missing files, unassigned runs, or swaps outside IL02
if missing:
    sys.exit(f"ERROR: {len(missing)} runs have no FASTQ pair, e.g. {missing[:5]}")

unapplied = swap_ids - applied
if unapplied:
    sys.exit(f"ERROR: {len(unapplied)} swap-list accessions never matched a row: "
             f"{sorted(unapplied)[:5]}")
assert len(applied) == len(swap_ids), "swap count mismatch"

bad_run = [r["_srr"] for r in rows if r["run"] in ("NA", "")]
if bad_run:
    sys.exit(f"ERROR: {len(bad_run)} rows have no sequencing run assigned "
             f"(would create a phantom batch): {bad_run[:5]}")

swap_runs = {r["run"] for r in rows if r["_swapped"] == "yes"}
if swap_runs != {"IL02"}:
    sys.exit(f"ERROR: swapped samples span unexpected runs: {sorted(swap_runs)}")

def write_sheet(path, subset):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sampleID", "forwardReads", "reverseReads", "run"])
        for r in subset:
            w.writerow([r["sampleID"], r["forwardReads"], r["reverseReads"], r["run"]])
    print(f"{len(subset):5d}  {os.path.relpath(path, PROJ)}")

write_sheet(os.path.join(PROJ, "assets/samplesheet_all.csv"), rows)
write_sheet(os.path.join(PROJ, "assets/samplesheet_swabs.csv"),
            [r for r in rows if r["_type"] == "VaginalSwab"])

cols = ["SubjectID","Cohort","Race","GW","GDcoll","Outcome","preterm","GWdel","DeliveryMode"]
with open(os.path.join(PROJ, "assets/metadata_ampliseq.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["ID","SRR","SampleType","Run","MatesSwapped"] + cols)
    for r in rows:
        m = meta.get(r["_alias"], {})
        w.writerow([r["sampleID"], r["_srr"], r["_type"], r["run"], r["_swapped"]] +
                   [m.get(c, "NA") for c in cols])

print("\nby sample type:")
for k, v in collections.Counter(r["_type"] for r in rows).most_common():
    print(f"  {k:12s} {v}")
print("\nby run:")
for k, v in sorted(collections.Counter(r["run"] for r in rows).items()):
    print(f"  {k:6s} {v}")
print(f"\nR1/R2 corrected: {len(applied)} samples (all IL02)")
