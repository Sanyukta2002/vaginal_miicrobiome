#!/usr/bin/env python3
"""Build nf-core/ampliseq samplesheet + metadata from ENA table and authors' files."""
import csv, glob, os, re, sys, collections

PROJ = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
ENA  = os.path.join(PROJ, "assets/SRP115697_ena_runs.tsv")
PS   = os.path.join(PROJ, "assets/callahan_allsamples_run_bodysite.tsv")
META = os.path.join(PROJ, "assets/callahan_sample_metadata.tsv")
FQ   = os.path.join(PROJ, "raw/fetchngs/fastq")

def read_tsv(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))

# 1. index FASTQ files by run accession
files = {}
for mate in (1, 2):
    for path in glob.glob(f"{FQ}/*_{mate}.fastq.gz"):
        m = re.search(r"(SRR\d+)", os.path.basename(path))
        if m:
            files.setdefault(m.group(1), {})[mate] = path

# 2. authors' Run / BodySite for every deposited sample
ps = {r["SampleID"]: r for r in read_tsv(PS)}
# 3. clinical metadata for the 2179 analysed swabs
meta = {r["SampleID"]: r for r in read_tsv(META)}

rows, missing = [], []
for r in read_tsv(ENA):
    srr, alias = r["run_accession"], r["sample_alias"]
    if srr not in files or 1 not in files[srr] or 2 not in files[srr]:
        missing.append(srr); continue
    info = ps.get(alias, {})
    body = info.get("BodySite", "NA")
    stype = "Mock" if body.startswith("mock") else body
    rows.append({
        "sampleID":     "S" + alias,
        "forwardReads": files[srr][1],
        "reverseReads": files[srr][2],
        "run":          info.get("Run", "NA"),
        "_alias": alias, "_srr": srr, "_type": stype,
    })

if missing:
    sys.exit(f"ERROR: {len(missing)} runs have no FASTQ pair, e.g. {missing[:5]}")

# 4. write samplesheets
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

# 5. metadata file keyed to sampleID
cols = ["SubjectID","Cohort","Race","GW","GDcoll","Outcome","preterm","GWdel","DeliveryMode"]
with open(os.path.join(PROJ, "assets/metadata_ampliseq.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["ID","SRR","SampleType","Run"] + cols)
    for r in rows:
        m = meta.get(r["_alias"], {})
        w.writerow([r["sampleID"], r["_srr"], r["_type"], r["run"]] +
                   [m.get(c, "NA") for c in cols])

print("\nby sample type:")
for k, v in collections.Counter(r["_type"] for r in rows).most_common():
    print(f"  {k:12s} {v}")
print("\nby run:")
for k, v in sorted(collections.Counter(r["run"] for r in rows).items()):
    print(f"  {k:6s} {v}")
