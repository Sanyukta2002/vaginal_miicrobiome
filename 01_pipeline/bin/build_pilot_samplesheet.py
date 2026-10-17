#!/usr/bin/env python3
"""Stratified pilot samplesheet for parameter selection (Steps 2-4).

~127 vaginal swabs across all 7 runs (floor 10/run), plus EC/NTC/mock controls.
IL02 is deliberately loaded with both R1/R2-swapped and normal samples so the
orientation correction is validated by the merge rate.
Deterministic: fixed seed.
"""
import csv, os, random, sys, collections

PROJ  = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
SHEET = os.path.join(PROJ, "assets/samplesheet_all.csv")
META  = os.path.join(PROJ, "assets/metadata_ampliseq.tsv")
ENA   = os.path.join(PROJ, "assets/SRP115697_ena_runs.tsv")
OUT   = os.path.join(PROJ, "assets/samplesheet_pilot.csv")
OUTM  = os.path.join(PROJ, "assets/metadata_pilot.tsv")

SEED         = 20260924
BASE_PER_RUN = 10      # floor of swabs per run
EXTRA_TOTAL  = 50      # distributed proportionally on top of the floor
IL02_SWAPPED = 10      # forced swapped IL02 swabs
IL02_NORMAL  = 10      # forced normal  IL02 swabs
MIN_DEPTH    = 10000   # swabs only; controls are exempt
EC_PER_RUN   = 2
NTC_PER_RUN  = 2
MOCK_PER_RUN = 2
FORCE_SRR    = {"SRR5970816"}   # IL03 extraction control, 90,157 reads — extreme outlier

random.seed(SEED)

def read_csv(p, delim=","):
    with open(p, newline="") as f:
        return [{k.strip(): (v.strip() if isinstance(v, str) else v)
                 for k, v in r.items()} for r in csv.DictReader(f, delimiter=delim)]

sheet = {r["sampleID"]: r for r in read_csv(SHEET)}
meta  = read_csv(META, "\t")
depth = {r["run_accession"]: int(r["read_count"] or 0) for r in read_csv(ENA, "\t")}

for m in meta:
    m["_depth"] = depth.get(m["SRR"], 0)

runs = sorted({m["Run"] for m in meta})
swabs = [m for m in meta if m["SampleType"] == "VaginalSwab"]
swab_n = collections.Counter(m["Run"] for m in swabs)
total_swabs = sum(swab_n.values())

# --- allocation per run -------------------------------------------------
alloc = {}
for r in runs:
    alloc[r] = BASE_PER_RUN + round(EXTRA_TOTAL * swab_n[r] / total_swabs)
alloc["IL02"] = IL02_SWAPPED + IL02_NORMAL

def pick(pool, n):
    pool = list(pool)
    random.shuffle(pool)
    return pool[:n]

chosen = []

# IL02: forced split so the swap fix is exercised
il02 = [m for m in swabs if m["Run"] == "IL02" and m["_depth"] >= MIN_DEPTH]
sw  = [m for m in il02 if m["MatesSwapped"] == "yes"]
nsw = [m for m in il02 if m["MatesSwapped"] == "no"]
if len(sw) < IL02_SWAPPED or len(nsw) < IL02_NORMAL:
    sys.exit(f"ERROR: not enough IL02 swabs above {MIN_DEPTH} reads "
             f"(swapped={len(sw)}, normal={len(nsw)})")
chosen += pick(sw, IL02_SWAPPED) + pick(nsw, IL02_NORMAL)

# other runs: swabs above depth floor
for r in runs:
    if r == "IL02":
        continue
    pool = [m for m in swabs if m["Run"] == r and m["_depth"] >= MIN_DEPTH]
    want = min(alloc[r], len(pool))
    if want < alloc[r]:
        print(f"  note: {r} only has {len(pool)} swabs above {MIN_DEPTH} reads "
              f"(wanted {alloc[r]})")
    chosen += pick(pool, want)

# controls and mocks: NO depth filter (near-empty NTCs are the point)
for r in runs:
    for stype, k in (("EC", EC_PER_RUN), ("NTC", NTC_PER_RUN), ("Mock", MOCK_PER_RUN)):
        pool = [m for m in meta if m["Run"] == r and m["SampleType"] == stype]
        chosen += pick(pool, min(k, len(pool)))

# forced inclusions (outliers we want eyes on early)
have = {m["ID"] for m in chosen}
for m in meta:
    if m["SRR"] in FORCE_SRR and m["ID"] not in have:
        chosen.append(m)
        print(f"  forced: {m['ID']} ({m['SRR']}, {m['SampleType']}, "
              f"{m['Run']}, {m['_depth']:,} reads)")

# --- write --------------------------------------------------------------
ids = [m["ID"] for m in chosen]
assert len(ids) == len(set(ids)), "duplicate sample in pilot"

with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sampleID", "forwardReads", "reverseReads", "run"])
    for i in ids:
        r = sheet[i]
        w.writerow([r["sampleID"], r["forwardReads"], r["reverseReads"], r["run"]])

keys = list(meta[0].keys())
keys.remove("_depth")
with open(OUTM, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(keys)
    for m in chosen:
        w.writerow([m[k] for k in keys])

# --- report -------------------------------------------------------------
print(f"\npilot: {len(chosen)} samples -> {os.path.relpath(OUT, PROJ)}")
print("\nby run x type:")
tab = collections.Counter((m["Run"], m["SampleType"]) for m in chosen)
types = ["VaginalSwab", "EC", "NTC", "Mock"]
print(f"  {'run':6s} " + "".join(f"{t:>13s}" for t in types) + f"{'total':>8s}")
for r in runs:
    row = [tab[(r, t)] for t in types]
    print(f"  {r:6s} " + "".join(f"{v:>13d}" for v in row) + f"{sum(row):>8d}")
print(f"  {'ALL':6s} " + "".join(f"{sum(tab[(r,t)] for r in runs):>13d}" for t in types)
      + f"{len(chosen):>8d}")

nsw = sum(1 for m in chosen if m["MatesSwapped"] == "yes")
print(f"\nR1/R2-corrected samples in pilot: {nsw}")
d = [m["_depth"] for m in chosen if m["SampleType"] == "VaginalSwab"]
print(f"swab depth: min={min(d):,}  median={sorted(d)[len(d)//2]:,}  max={max(d):,}")
