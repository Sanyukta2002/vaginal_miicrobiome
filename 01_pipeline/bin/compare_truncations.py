#!/usr/bin/env python3
"""Truncation sensitivity: merge rates per setting, and Mantel test of Bray-Curtis matrices against 245/235.

Run from 01_pipeline/ after the pilot and truncation runs.
"""
import csv, itertools, os
import numpy as np

REF   = ("results/pilot",            "245/235 (paper)")
ALTS  = [("results/pilot_t240_220",  "240/220"),
         ("results/pilot_t220_180",  "220/180")]
NPERM = 999
RNG   = np.random.default_rng(20260924)

def load_table(d):
    p = os.path.join(d, "dada2/ASV_table.tsv")
    if not os.path.exists(p):
        return None, None, None
    rows = list(csv.reader(open(p), delimiter="\t"))
    samples = rows[0][1:]
    ids = [r[0] for r in rows[1:]]
    mat = np.array([[float(x) for x in r[1:]] for r in rows[1:]], float)
    return samples, mat.T, ids

def load_seqs(d):
    out, name = {}, None
    for ln in open(os.path.join(d, "dada2/ASV_seqs.fasta")):
        if ln.startswith(">"):
            name = ln[1:].split()[0]; out[name] = ""
        else:
            out[name] += ln.strip()
    return out

def load_stats(d):
    out = {}
    with open(os.path.join(d, "dada2/DADA2_stats.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out[r[list(r)[0]]] = r
    return out

def val(r, *names):
    for n in names:
        if n in r and r[n] not in ("", "NA", None):
            try: return float(r[n])
            except ValueError: pass
    return None

def bray_curtis(X):
    s = X.sum(axis=1, keepdims=True)
    s[s == 0] = 1.0
    rel = X / s
    n = rel.shape[0]
    D = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        D[i, j] = D[j, i] = np.abs(rel[i] - rel[j]).sum() / 2.0
    return D

def mantel(A, B, nperm=NPERM):
    n = A.shape[0]
    iu = np.triu_indices(n, 1)
    a, b = A[iu], B[iu]
    r = np.corrcoef(a, b)[0, 1]
    ge = sum(1 for _ in range(nperm)
             if abs(np.corrcoef(A[np.ix_(p := RNG.permutation(n), p)][iu], b)[0, 1]) >= abs(r))
    return r, (ge + 1) / (nperm + 1)

meta = {r["ID"]: r for r in csv.DictReader(open("assets/metadata_pilot.tsv"), delimiter="\t")}

def merged_pct(stats, sid):
    r = stats.get(sid)
    if not r: return None
    inp = val(r, "DADA2_input", "input", "reads.in")
    mrg = val(r, "merged")
    return 100 * mrg / inp if inp and mrg is not None and inp > 0 else None

def core_frac(d):
    """Percent of reads in ASVs of on-target length (252-254 bp)."""
    seqs = load_seqs(d)
    rows = list(csv.reader(open(os.path.join(d, "dada2/ASV_table.tsv")), delimiter="\t"))
    tot = core = 0.0
    for r in rows[1:]:
        n = sum(float(x) for x in r[1:])
        tot += n
        if 252 <= len(seqs.get(r[0], "")) <= 254:
            core += n
    return 100 * core / tot if tot else 0.0

ref_s, ref_m, _ = load_table(REF[0])
if ref_s is None:
    raise SystemExit(f"reference table missing: {REF[0]}")
ref_stats = load_stats(REF[0])

settings = [REF] + [a for a in ALTS if load_table(a[0])[0] is not None]

print("=" * 96)
print(f"{'setting':16s} {'ASVs':>6s} {'median':>8s} {'min':>7s} "
      f"{'on-target':>11s} {'Mantel r':>10s} {'p':>8s}")
print(f"{'':16s} {'':>6s} {'merged%':>8s} {'merged%':>7s} {'reads %':>11s}")
print("-" * 96)

for d, lbl in settings:
    s, m, _ = load_table(d)
    st = load_stats(d)
    swabs = [x for x in s if meta.get(x, {}).get("SampleType") == "VaginalSwab"]
    v = sorted(p for p in (merged_pct(st, x) for x in swabs) if p is not None)
    cf = core_frac(d)
    if d == REF[0]:
        r_txt, p_txt = "—", "—"
    else:
        shared = [x for x in ref_s if x in set(s)]
        Dr = bray_curtis(ref_m[[ref_s.index(x) for x in shared]])
        Da = bray_curtis(m[[s.index(x) for x in shared]])
        r, p = mantel(Dr, Da)
        r_txt, p_txt = f"{r:.4f}", f"{p:.4f}"
    print(f"{lbl:16s} {m.shape[1]:>6d} {np.median(v):>8.1f} {v[0]:>7.1f} "
          f"{cf:>10.2f}% {r_txt:>10s} {p_txt:>8s}")

print("=" * 96)
print("\nmedian merged% by sequencing run (swabs only):\n")
runs = sorted({m_["Run"] for m_ in meta.values()})
print(f"  {'run':6s}" + "".join(f"{l:>14s}" for _, l in settings))
for rn in runs:
    cells = []
    for d, _ in settings:
        st = load_stats(d)
        v = [p for sid in st
             if (mm := meta.get(sid)) and mm["Run"] == rn
             and mm["SampleType"] == "VaginalSwab"
             and (p := merged_pct(st, sid)) is not None]
        cells.append(np.median(v) if v else float("nan"))
    best = max((c for c in cells if c == c), default=None)
    print(f"  {rn:6s}" + "".join(
        f"{c:>13.1f}{'*' if best and abs(c-best) < 1e-9 else ' '}" for c in cells))
print("\n  * = highest median for that run")
