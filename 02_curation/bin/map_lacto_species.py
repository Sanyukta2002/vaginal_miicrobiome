#!/usr/bin/env python3
"""Apply the authors' hand curation of Lactobacillus V4 species to our ASVs.

Why: SILVA 138's species-level assignment cannot separate the L. acidophilus
complex on V4. Our run calls 19.98% of reads "Lactobacillus acidophilus" and
0.00% "L. crispatus", which is backwards for vaginal communities -- the signal
the paper is built on. Callahan et al. solved this by curating the V4 sequences
of the Lactobacillus species that actually occur here (LactoSpecies_ByHand.rds,
40 sequences / 14 species) and overwriting the Species column with it.

Their curated sequences are 235 bp; ours are 253 bp, because their trimming
started ~10 bases further in and ended ~8 earlier. So the match is exact
containment of their sequence inside ours, not equality.

This script reports before deciding anything: offsets, coverage, ambiguity, and
the resulting shift in species abundance. It writes the ASV -> species map but
does NOT modify the ASV table or taxonomy file; that is a separate step once
these numbers look right.
"""
import os
import sys
from collections import Counter, defaultdict

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(PROJ, "results/main/dada2")
FASTA = os.path.join(D, "ASV_seqs.fasta")
TAX = os.path.join(D, "ASV_tax_species.silva_138.tsv")
TABLE = os.path.join(D, "ASV_table.tsv")
CURATED = os.path.join(PROJ, "assets/LactoSpecies_ByHand.tsv")
OUT = os.path.join(PROJ, "assets/lacto_species_map.tsv")

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def read_tax(path):
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        for col in ("ASV_ID", "Genus", "Species"):
            if col not in idx:
                sys.exit("taxonomy file lacks column %r; has %s" % (col, hdr))
        rows = {}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(hdr):
                p += [""] * (len(hdr) - len(p))
            rows[p[idx["ASV_ID"]]] = (p[idx["Genus"]], p[idx["Species"]])
    return rows


def read_totals(path):
    """Total reads per ASV, summed across all samples."""
    totals = {}
    with open(path) as f:
        f.readline()                       # header: ASV_ID + sample columns
        for line in f:
            p = line.rstrip("\n").split("\t")
            totals[p[0]] = sum(int(x) for x in p[1:] if x)
    return totals


def read_curated(path):
    out = []
    with open(path) as f:
        f.readline()
        for line in f:
            sp, _len, seq = line.rstrip("\n").split("\t")
            out.append((sp, seq))
    return out


def pct(n, d):
    return 0.0 if not d else 100.0 * n / d


def main():
    asv_seq = read_fasta(FASTA)
    tax = read_tax(TAX)
    totals = read_totals(TABLE)
    curated = read_curated(CURATED)

    print("inputs")
    print("  ASVs in fasta        : %d" % len(asv_seq))
    print("  ASVs in taxonomy     : %d" % len(tax))
    print("  ASVs in count table  : %d" % len(totals))
    print("  curated sequences    : %d (%d species)"
          % (len(curated), len({s for s, _ in curated})))
    grand = sum(totals.values())
    print("  total reads          : %d" % grand)

    missing = set(asv_seq) - set(tax)
    if missing:
        print("  WARNING %d ASVs in fasta absent from taxonomy" % len(missing))

    # ---- match ---------------------------------------------------------
    hits = defaultdict(list)               # asv_id -> [(species, offset)]
    rc_only = []
    offsets = Counter()
    for asv, seq in asv_seq.items():
        for sp, cseq in curated:
            off = seq.find(cseq)
            if off >= 0:
                hits[asv].append((sp, off))
                offsets[off] += 1
            elif seq.find(rc(cseq)) >= 0:
                rc_only.append((asv, sp))

    matched_curated = {sp for v in hits.values() for sp, _ in v}
    unmatched_curated = [(sp, s) for sp, s in curated
                         if not any(s in q for q in asv_seq.values())]

    print("\nmatching")
    print("  ASVs with >=1 curated hit : %d" % len(hits))
    print("  offset distribution       : %s"
          % ", ".join("%d bp x%d" % (o, n) for o, n in sorted(offsets.items())))
    print("  curated seqs matched      : %d / %d"
          % (len(curated) - len(unmatched_curated), len(curated)))
    if unmatched_curated:
        print("  unmatched curated sequences:")
        for sp, s in unmatched_curated:
            print("    %-20s %s..." % (sp, s[:40]))
    if rc_only:
        print("  WARNING %d reverse-complement-only matches (orientation issue):"
              % len(rc_only))
        for asv, sp in rc_only[:5]:
            print("    %s  %s" % (asv, sp))

    ambiguous = {a: v for a, v in hits.items() if len({sp for sp, _ in v}) > 1}
    if ambiguous:
        print("  WARNING %d ASVs match more than one species:" % len(ambiguous))
        for a, v in sorted(ambiguous.items(),
                           key=lambda kv: -totals.get(kv[0], 0))[:10]:
            print("    %s  %8d reads  %s"
                  % (a, totals.get(a, 0), ", ".join(sorted({sp for sp, _ in v}))))

    # ---- sanity: did we match things SILVA did not call Lactobacillus? --
    not_lacto = [(a, tax.get(a, ("", ""))[0], totals.get(a, 0))
                 for a in hits if tax.get(a, ("", ""))[0] != "Lactobacillus"]
    print("\nconsistency with SILVA genus")
    print("  matched ASVs SILVA did NOT call Lactobacillus : %d" % len(not_lacto))
    for a, g, n in sorted(not_lacto, key=lambda x: -x[2])[:10]:
        print("    %s  %8d reads  SILVA genus=%r" % (a, n, g or "unassigned"))

    lacto_asvs = [a for a, (g, _) in tax.items() if g == "Lactobacillus"]
    lacto_reads = sum(totals.get(a, 0) for a in lacto_asvs)
    covered = sum(totals.get(a, 0) for a in hits
                  if tax.get(a, ("", ""))[0] == "Lactobacillus")
    print("  genus Lactobacillus : %d ASVs, %d reads (%.2f%% of all)"
          % (len(lacto_asvs), lacto_reads, pct(lacto_reads, grand)))
    print("  of those, curated   : %d ASVs, %d reads (%.2f%% of Lactobacillus)"
          % (len([a for a in hits if tax.get(a, ("", ""))[0] == "Lactobacillus"]),
             covered, pct(covered, lacto_reads)))

    # ---- before / after -------------------------------------------------
    before = Counter()
    after = Counter()
    for a, n in totals.items():
        g, s = tax.get(a, ("", ""))
        b = ("%s %s" % (g, s)).strip() if s else (g or "unassigned")
        before[b] += n
        if a in hits:
            sp = sorted({x for x, _ in hits[a]})
            after["Lactobacillus " + ("/".join(sp))] += n
        else:
            after[b] += n

    print("\nLactobacillus species abundance, SILVA vs curated")
    print("  %-42s %12s %7s   %12s %7s" % ("label", "before", "%", "after", "%"))
    labels = sorted({k for k in list(before) + list(after)
                     if k.startswith("Lactobacillus")},
                    key=lambda k: -max(before.get(k, 0), after.get(k, 0)))
    for k in labels:
        b, a_ = before.get(k, 0), after.get(k, 0)
        if b == 0 and a_ == 0:
            continue
        print("  %-42s %12d %6.2f%%   %12d %6.2f%%"
              % (k[:42], b, pct(b, grand), a_, pct(a_, grand)))

    # ---- write map ------------------------------------------------------
    with open(OUT, "w") as out:
        out.write("ASV_ID\tcurated_species\toffset\treads\tsilva_genus\tsilva_species\n")
        for a in sorted(hits, key=lambda x: -totals.get(x, 0)):
            sp = "/".join(sorted({x for x, _ in hits[a]}))
            off = hits[a][0][1]
            g, s = tax.get(a, ("", ""))
            out.write("%s\t%s\t%d\t%d\t%s\t%s\n"
                      % (a, sp, off, totals.get(a, 0), g, s))
    print("\nwrote %s (%d ASVs)" % (OUT, len(hits)))


if __name__ == "__main__":
    main()
