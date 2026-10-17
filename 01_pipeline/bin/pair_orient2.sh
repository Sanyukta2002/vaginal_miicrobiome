#!/bin/bash
srr="$1"
r1=$(ls raw/fetchngs/fastq/*${srr}_1.fastq.gz 2>/dev/null | head -1)
r2=$(ls raw/fetchngs/fastq/*${srr}_2.fastq.gz 2>/dev/null | head -1)
orient () {
  zcat "$1" 2>/dev/null | head -80000 | awk '
    NR%4==2 { t++; if ($0 ~ /^CCT/) r++; else if ($0 ~ /^TAC/) f++ }
    END { if (t==0) { print "EMPTY\t0\t0"; exit }
          printf "%s\t%.1f\t%.1f\n", (r/t>0.8 ? "REV" : (f/t>0.8 ? "FWD" : "MIX")), 100*f/t, 100*r/t }'
}
read o1 f1 r1p <<< "$(orient "$r1")"
read o2 f2 r2p <<< "$(orient "$r2")"
if   [ "$o1" = "FWD" ] && [ "$o2" = "REV" ]; then cls="NORMAL"
elif [ "$o1" = "REV" ] && [ "$o2" = "FWD" ]; then cls="SWAPPED"
else cls="WEIRD"; fi
printf "%s\t%s(f%.0f/r%.0f)\t%s(f%.0f/r%.0f)\t%s\n" "$srr" "$o1" "$f1" "$r1p" "$o2" "$f2" "$r2p" "$cls"
