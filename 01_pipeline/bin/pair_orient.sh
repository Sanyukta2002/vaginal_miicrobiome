#!/bin/bash
srr="$1"
r1=$(ls raw/fetchngs/fastq/*${srr}_1.fastq.gz 2>/dev/null | head -1)
r2=$(ls raw/fetchngs/fastq/*${srr}_2.fastq.gz 2>/dev/null | head -1)
orient () {
  zcat "$1" 2>/dev/null | head -8000 | awk '
    NR%4==2 { t++; if ($0 ~ /^.CTGTT/) r++; else if ($0 ~ /^[ACGTN][AGT][CT][ACGT]/) f++ }
    END { if (t==0) { print "EMPTY"; exit }
          if (r/t > 0.8) print "REV"; else if (f/t > 0.7) print "FWD"; else print "MIX" }'
}
o1=$(orient "$r1"); o2=$(orient "$r2")
if   [ "$o1" = "FWD" ] && [ "$o2" = "REV" ]; then cls="NORMAL"
elif [ "$o1" = "REV" ] && [ "$o2" = "FWD" ]; then cls="SWAPPED"
else cls="WEIRD"; fi
printf "%s\t%s\t%s\t%s\n" "$srr" "$o1" "$o2" "$cls"
