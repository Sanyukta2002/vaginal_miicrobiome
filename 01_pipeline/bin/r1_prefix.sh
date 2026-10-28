#!/bin/bash
# Print the first 12 bases of the first read in a FASTQ file.
f="$1"
acc=$(basename "$f" | grep -oE 'SRR[0-9]+')
p=$(zcat "$f" 2>/dev/null | head -2 | tail -1 | cut -c1-12)
printf "%s\t%s\n" "$acc" "$p"
