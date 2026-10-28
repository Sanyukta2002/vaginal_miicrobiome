#!/bin/bash
# Re-download runs that fetchngs failed on (malformed ENA FTP paths): https, trailing slash removed, md5 checked.
# Run from 01_pipeline/.
set -u
TABLE=assets/SRP115697_ena_runs.tsv
OUTDIR=raw/fetchngs/fastq

awk -F'\t' 'NR>1{print $1}' "$TABLE" | sort > /tmp/all_runs.txt
ls "$OUTDIR" | grep -o 'SRR[0-9]*' | sort -u > /tmp/have_runs.txt
comm -23 /tmp/all_runs.txt /tmp/have_runs.txt > assets/missing_runs.txt
echo "missing runs: $(wc -l < assets/missing_runs.txt)"

while read -r RUN; do
  [ -z "$RUN" ] && continue
  FTP=$(awk -F'\t' -v r="$RUN" 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next} $h["run_accession"]==r {print $h["fastq_ftp"]}' "$TABLE")
  MD5=$(awk -F'\t' -v r="$RUN" 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next} $h["run_accession"]==r {print $h["fastq_md5"]}' "$TABLE")
  i=1
  for U in $(echo "$FTP" | tr ';' ' '); do
    U=${U%/}
    OUT="$OUTDIR/${RUN}_${i}.fastq.gz"
    wget -q -O "$OUT" "https://${U}"
    EXPECT=$(echo "$MD5" | tr ';' ' ' | awk -v n=$i '{print $n}')
    GOT=$(md5sum "$OUT" | cut -d' ' -f1)
    if [ "$EXPECT" = "$GOT" ]; then echo "OK   $OUT"; else echo "FAIL $OUT"; fi
    i=$((i+1))
  done
done < assets/missing_runs.txt
