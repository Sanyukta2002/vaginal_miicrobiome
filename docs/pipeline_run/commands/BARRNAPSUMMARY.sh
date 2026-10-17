#!/usr/bin/env bash -C -e -u -o pipefail
summarize_barrnap.py rrna.arc.gff rrna.bac.gff rrna.euk.gff rrna.mito.gff

if [[ $(wc -l < summary.tsv ) -le 1 ]]; then
    touch WARNING_no_rRNA_found_warning.txt
else
    touch no_warning.txt
fi

cat <<-END_VERSIONS > versions.yml
"NFCORE_AMPLISEQ:AMPLISEQ:BARRNAPSUMMARY":
    python: $( python --version )
END_VERSIONS
