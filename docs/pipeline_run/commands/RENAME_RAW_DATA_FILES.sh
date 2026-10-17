#!/usr/bin/env bash -C -e -u -o pipefail
[ -f "S2007001548_1.fastq.gz" ] || cp "SRX3103526_SRR5944962_2.fastq.gz" "S2007001548_1.fastq.gz"
[ -f "S2007001548_2.fastq.gz" ] || cp "SRX3103526_SRR5944962_1.fastq.gz" "S2007001548_2.fastq.gz"

cat <<-END_VERSIONS > versions.yml
"NFCORE_AMPLISEQ:AMPLISEQ:RENAME_RAW_DATA_FILES":
    sed: $(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
END_VERSIONS
