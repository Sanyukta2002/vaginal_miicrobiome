#!/usr/bin/env bash -C -e -u -o pipefail
IFS=',' read -r -a kingdom <<< "bac,arc,mito,euk"

for KINGDOM in "${kingdom[@]}"
do
    barrnap \
        --threads 2 \
        --quiet --reject 0.1 \
        --kingdom $KINGDOM \
        --outseq Filtered.${KINGDOM}.fasta \
        < ASV_seqs.fasta \
        > rrna.${KINGDOM}.gff

    #this fails when processing an empty file, so it requires a workaround!
    if [ -s Filtered.${KINGDOM}.fasta ]; then
        grep -h '>' Filtered.${KINGDOM}.fasta | sed 's/^>//' | sed 's/:\+/\t/g' | awk '{print $2}' | sort -u >${KINGDOM}.matches.txt
    else
        touch ${KINGDOM}.matches.txt
    fi
done

cat <<-END_VERSIONS > versions.yml
"NFCORE_AMPLISEQ:AMPLISEQ:BARRNAP":
    barrnap: $(echo $(barrnap --version 2>&1) | sed "s/^.*barrnap //g")
END_VERSIONS
