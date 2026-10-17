#!/usr/bin/env bash -C -e -u -o pipefail
taxref_reformat_standard.sh \

#Giving out information
echo -e "--dada_ref_taxonomy: silva=138\n" >ref_taxonomy.silva_138.txt
echo -e "Title: Silva 138.1 prokaryotic SSU\n" >>ref_taxonomy.silva_138.txt
echo -e "Citation: Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Res. 2013 Jan;41(Database issue):D590-6. doi: 10.1093/nar/gks1219. Epub 2012 Nov 28. PMID: 23193283; PMCID: PMC3531112.\n" >>ref_taxonomy.silva_138.txt
echo "All entries: [title:Silva 138.1 prokaryotic SSU, file:[https://zenodo.org/records/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz, https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz], citation:Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Res. 2013 Jan;41(Database issue):D590-6. doi: 10.1093/nar/gks1219. Epub 2012 Nov 28. PMID: 23193283; PMCID: PMC3531112., fmtscript:taxref_reformat_standard.sh, dbversion:SILVA v138.1 (https://zenodo.org/records/4587955)]" >>ref_taxonomy.silva_138.txt

cat <<-END_VERSIONS > versions.yml
"NFCORE_AMPLISEQ:AMPLISEQ:FORMAT_TAXONOMY":
    sed: $(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
END_VERSIONS
