#!/usr/bin/env Rscript
library(rmarkdown)

# Work around  https://github.com/rstudio/rmarkdown/issues/1508
# If the symbolic link is not replaced by a physical file
# output- and temporary files will be written to the original directory.
file.copy("./report_template.Rmd", "./template.Rmd", overwrite = TRUE)

rmarkdown::render("template.Rmd", output_file = "summary_report.html", params = list(css='nf-core_style.css',report_logo='nf-core-ampliseq_logo_light_long.png',workflow_manifest_version='2.18.0',workflow_scriptid='be54f69b2e',report_title='Summary of analysis results',metadata='metadata_ampliseq.tsv',input_samplesheet='samplesheet_all.csv',mqc_plot='multiqc_plots/svg/fastqc_per_sequence_quality_scores_plot.svg',truncq=2,trunclenf='245',trunclenr='235',max_ee=2,dada_filtntrim_args='filterAndTrim.args.txt',dada_sample_inference='independent',dada_err_path='IL02_1.err.svg,IL02_2.err.svg,IL03_1.err.svg,IL03_2.err.svg,IL05_1.err.svg,IL05_2.err.svg,IL07_1.err.svg,IL07_2.err.svg,IL08_1.err.svg,IL08_2.err.svg,IL06_1.err.svg,IL06_2.err.svg,IL09_1.err.svg,IL09_2.err.svg',dada_err_run='IL02,IL03,IL05,IL07,IL08,IL06,IL09',asv_table_path='ASV_table.tsv',path_asv_fa='ASV_seqs.fasta',path_dada2_tab='DADA2_table.tsv',dada_stats_path='DADA2_stats.tsv',mergepairs_strategy='merge',path_barrnap_sum='summary.tsv',min_len_asv=0,max_len_asv=0,dada_min_boot=50,dada2_taxonomy='ASV_tax_species.silva_138.tsv',dada2_ref_tax_title='Silva 138.1 prokaryotic SSU',dada2_ref_tax_file='[https://zenodo.org/records/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz, https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz]',dada2_ref_tax_citation='Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Res. 2013 Jan;41(Database issue):D590-6. doi: 10.1093/nar/gks1219. Epub 2012 Nov 28. PMID: 23193283; PMCID: PMC3531112.',phyloseq='phyloseq,dada2_phyloseq.rds',tse='tse,dada2_TreeSummarizedExperiment.rds'), envir = new.env())

writeLines(c("\"NFCORE_AMPLISEQ:AMPLISEQ:SUMMARY_REPORT\":",
    paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
    paste0("    rmarkdown: ", packageVersion("rmarkdown")),
    paste0("    knitr: ", packageVersion("knitr")) ),
    "versions.yml")
