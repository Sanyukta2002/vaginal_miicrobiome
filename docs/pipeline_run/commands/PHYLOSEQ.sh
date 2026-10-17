#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(phyloseq))

otu_df  <- read.table("ASV_table.tsv", sep="\t", header=TRUE, row.names=1)
tax_df  <- read.table("ASV_tax_species.silva_138.tsv", sep="\t", header=TRUE, row.names=1)
otu_mat <- as.matrix(otu_df)
tax_mat <- as.matrix(tax_df)

OTU     <- otu_table(otu_mat, taxa_are_rows=TRUE)
TAX     <- tax_table(tax_mat)
phy_obj <- phyloseq(OTU, TAX)

if (file.exists("metadata_ampliseq.tsv")) {
    sam_df  <- read.table("metadata_ampliseq.tsv", sep="\t", header=TRUE, row.names=1)
    SAM     <- sample_data(sam_df)
    phy_obj <- merge_phyloseq(phy_obj, SAM)
}

if (file.exists("")) {
    TREE    <- read_tree("")
    phy_obj <- merge_phyloseq(phy_obj, TREE)
}

saveRDS(phy_obj, file = paste0("dada2", "_phyloseq.rds"))

# Version information
writeLines(c("\"NFCORE_AMPLISEQ:AMPLISEQ:ROBJECT_WORKFLOW:PHYLOSEQ\":",
    paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
    paste0("    phyloseq: ", packageVersion("phyloseq"))),
    "versions.yml"
)
