#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(TreeSummarizedExperiment))

# Read otu table. It must be in a SimpleList as a matrix where rows
# represent taxa and columns samples.
otu_mat  <- read.table("ASV_table.tsv", sep="\t", header=TRUE, row.names=1)
otu_mat <- as.matrix(otu_mat)
assays <- SimpleList(counts = otu_mat)
# Read taxonomy table. Correct format for it is DataFrame.
taxonomy_table  <- read.table("ASV_tax_species.silva_138.tsv", sep="\t", header=TRUE, row.names=1)
taxonomy_table <- DataFrame(taxonomy_table)

# Match rownames between taxonomy table and abundance matrix.
taxonomy_table <- taxonomy_table[match(rownames(otu_mat), rownames(taxonomy_table)), ]

# Create TreeSE object.
tse <- TreeSummarizedExperiment(
    assays = assays,
    rowData = taxonomy_table
)

# If taxonomy table contains sequences, move them to referenceSeq slot
if (!is.null(rowData(tse)[["sequence"]])) {
    referenceSeq(tse) <- DNAStringSet( rowData(tse)[["sequence"]] )
    rowData(tse)[["sequence"]] <- NULL
}

# If provided, we add sample metadata as DataFrame object. rownames of
# sample metadata must match with colnames of abundance matrix.
if (file.exists("metadata_ampliseq.tsv")) {
    sample_meta  <- read.table("metadata_ampliseq.tsv", sep="\t", header=TRUE, row.names=1)
    sample_meta <- sample_meta[match(colnames(tse), rownames(sample_meta)), ]
    sample_meta  <- DataFrame(sample_meta)
    colData(tse) <- sample_meta
}

# If provided, we add phylogeny. The rownames in abundance matrix must match
# with node labels in phylogeny.
if (file.exists("")) {
    phylogeny <- ape::read.tree("")
    rowTree(tse) <- phylogeny
}

saveRDS(tse, file = paste0("dada2", "_TreeSummarizedExperiment.rds"))

# Version information
writeLines(c("\"NFCORE_AMPLISEQ:AMPLISEQ:ROBJECT_WORKFLOW:TREESUMMARIZEDEXPERIMENT\":",
    paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
    paste0("    TreeSummarizedExperiment: ", packageVersion("TreeSummarizedExperiment"))),
    "versions.yml"
)
