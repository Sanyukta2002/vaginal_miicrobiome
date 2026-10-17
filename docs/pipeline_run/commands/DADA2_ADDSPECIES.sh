#!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    set.seed(100) # Initialize random number generator for reproducibility

    #add "Species" if not already in taxlevels
    taxlevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    if ( !"Species" %in% taxlevels ) { taxlevels <- c(taxlevels,"Species") }

    taxtable <- readRDS("ASV_seqs.2.ASV_tax.silva_138.rds")

    #remove Species annotation from assignTaxonomy
    taxa_nospecies <- taxtable[,!colnames(taxtable) %in% 'Species']

    tx <- addSpecies(taxa_nospecies, "addSpecies.fna.gz", n = 1e5,,tryRC = FALSE, verbose=TRUE)

    # Create a table with specified column order
    tmp <- data.frame(row.names(tx)) # To separate ASV_ID from sequence
    expected_order <- c("ASV_ID",taxlevels,"confidence")
    taxa <- as.data.frame( subset(tx, select = expected_order) )
    taxa$sequence <- tmp[,1]
    row.names(taxa) <- row.names(tmp)

    #rename Species annotation to Species_exact
    colnames(taxa)[which(names(taxa) == "Species")] <- "Species_exact"

    #add Species annotation from assignTaxonomy again, after "Genus" column
    if ( "Species" %in% colnames(taxtable) ) {
        taxtable <- data.frame(taxtable)
        taxa_export <- data.frame(append(taxa, list(Species=taxtable$Species), after=match("Genus", names(taxa))))
    } else {
        taxa_export <- taxa
    }

    write.table(taxa_export, file = "ASV_seqs.2.ASV_tax.silva_138.species.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

    write.table('addSpecies	n = 1e5,,tryRC = FALSE
taxlevels	c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
seed	100', file = "addSpecies.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    writeLines(c("\"NFCORE_AMPLISEQ:AMPLISEQ:DADA2_TAXONOMY_WF:DADA2_ADDSPECIES\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
