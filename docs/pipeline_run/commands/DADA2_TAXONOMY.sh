#!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    set.seed(100) # Initialize random number generator for reproducibility

    taxlevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

    seq <- getSequences("ASV_seqs.2.fasta", collapse = TRUE, silence = FALSE)
    taxa <- assignTaxonomy(seq, "assignTaxonomy.fna", taxLevels = taxlevels, minBoot = 50,tryRC = FALSE, multithread = 10, verbose=TRUE, outputBootstraps = TRUE)

    # (1) Make a data frame, add ASV_ID from seq
    tx <- data.frame(ASV_ID = names(seq), taxa, sequence = row.names(taxa$tax), row.names = names(seq))

    # (2) Set confidence to the bootstrap for the most specific taxon
    # extract columns with taxonomic values
    tax <- tx[ , grepl( "tax." , names( tx ) ) ]
    # find first occurrence of NA
    res <- max.col(is.na(tax), ties = "first")
    # correct if no NA is present in column to NA
    if(any(res == 1)) is.na(res) <- (res == 1) & !is.na(tax[[1]])
    # find index of last entry before NA, which is the bootstrap value
    res <- res-1
    # if NA choose last entry
    res[is.na(res)] <- ncol(tax)
    # extract bootstrap values
    boot <- tx[ , grepl( "boot." , names( tx ) ) ]
    boot$last_tax <- res
    valid_boot <- apply(boot,1,function(x) x[x[length(x)]][1]/100 )
    # replace missing bootstrap values (NA) with 0
    valid_boot[is.na(valid_boot)] <- 0
    # add bootstrap values to column confidence
    tx$confidence <- valid_boot

    # (3) Reorder columns before writing to file
    expected_order <- c("ASV_ID",paste0("tax.",taxlevels),"confidence","sequence")
    expected_order <- intersect(expected_order,colnames(tx))
    taxa_export <- subset(tx, select = expected_order)
    colnames(taxa_export) <- sub("tax.", "", colnames(taxa_export))
    rownames(taxa_export) <- names(seq)

    write.table(taxa_export, file = "ASV_seqs.2.ASV_tax.silva_138.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

    # Save a version with rownames for addSpecies
    taxa_export <- cbind( ASV_ID = tx$ASV_ID, taxa$tax, confidence = tx$confidence)
    saveRDS(taxa_export, "ASV_seqs.2.ASV_tax.silva_138.rds")

    write.table('assignTaxonomy	minBoot = 50,tryRC = FALSE
taxlevels	c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
seed	100', file = "assignTaxonomy.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    writeLines(c("\"NFCORE_AMPLISEQ:AMPLISEQ:DADA2_TAXONOMY_WF:DADA2_TAXONOMY\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
