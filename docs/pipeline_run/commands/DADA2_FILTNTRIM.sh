#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

out <- filterAndTrim("S2007001548_1.fastq.gz", "S2007001548_1.filt.fastq.gz", "S2007001548_2.fastq.gz", "S2007001548_2.filt.fastq.gz",
    truncLen = c(245, 235),
    maxN = 0, trimRight = 0, minQ = 0, rm.lowcomplex = 0, orient.fwd = NULL, matchIDs = FALSE, id.sep = "\\s", id.field = NULL, n = 1e+05, OMP = TRUE,qualityType = "Auto",truncQ = 2,minLen = 50,maxLen = Inf,trimLeft = 0,maxEE = c(2, 2),rm.phix = TRUE,
    compress = TRUE,
    multithread = 2,
    verbose = TRUE)
out <- cbind(out, ID = row.names(out))

# If no reads passed the filter, write an empty GZ file
if(out[2] == '0'){
    for(fp in c("S2007001548_1.filt.fastq.gz", "S2007001548_2.filt.fastq.gz")){
        print(paste("Writing out an empty file:", fp))
        handle <- gzfile(fp, "w")
        write("", handle)
        close(handle)
    }
}

write.table( out, file = "S2007001548.filter_stats.tsv", sep = "\t", row.names = FALSE, quote = FALSE, na = '')
write.table(paste('filterAndTrim	truncLen = c(245, 235)','maxN = 0, trimRight = 0, minQ = 0, rm.lowcomplex = 0, orient.fwd = NULL, matchIDs = FALSE, id.sep = "\\s", id.field = NULL, n = 1e+05, OMP = TRUE,qualityType = "Auto",truncQ = 2,minLen = 50,maxLen = Inf,trimLeft = 0,maxEE = c(2, 2),rm.phix = TRUE',sep=","), file = "filterAndTrim.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
writeLines(c("\"NFCORE_AMPLISEQ:AMPLISEQ:DADA2_PREPROCESSING:DADA2_FILTNTRIM\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
