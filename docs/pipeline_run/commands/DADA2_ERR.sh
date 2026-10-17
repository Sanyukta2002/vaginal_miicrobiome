#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
set.seed(100) # Initialize random number generator for reproducibility



fnFs <- list.files(".", pattern = "_1.filt.fastq.gz", full.names = TRUE)
fnRs <- sub("_1.filt.fastq.gz", "_2.filt.fastq.gz", fnFs)

sink(file = "IL02.err.log")
errF <- learnErrors(fnFs, nbases = 1e8, nreads = NULL, randomize = TRUE, MAX_CONSIST = 10, OMEGA_C = 0,qualityType = "Auto",errorEstimationFunction = loessErrfun, multithread = 6, verbose = TRUE)
saveRDS(errF, "IL02_1.err.rds")
errR <- learnErrors(fnRs, nbases = 1e8, nreads = NULL, randomize = TRUE, MAX_CONSIST = 10, OMEGA_C = 0,qualityType = "Auto",errorEstimationFunction = loessErrfun, multithread = 6, verbose = TRUE)
saveRDS(errR, "IL02_2.err.rds")
sink(file = NULL)

pdf("IL02_1.err.pdf")
plotErrors(errF, nominalQ = TRUE)
dev.off()
svg("IL02_1.err.svg")
plotErrors(errF, nominalQ = TRUE)
dev.off()

pdf("IL02_2.err.pdf")
plotErrors(errR, nominalQ = TRUE)
dev.off()
svg("IL02_2.err.svg")
plotErrors(errR, nominalQ = TRUE)
dev.off()

sink(file = "IL02_1.err.convergence.txt")
dada2:::checkConvergence(errF)
sink(file = NULL)

sink(file = "IL02_2.err.convergence.txt")
dada2:::checkConvergence(errR)
sink(file = NULL)

write.table('learnErrors	nbases = 1e8, nreads = NULL, randomize = TRUE, MAX_CONSIST = 10, OMEGA_C = 0,qualityType = "Auto",errorEstimationFunction = loessErrfun', file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
writeLines(c("\"NFCORE_AMPLISEQ:AMPLISEQ:DADA2_ERR\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
