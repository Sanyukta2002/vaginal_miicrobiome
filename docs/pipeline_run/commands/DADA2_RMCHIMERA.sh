#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

seqtab = readRDS("IL02.seqtab.rds")

#remove chimera
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", minSampleFraction = 0.9, ignoreNNegatives = 1, minFoldParentOverAbundance = 2, minParentAbundance = 8, allowOneOff = FALSE, minOneOffParentDistance = 4, maxShift = 16, multithread=6, verbose=TRUE)
if ( 182 == 1 ) { rownames(seqtab.nochim) <- "S1000003920" }
saveRDS(seqtab.nochim,"IL02.ASVtable.rds")

write.table('removeBimeraDenovo	method="consensus", minSampleFraction = 0.9, ignoreNNegatives = 1, minFoldParentOverAbundance = 2, minParentAbundance = 8, allowOneOff = FALSE, minOneOffParentDistance = 4, maxShift = 16', file = "removeBimeraDenovo.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
writeLines(c("\"NFCORE_AMPLISEQ:AMPLISEQ:DADA2_RMCHIMERA\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
