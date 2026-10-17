#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

errF <- readRDS("IL02_1.err.rds")
errR <- readRDS("IL02_2.err.rds")

filtFs <- list.files("./filtered/", pattern = "_1.filt.fastq.gz", full.names = TRUE)
filtRs <- sub("_1.filt.fastq.gz", "_2.filt.fastq.gz", filtFs)

#denoising
sink(file = "IL02.dada.log")
if ("Auto" == "Auto") {
    # Avoid using memory-inefficient derepFastq() if not necessary
    dadaFs <- dada(filtFs, err = errF, selfConsist = FALSE, priors = character(0), DETECT_SINGLETONS = FALSE, GAPLESS = TRUE, GAP_PENALTY = -8, GREEDY = TRUE, KDIST_CUTOFF = 0.42, MATCH = 5, MAX_CLUST = 0, MAX_CONSIST = 10, MIN_ABUNDANCE = 1, MIN_FOLD = 1, MIN_HAMMING = 1, MISMATCH = -4, OMEGA_A = 1e-40, OMEGA_C = 1e-40, OMEGA_P = 1e-4, PSEUDO_ABUNDANCE = Inf, PSEUDO_PREVALENCE = 2, SSE = 2, USE_KMERS = TRUE, USE_QUALS = TRUE, VECTORIZED_ALIGNMENT = TRUE,BAND_SIZE = 16, HOMOPOLYMER_GAP_PENALTY = NULL,pool = FALSE, multithread = 10)
    dadaRs <- dada(filtRs, err = errR, selfConsist = FALSE, priors = character(0), DETECT_SINGLETONS = FALSE, GAPLESS = TRUE, GAP_PENALTY = -8, GREEDY = TRUE, KDIST_CUTOFF = 0.42, MATCH = 5, MAX_CLUST = 0, MAX_CONSIST = 10, MIN_ABUNDANCE = 1, MIN_FOLD = 1, MIN_HAMMING = 1, MISMATCH = -4, OMEGA_A = 1e-40, OMEGA_C = 1e-40, OMEGA_P = 1e-4, PSEUDO_ABUNDANCE = Inf, PSEUDO_PREVALENCE = 2, SSE = 2, USE_KMERS = TRUE, USE_QUALS = TRUE, VECTORIZED_ALIGNMENT = TRUE,BAND_SIZE = 16, HOMOPOLYMER_GAP_PENALTY = NULL,pool = FALSE, multithread = 10)
} else {
    derepFs <- derepFastq(filtFs, qualityType="Auto")
    dadaFs <- dada(derepFs, err = errF, selfConsist = FALSE, priors = character(0), DETECT_SINGLETONS = FALSE, GAPLESS = TRUE, GAP_PENALTY = -8, GREEDY = TRUE, KDIST_CUTOFF = 0.42, MATCH = 5, MAX_CLUST = 0, MAX_CONSIST = 10, MIN_ABUNDANCE = 1, MIN_FOLD = 1, MIN_HAMMING = 1, MISMATCH = -4, OMEGA_A = 1e-40, OMEGA_C = 1e-40, OMEGA_P = 1e-4, PSEUDO_ABUNDANCE = Inf, PSEUDO_PREVALENCE = 2, SSE = 2, USE_KMERS = TRUE, USE_QUALS = TRUE, VECTORIZED_ALIGNMENT = TRUE,BAND_SIZE = 16, HOMOPOLYMER_GAP_PENALTY = NULL,pool = FALSE, multithread = 10)
    derepRs <- derepFastq(filtRs, qualityType="Auto")
    dadaRs <- dada(derepRs, err = errR, selfConsist = FALSE, priors = character(0), DETECT_SINGLETONS = FALSE, GAPLESS = TRUE, GAP_PENALTY = -8, GREEDY = TRUE, KDIST_CUTOFF = 0.42, MATCH = 5, MAX_CLUST = 0, MAX_CONSIST = 10, MIN_ABUNDANCE = 1, MIN_FOLD = 1, MIN_HAMMING = 1, MISMATCH = -4, OMEGA_A = 1e-40, OMEGA_C = 1e-40, OMEGA_P = 1e-4, PSEUDO_ABUNDANCE = Inf, PSEUDO_PREVALENCE = 2, SSE = 2, USE_KMERS = TRUE, USE_QUALS = TRUE, VECTORIZED_ALIGNMENT = TRUE,BAND_SIZE = 16, HOMOPOLYMER_GAP_PENALTY = NULL,pool = FALSE, multithread = 10)
}
saveRDS(dadaFs, "IL02_1.dada.rds")
saveRDS(dadaRs, "IL02_2.dada.rds")
sink(file = NULL)

# merge
if ("merge" == "consensus") {
    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, homo_gap = NULL, endsfree = TRUE, vec = FALSE, propagateCol = character(0), trimOverhang = FALSE, returnRejects = TRUE,justConcatenate = FALSE, match = 1, mismatch = -64, gap = -64, minOverlap = 12, maxMismatch = 0, justConcatenate = FALSE, verbose=TRUE)
    concats <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, homo_gap = NULL, endsfree = TRUE, vec = FALSE, propagateCol = character(0), trimOverhang = FALSE, returnRejects = TRUE,justConcatenate = FALSE, match = 1, mismatch = -64, gap = -64, minOverlap = 12, maxMismatch = 0, justConcatenate = TRUE, verbose=TRUE)

    # in case there is only one sample in the entire run
    if (is.data.frame(mergers)) {
        mergers <- list(sample = mergers)
        concats <- list(sample = concats)
    }

    # define the overlap threshold to decide if concatenation or not
    min_overlap_obs <- lapply(mergers, function(X) {
        mergers_accepted <- X[["accept"]]
        if (sum(mergers_accepted) > 0) {
            min_overlap_obs <- X[["nmatch"]][mergers_accepted] + X[["nmismatch"]][mergers_accepted]
            rep(min_overlap_obs, X[["abundance"]][mergers_accepted])
        } else {
            NA
        }
    })

    min_overlap_obs <- Reduce(c, min_overlap_obs)
    min_overlap_obs <- min_overlap_obs[!is.na(min_overlap_obs)]
    min_overlap_obs <- quantile(min_overlap_obs, 0.001)

    for (x in names(mergers)) {
        to_concat <- !mergers[[x]][["accept"]] & (mergers[[x]][["nmismatch"]] + mergers[[x]][["nmatch"]]) < min_overlap_obs

        if (sum(to_concat) > 0) {
            mergers[[x]][to_concat, ] <- concats[[x]][to_concat, ]
            # filter out unaccepted non concatenated sequences
            mergers[[x]] <- mergers[[x]][mergers[[x]][["accept"]], ]
        }

    }

    # if one sample, need to convert back to df for next steps

    if(length(mergers) == 1) {
        mergers <- mergers[[1]]
    }

} else {
    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, homo_gap = NULL, endsfree = TRUE, vec = FALSE, propagateCol = character(0), trimOverhang = FALSE, returnRejects = TRUE,justConcatenate = FALSE, match = 1, mismatch = -64, gap = -64, minOverlap = 12, maxMismatch = 0, verbose=TRUE)
}

saveRDS(mergers, "IL02.mergers.rds")

# make table
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "IL02.seqtab.rds")

write.table('dada	selfConsist = FALSE, priors = character(0), DETECT_SINGLETONS = FALSE, GAPLESS = TRUE, GAP_PENALTY = -8, GREEDY = TRUE, KDIST_CUTOFF = 0.42, MATCH = 5, MAX_CLUST = 0, MAX_CONSIST = 10, MIN_ABUNDANCE = 1, MIN_FOLD = 1, MIN_HAMMING = 1, MISMATCH = -4, OMEGA_A = 1e-40, OMEGA_C = 1e-40, OMEGA_P = 1e-4, PSEUDO_ABUNDANCE = Inf, PSEUDO_PREVALENCE = 2, SSE = 2, USE_KMERS = TRUE, USE_QUALS = TRUE, VECTORIZED_ALIGNMENT = TRUE,BAND_SIZE = 16, HOMOPOLYMER_GAP_PENALTY = NULL,pool = FALSE', file = "dada.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
write.table('mergePairs	homo_gap = NULL, endsfree = TRUE, vec = FALSE, propagateCol = character(0), trimOverhang = FALSE, returnRejects = TRUE,justConcatenate = FALSE, match = 1, mismatch = -64, gap = -64, minOverlap = 12, maxMismatch = 0', file = "mergePairs.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
writeLines(c("\"NFCORE_AMPLISEQ:AMPLISEQ:DADA2_DENOISING\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
