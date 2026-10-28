# Collapse our ASVs to the original study's 235-bp window and match them to its deposited variants.
# In:  data/pipeline (ampliseq), data/external (study deposit)
# Out: data/processed/{our_col_collapsed_swabs, asv2core_map}.rds, results/01_match_variants/

library(here)
library(phyloseq)

out <- here("03_analysis/results/01_match_variants")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

load(here("data/external/processed.rda"))
ps_all   <- readRDS(here("data/external/ps_allsams.rds"))
our_tab  <- read.delim(gzfile(here("data/pipeline/ASV_table.tsv.gz")), row.names = 1, check.names = FALSE)
our_tax  <- read.delim(here("data/pipeline/ASV_tax.silva_138.tsv"), row.names = 1, check.names = FALSE)
our_meta <- read.delim(here("data/pipeline/metadata_ampliseq.tsv"), check.names = FALSE)

# Vaginal swabs only, ASVs with reads in them
swab_ids <- our_meta$ID[our_meta$SampleType == "VaginalSwab"]
our_swab <- our_tab[, swab_ids]
stopifnot(setequal(sub("^S", "", colnames(our_swab)), rownames(st)))
our_swab <- our_swab[rowSums(our_swab) > 0, ]

# Trim to the 235-bp window (bases 11 to 245, or length - 10 for long ASVs) and merge identical cores
our_seq  <- our_tax[rownames(our_swab), "sequence"]
our_core <- substr(our_seq, 11, pmax(245, nchar(our_seq) - 10))
asv2core <- data.frame(ASV_ID = rownames(our_swab), core = our_core)
our_col  <- rowsum(as.matrix(our_swab), group = our_core)
stopifnot(sum(our_col) == sum(our_swab))

# Variant overlap with the deposited table
their_seqs  <- colnames(st)
our_seqs    <- rownames(our_col)
their_reads <- colSums(st)
our_reads   <- rowSums(our_col)
in_ours     <- their_seqs %in% our_seqs
in_theirs   <- our_seqs %in% their_seqs

# Study filter (>= 3 samples and >= 100 reads), counted over swabs and over all 6,088 deposited samples
keep_rule <- function(m) rowSums(m > 0) >= 3 & rowSums(m) >= 100
otu <- as(otu_table(ps_all), "matrix")
if (taxa_are_rows(ps_all)) otu <- t(otu)
pass_all <- colSums(otu > 0) >= 3 & colSums(otu) >= 100

summary_tab <- data.frame(
  metric = c("ASVs, all samples", "ASVs with reads in swabs", "variants after collapsing to 235 bp",
             "deposited variants", "deposited variants found in ours", "% deposited reads found in ours",
             "our variants found in deposit", "% our reads found in deposit",
             "our variants passing filter (swabs)", "deposited variants passing filter (swabs)",
             "deposited variants passing filter (all 6,088 samples)"),
  value = c(nrow(our_tab), nrow(our_swab), nrow(our_col),
            length(their_seqs), sum(in_ours), round(100 * sum(their_reads[in_ours]) / sum(their_reads), 3),
            sum(in_theirs), round(100 * sum(our_reads[in_theirs]) / sum(our_reads), 3),
            sum(keep_rule(our_col)), sum(keep_rule(t(st))), sum(pass_all[their_seqs], na.rm = TRUE)))
write.csv(summary_tab, file.path(out, "variant_matching_summary.csv"), row.names = FALSE)

# Gardnerella variants: deposited (abundance order) vs ours, by exact containment
gg <- which(tax[, "Genus"] %in% "Gardnerella")
our_g <- rownames(our_tax)[our_tax$Genus %in% "Gardnerella"]
our_gard <- data.frame(ASV_ID = our_g, our_reads = rowSums(our_swab[our_g, ]), our_seq = our_tax[our_g, "sequence"])
our_gard <- our_gard[order(-our_gard$our_reads), ]
our_gard$our_rank <- seq_len(nrow(our_gard))

gard_match <- do.call(rbind, lapply(seq_along(gg), function(i) {
  s   <- their_seqs[gg[i]]
  hit <- which(grepl(s, our_gard$our_seq, fixed = TRUE))
  data.frame(label = c("G1", "G2", "G3", rep("", length(gg) - 3))[i], their_rank = i,
             their_reads = unname(their_reads[gg[i]]), our_rank = our_gard$our_rank[hit],
             ASV_ID = our_gard$ASV_ID[hit], our_reads = our_gard$our_reads[hit],
             offset = regexpr(s, our_gard$our_seq[hit], fixed = TRUE) - 1)
}))
write.csv(gard_match, file.path(out, "gardnerella_variant_match.csv"), row.names = FALSE)

saveRDS(our_col,  here("data/processed/our_col_collapsed_swabs.rds"))
saveRDS(asv2core, here("data/processed/asv2core_map.rds"))
