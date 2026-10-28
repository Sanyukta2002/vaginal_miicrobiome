# Build the analysis table: 2,179 swabs x the variants shared with the deposit, with taxonomy and metadata.
# In:  data/processed (from 01), data/pipeline, data/external
# Out: data/processed/{counts_swabs_x_variants, taxonomy_variants, sample_metadata}.rds,
#      data/processed/genus_crosswalk_v123_vs_138.tsv

library(here)

load(here("data/external/processed.rda"))
our_col  <- readRDS(here("data/processed/our_col_collapsed_swabs.rds"))
asv2core <- readRDS(here("data/processed/asv2core_map.rds"))
our_tax  <- read.delim(here("data/pipeline/ASV_tax_species.silva_138.tsv"), row.names = 1, check.names = FALSE)
cur      <- read.delim(here("data/external/LactoSpecies_ByHand.tsv"), stringsAsFactors = FALSE)

# Count table: samples x shared variants, in the deposit's sample and variant order
shared <- colnames(st)[colnames(st) %in% rownames(our_col)]
our_counts <- t(our_col[shared, ])
rownames(our_counts) <- sub("^S", "", rownames(our_counts))
our_counts <- our_counts[rownames(st), ]
stopifnot(identical(rownames(our_counts), rownames(st)), min(rowSums(our_counts)) > 0)

# SILVA 138 genus per variant: first non-blank genus among its member ASVs
mem <- asv2core[asv2core$core %in% shared, ]
mem$Genus138 <- our_tax[mem$ASV_ID, "Genus"]
genus138 <- tapply(mem$Genus138, mem$core, function(g) {
  g <- g[!is.na(g) & g != ""]
  if (length(g) == 0) NA else g[1]
})[shared]

# Genus crosswalk: SILVA v123 (deposit) vs SILVA 138 (ours), classified by kind of difference
norm <- function(g) {
  g[!is.na(g) & g == ""] <- NA
  sub(" [0-9]+$", "", gsub("/", "-", gsub("_", " ", g)))
}
a <- tax[shared, "Genus"]; a[!is.na(a) & a == ""] <- NA
b <- genus138;            b[!is.na(b) & b == ""] <- NA
kind <- ifelse(is.na(a) & is.na(b), "both unassigned",
        ifelse(is.na(a), "new name in 138 (v123 had none)",
        ifelse(is.na(b), "lost in 138 (138 has none)",
        ifelse(a == b, "identical",
        ifelse(norm(a) == norm(b), "formatting only", "real rename")))))

full <- data.frame(kind = kind, Genus_v123 = a, Genus_138 = b, reads = colSums(our_counts))
full <- full[full$kind != "identical", ]
full$Genus_v123[is.na(full$Genus_v123)] <- "(none)"
full$Genus_138[is.na(full$Genus_138)]   <- "(none)"
pairs_all <- aggregate(reads ~ kind + Genus_v123 + Genus_138, full, sum)
pairs_all$n_variants <- aggregate(reads ~ kind + Genus_v123 + Genus_138, full, length)$reads
pairs_all <- pairs_all[order(pairs_all$kind, -pairs_all$reads), ]
write.table(pairs_all, here("data/processed/genus_crosswalk_v123_vs_138.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Taxonomy table: one representative ASV per variant (prefer one with a genus name)
rep_idx <- tapply(seq_len(nrow(mem)), mem$core, function(i) {
  g <- mem$Genus138[i]
  j <- i[!is.na(g) & g != ""]
  if (length(j)) j[1] else i[1]
})
rep_id <- mem$ASV_ID[rep_idx[shared]]

ranks  <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
taxtab <- our_tax[rep_id, ranks]
rownames(taxtab) <- shared
taxtab[taxtab == ""] <- NA
taxtab$Genus_v123       <- tax[shared, "Genus"]
taxtab$Species_SILVA138 <- our_tax[rep_id, "Species"]
taxtab$Species_SILVA138[taxtab$Species_SILVA138 == ""] <- NA

# Curated Lactobacillus species (exact sequence match) and Gardnerella G1-G3
taxtab$Species <- cur$species[match(shared, cur$sequence)]
gg_seq <- colnames(st)[which(tax[, "Genus"] %in% "Gardnerella")]
taxtab$Gvariant <- NA
taxtab[gg_seq[1:3][gg_seq[1:3] %in% shared], "Gvariant"] <- c("G1", "G2", "G3")
stopifnot(identical(rownames(taxtab), colnames(our_counts)))

# Sample metadata from the deposit, matched by SampleID
rownames(df) <- as.character(df$SampleID)
meta <- df[rownames(our_counts), ]
stopifnot(identical(rownames(meta), rownames(our_counts)))

saveRDS(our_counts, here("data/processed/counts_swabs_x_variants.rds"))
saveRDS(taxtab,     here("data/processed/taxonomy_variants.rds"))
saveRDS(meta,       here("data/processed/sample_metadata.rds"))
