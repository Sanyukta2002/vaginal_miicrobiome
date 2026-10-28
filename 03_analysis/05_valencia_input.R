# Map each variant to a VALENCIA taxon name and write the VALENCIA input table.
# In:  data/processed (from 02), data/external, tools/VALENCIA
# Out: data/processed/valencia/{valencia_input, valencia_taxon_map}.csv

library(here)

d11 <- here("data/processed/valencia")
dir.create(d11, recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(here("data/processed/counts_swabs_x_variants.rds"))
taxtab <- readRDS(here("data/processed/taxonomy_variants.rds"))
meta   <- readRDS(here("data/processed/sample_metadata.rds"))
refs   <- read.delim(here("data/external/reference_16S_BVAB_Mgirerdii.tsv"), stringsAsFactors = FALSE)
cent   <- names(read.csv(here("tools/VALENCIA/CST_centroids_012920.csv"), check.names = FALSE, nrows = 1))[-1]

cl <- function(x) gsub("[^A-Za-z0-9_]", ".", x)
ok <- function(x) !is.na(x) & x != ""

# Genus label (deposit genus, else SILVA 138 genus), else the deepest available rank
gen <- ifelse(ok(taxtab$Genus_v123), taxtab$Genus_v123, ifelse(ok(taxtab$Genus), taxtab$Genus, NA))
lab <- ifelse(ok(gen), paste0("g_", cl(gen)), NA)
for (rk in c("Family", "Order", "Class", "Phylum", "Kingdom")) {
  fill <- is.na(lab) & ok(taxtab[[rk]])
  lab[fill] <- paste0(c(Family = "f_", Order = "o_", Class = "c_", Phylum = "p_", Kingdom = "d_")[rk],
                      cl(taxtab[[rk]][fill]))
}
lab[is.na(lab)] <- "d_NA"

# Species-level names used by VALENCIA
lab[grepl("^g_Prevotella(_[0-9]+)?$", lab)] <- "g_Prevotella"
lab[lab == "g_Gardnerella"] <- "Gardnerella_vaginalis"
is_l <- taxtab$Genus_v123 %in% "Lactobacillus" & ok(taxtab$Species)
sp_l <- paste0("Lactobacillus_", c("rhamnosus-casei" = "paracasei_Lactobacillus_rhamnosus",
                                   "reuteri-vaginalis" = "reuteri")[taxtab$Species])
sp_l <- ifelse(is.na(sp_l) | sp_l == "Lactobacillus_NA", paste0("Lactobacillus_", taxtab$Species), sp_l)
lab[is_l & sp_l %in% cent] <- sp_l[is_l & sp_l %in% cent]
lab[taxtab$Genus_v123 %in% "Trichomonas" | taxtab$Kingdom %in% "Eukaryota"] <- "d_Eukaryota"
lab[taxtab$Order %in% "Chloroplast"] <- "c_Chloroplast"

# BVAB1-3 by exact match inside the reference 16S
for (b in c("BVAB1", "BVAB2", "BVAB3")) {
  hit <- vapply(colnames(counts), grepl, logical(1), x = refs$sequence[refs$name == b], fixed = TRUE)
  stopifnot(sum(hit) == 1)
  lab[hit] <- if (b == "BVAB1") "BVAB1" else "BVAB2/BVAB3"
}

vc  <- t(rowsum(t(counts), lab))
vin <- data.frame(sampleID = as.character(meta$SampleID), read_count = rowSums(counts), vc, check.names = FALSE)
write.csv(vin, file.path(d11, "valencia_input.csv"), row.names = FALSE)
write.csv(data.frame(variant = colnames(counts), genus_v123 = taxtab$Genus_v123, species = taxtab$Species,
                     valencia_taxon = lab, reads = colSums(counts)),
          file.path(d11, "valencia_taxon_map.csv"), row.names = FALSE)
