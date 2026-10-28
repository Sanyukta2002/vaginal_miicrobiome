# Community ecology: cohort composition, L. crispatus / L. iners vs Gardnerella co-occurrence,
# week-to-week instability, Bray-Curtis MDS, time-window tests and STIs.
# In:  data/processed (from 02), data/external
# Out: results/04_community_ecology/

library(here)
library(ggplot2)
library(vegan)
source(here("03_analysis/functions.R"))

out <- here("03_analysis/results/04_community_ecology")
for (d in c("figures", "tables")) dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)
save_fig <- function(p, name, w, h) for (ext in c("pdf", "png"))
  ggsave(file.path(out, "figures", paste0(name, ".", ext)), p, width = w, height = h, dpi = 300)

counts <- readRDS(here("data/processed/counts_swabs_x_variants.rds"))
taxtab <- readRDS(here("data/processed/taxonomy_variants.rds"))
meta   <- readRDS(here("data/processed/sample_metadata.rds"))
stopifnot(nrow(counts) == nrow(meta), identical(colnames(counts), rownames(taxtab)),
          all(c("GW", "GDcoll") %in% names(meta)))
A <- new.env()
load(here("data/external/processed.rda"), envir = A)
refs <- read.delim(here("data/external/reference_16S_BVAB_Mgirerdii.tsv"), stringsAsFactors = FALSE)

ft    <- counts / rowSums(counts)
genus <- taxtab$Genus_v123
spc   <- paste(genus, taxtab$Species)
subj  <- as.character(meta$SubjectID)

# Per-swab frequencies of the key taxa
key <- data.frame(
  SampleID = rownames(ft), SubjectID = subj, Cohort = meta$Cohort,
  Outcome = factor(meta$Outcome, levels = c("Preterm", "Term")),
  Lactobacillus = rowSums(ft[, genus %in% "Lactobacillus", drop = FALSE]),
  crispatus   = rowSums(ft[, spc %in% "Lactobacillus crispatus", drop = FALSE]),
  iners       = rowSums(ft[, spc %in% "Lactobacillus iners", drop = FALSE]),
  jensenii    = rowSums(ft[, spc %in% "Lactobacillus jensenii", drop = FALSE]),
  gasseri     = rowSums(ft[, spc %in% "Lactobacillus gasseri", drop = FALSE]),
  Gardnerella = rowSums(ft[, genus %in% "Gardnerella", drop = FALSE]),
  G1 = rowSums(ft[, taxtab$Gvariant %in% "G1", drop = FALSE]),
  G2 = rowSums(ft[, taxtab$Gvariant %in% "G2", drop = FALSE]),
  G3 = rowSums(ft[, taxtab$Gvariant %in% "G3", drop = FALSE]))

# Per-woman means
subdf <- meta[!duplicated(subj), c("SubjectID", "Cohort", "Race", "preterm", "Outcome")]
subdf$SubjectID <- as.character(subdf$SubjectID)
rownames(subdf) <- subdf$SubjectID
subdf$Outcome <- factor(subdf$Outcome, levels = c("Preterm", "Term"))
for (v in c("Lactobacillus", "crispatus", "iners", "Gardnerella", "G1", "G2"))
  subdf[[v]] <- tapply(key[[v]], subj, mean)[subdf$SubjectID]

# Cohort composition
key$L_other <- key$Lactobacillus - key$crispatus - key$iners
key$G_other <- key$Gardnerella - key$G1 - key$G2
comp <- aggregate(key[, c("crispatus", "iners", "L_other", "G1", "G2", "G_other")], list(Cohort = key$Cohort), mean)
comp$Gardnerella_total <- rowSums(comp[, c("G1", "G2", "G_other")])
write.csv(comp, file.path(out, "tables/cohort_composition.csv"), row.names = FALSE)

taxa3a <- c("L. crispatus", "L. iners", "L. (other)", "G1", "G2", "G (other)")
c3 <- data.frame(Cohort = rep(comp$Cohort, 6), Taxa = factor(rep(taxa3a, each = 2), levels = rev(taxa3a)),
                 Frequency = unlist(comp[, 2:7]))
fig3a <- ggplot(c3, aes(Cohort, Frequency, fill = Taxa)) + geom_col(width = 0.8) +
  scale_fill_manual(values = c("L. crispatus" = "#F8766D", "L. iners" = "#B79F00", "L. (other)" = "#00BA38",
                               "G1" = "#00BFC4", "G2" = "#619CFF", "G (other)" = "#F564E3")) +
  labs(y = "Mean frequency (pooled swabs)", x = NULL) + theme_bw()
save_fig(fig3a, "cohort_composition", 3.35, 3.9)

# Co-occurrence: Lactobacillus : Gardnerella ratio vs their summed frequency
co <- rbind(
  data.frame(Cohort = key$Cohort, Lacto = "L. crispatus", sum = key$crispatus + key$Gardnerella,
             ratio = key$crispatus / key$Gardnerella),
  data.frame(Cohort = key$Cohort, Lacto = "L. iners", sum = key$iners + key$Gardnerella,
             ratio = key$iners / key$Gardnerella))
co$coexist <- co$ratio >= 0.1 & co$ratio <= 10
co_tab <- aggregate(cbind(n_swabs = sum > 0.5, coexist = sum > 0.5 & coexist) ~ Lacto + Cohort, data = co, FUN = sum)
co_tab$pct_coexist <- round(100 * co_tab$coexist / co_tab$n_swabs, 1)
write.csv(co_tab, file.path(out, "tables/cooccurrence_coexistence.csv"), row.names = FALSE)

fig3b <- ggplot(subset(co, is.finite(ratio) & ratio > 0), aes(sum, ratio)) +
  geom_point(shape = 16, alpha = 0.3) + facet_grid(Lacto ~ Cohort) +
  scale_y_log10(breaks = 10^c(-4, -2, 0, 2, 4), labels = c("1:10000", "1:100", "1:1", "100:1", "10000:1")) +
  labs(x = "Summed frequency (Lactobacillus species + Gardnerella)", y = "Lactobacillus : Gardnerella") +
  theme_bw() + theme(aspect.ratio = 1)
save_fig(fig3b, "cooccurrence", 6, 5.5)

# Instability: mean Bray-Curtis between consecutive swabs 4-10 days apart
inst <- sapply(split(seq_len(nrow(meta)), subj), function(ix) {
  ix <- ix[order(meta$GDcoll[ix])]
  if (length(ix) < 2) return(NA)
  a <- ix[-length(ix)]; b <- ix[-1]
  keep <- (meta$GDcoll[b] - meta$GDcoll[a]) %in% 4:10
  if (!any(keep)) return(NA)
  mean(0.5 * rowSums(abs(ft[a[keep], , drop = FALSE] - ft[b[keep], , drop = FALSE])))
})
subdf$AveInst <- inst[subdf$SubjectID]

inst_tab <- data.frame(
  Cohort = c("Stanford", "UAB"),
  women_with_value = sapply(c("Stanford", "UAB"), function(co) sum(!is.na(subdf$AveInst[subdf$Cohort == co]))),
  p_PTB = sapply(c("Stanford", "UAB"), function(co)
    with(subdf[subdf$Cohort == co, ], wilcox_p(AveInst[!preterm], AveInst[preterm]))),
  rho_vs_Lacto = sapply(c("Stanford", "UAB"), function(co) with(subdf[subdf$Cohort == co, ],
    cor(Lactobacillus, AveInst, method = "spearman", use = "complete.obs"))),
  paper_p = c(0.0018, NA), paper_rho = c(-0.76, NA))
write.csv(inst_tab, file.path(out, "tables/instability.csv"), row.names = FALSE)

# Bray-Curtis MDS with landmark swabs (highest L. crispatus, L. iners, G1, G2)
mds <- cmdscale(vegdist(counts, "bray"), k = 2, eig = TRUE)
rel <- round(100 * mds$eig[1:2] / sum(mds$eig[mds$eig > 0]), 1)
lm_i <- c("L. crispatus" = which.max(key$crispatus), "L. iners" = which.max(key$iners),
          G1 = which.max(key$G1), G2 = which.max(key$G2))
ord <- data.frame(Axis.1 = mds$points[, 1], Axis.2 = mds$points[, 2], Cohort = key$Cohort)
lmk <- do.call(rbind, lapply(c("Stanford", "UAB"), function(co)
  data.frame(ord[lm_i, 1:2], Cohort = co, label = names(lm_i))))
lmk$hjust <- ifelse(lmk$label == "L. crispatus", 0.15, 0.8)
fig2c <- ggplot(ord, aes(Axis.1, Axis.2)) +
  geom_point(shape = 16, alpha = 0.2, size = 0.8) +
  geom_point(data = lmk, colour = "magenta", size = 2) +
  geom_text(data = lmk, aes(label = label, hjust = hjust), size = 3, vjust = -0.6) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.12))) +
  facet_grid(Cohort ~ .) + coord_equal() +
  labs(x = paste0("Axis 1 [", rel[1], "%]"), y = paste0("Axis 2 [", rel[2], "%]")) + theme_bw()
save_fig(fig2c, "mds_landmarks", 5, 7)

# Tests within gestational-age windows
meta$Interval <- cut(meta$GW, c(12, 18, 24, 30, 36))
tw <- do.call(rbind, lapply(c("Lactobacillus", "crispatus", "iners", "Gardnerella", "G1", "G2"), function(tx)
  do.call(rbind, lapply(levels(meta$Interval), function(int) {
    ii  <- meta$Interval %in% int
    f   <- tapply(key[[tx]][ii], subj[ii], mean)[subdf$SubjectID]
    alt <- if (tx %in% c("Lactobacillus", "crispatus", "iners")) "less" else "greater"
    data.frame(taxa = tx, weeks = int, cohort = c("Stanford", "UAB"),
               pval = sapply(c("Stanford", "UAB"), function(co)
                 wilcox_p(f[subdf$Cohort == co & subdf$preterm], f[subdf$Cohort == co & !subdf$preterm], alt)))
  }))))
tw$weeks <- factor(tw$weeks, levels = levels(meta$Interval), labels = c("13-18", "19-24", "25-30", "31-36"))
write.csv(tw, file.path(out, "tables/time_windows.csv"), row.names = FALSE)

figS3 <- ggplot(tw, aes(weeks, pval, colour = cohort)) + geom_point(size = 2) + facet_wrap(~taxa) +
  geom_hline(yintercept = 0.05, linetype = "dashed") + scale_y_log10() +
  scale_colour_manual(values = c(Stanford = "#8C1515", UAB = "#F4C300")) +
  labs(x = "Gestational weeks", y = "p-value (one-sided Wilcoxon)") + theme_bw()
save_fig(figS3, "time_windows", 7, 5)

# Sexually transmitted infections
gon <- colnames(A$st)[A$tax[colnames(A$st), "Genus"] %in% "Neisseria"][3]
myc <- refs$sequence[refs$name == "Mgirerdii"]
sti <- list(Chlamydia   = genus %in% "Chlamydia",
            Gonorrhea   = colnames(ft) %in% gon,
            Trichomonas = genus %in% "Trichomonas",
            M.girerdii  = vapply(colnames(ft), function(v)
              grepl(v, myc, fixed = TRUE) || grepl(myc, v, fixed = TRUE), logical(1)))
sti_tab <- do.call(rbind, lapply(names(sti), function(s) {
  f <- tapply(rowSums(ft[, sti[[s]], drop = FALSE]), subj, mean)[subdf$SubjectID]
  data.frame(STI = s, n_variants = sum(sti[[s]]),
             women_Stanford = sum(f[subdf$Cohort == "Stanford"] > 0),
             women_UAB = sum(f[subdf$Cohort == "UAB"] > 0),
             p_PTB = wilcox_p(f[subdf$preterm], f[!subdf$preterm], "greater"))
}))
write.csv(sti_tab, file.path(out, "tables/stis.csv"), row.names = FALSE)
