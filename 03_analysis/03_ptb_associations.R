# Continuous-frequency associations with preterm birth: one value per woman, one-sided Wilcoxon within cohort.
# Each test is also run on the deposited table and compared with the published values (pipeline validation).
# In:  data/processed (from 02), data/external
# Out: results/03_ptb_associations/, data/processed/per_woman_frequencies.rds

library(here)
library(ggplot2)
library(patchwork)
library(epitools)
source(here("03_analysis/functions.R"))

out <- here("03_analysis/results/03_ptb_associations")
for (d in c("figures", "tables")) dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)
save_fig <- function(p, name, w, h) for (ext in c("pdf", "png"))
  ggsave(file.path(out, "figures", paste0(name, ".", ext)), p, width = w, height = h, dpi = 300)

counts <- readRDS(here("data/processed/counts_swabs_x_variants.rds"))
taxtab <- readRDS(here("data/processed/taxonomy_variants.rds"))
meta   <- readRDS(here("data/processed/sample_metadata.rds"))
stopifnot(identical(rownames(counts), rownames(meta)), identical(colnames(counts), rownames(taxtab)))
A <- new.env()
load(here("data/external/processed.rda"), envir = A)
refs <- read.delim(here("data/external/reference_16S_BVAB_Mgirerdii.tsv"), stringsAsFactors = FALSE)

ft   <- counts / rowSums(counts)
subj <- as.character(meta$SubjectID)
cols <- c(Preterm = "#C0392B", Term = "#2E86C1")

# Per-woman mean frequency of every genus and the four Lactobacillus species
subdf <- meta[!duplicated(subj), c("SubjectID", "Cohort", "Race", "preterm", "Outcome")]
subdf$SubjectID <- as.character(subdf$SubjectID)
rownames(subdf) <- subdf$SubjectID
subdf$Outcome <- factor(subdf$Outcome, levels = c("Preterm", "Term"))

genus  <- taxtab$Genus_v123
lactos <- c("Lactobacillus crispatus", "Lactobacillus iners", "Lactobacillus jensenii", "Lactobacillus gasseri")
genera <- c(unique(genus[!is.na(genus)]), lactos)
for (g in genera) {
  is_g <- genus %in% g | paste(genus, taxtab$Species) %in% g
  subdf[[g]] <- tapply(rowSums(ft[, is_g, drop = FALSE]), subj, mean)[subdf$SubjectID]
}
saveRDS(subdf, here("data/processed/per_woman_frequencies.rds"))

grp <- function(cohort, ptb) subdf$Cohort == cohort & subdf$preterm == ptb
wil <- function(var, cohort, alt) {
  x <- subdf[subdf$Cohort == cohort &  subdf$preterm, var]
  y <- subdf[subdf$Cohort == cohort & !subdf$preterm, var]
  wilcox_p(x, y, alt)
}
wil_cols <- function(m, cohort) apply(m, 2, function(v) wilcox_p(v[grp(cohort, TRUE)], v[grp(cohort, FALSE)], "greater"))

# Three hypothesised genera (Lactobacillus, Gardnerella, Ureaplasma)
hyp <- data.frame(genus     = c("Lactobacillus", "Gardnerella", "Ureaplasma"),
                  direction = c("lower in PTB", "higher in PTB", "higher in PTB"),
                  alt       = c("less", "greater", "greater"))
hyp$p_Stanford     <- mapply(wil, hyp$genus, "Stanford", hyp$alt)
hyp$p_UAB          <- mapply(wil, hyp$genus, "UAB", hyp$alt)
hyp$paper_Stanford <- c(0.009282, 0.007006, 0.3022)
hyp$paper_UAB      <- c(0.2478, 0.1931, 0.7728)
write.csv(hyp[, c("genus", "direction", "p_Stanford", "paper_Stanford", "p_UAB", "paper_UAB")],
          file.path(out, "tables/hypothesis_genera.csv"), row.names = FALSE)

# Figure: hypothesised genera by outcome
floor <- 1e-6
long <- do.call(rbind, lapply(hyp$genus, function(g)
  data.frame(Cohort = subdf$Cohort, Outcome = subdf$Outcome, Genus = g, Freq = subdf[[g]])))
plab <- rbind(data.frame(Genus = hyp$genus, Cohort = "Stanford", p = hyp$p_Stanford),
              data.frame(Genus = hyp$genus, Cohort = "UAB", p = hyp$p_UAB))
plab$label <- ifelse(plab$p < 0.001, "p < 0.001", sprintf("p = %.3f", plab$p))

panel <- function(g, log_scale) {
  d <- long[long$Genus == g, ]
  d$y <- if (log_scale) pmax(d$Freq, floor) else d$Freq
  lab <- plab[plab$Genus == g, ]
  lab$y <- if (log_scale) 3 else 1.07
  p <- ggplot(d, aes(Outcome, y, colour = Outcome)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.35) +
    geom_jitter(width = 0.15, height = 0, size = 1.1, alpha = 0.7) +
    geom_text(data = lab, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 3) +
    facet_grid(Cohort ~ .) +
    scale_colour_manual(values = cols) +
    labs(title = g, x = NULL, y = "Mean gestational frequency (per woman)") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", plot.title = element_text(face = "italic", hjust = 0.5),
          panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
  if (log_scale) {
    p + scale_y_log10(breaks = c(1e-6, 1e-4, 1e-2, 1), labels = c("0", "0.0001", "0.01", "1")) +
      coord_cartesian(ylim = c(floor / 2, 6))
  } else {
    p + scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1.12))
  }
}

fig1 <- panel("Lactobacillus", FALSE) +
  panel("Gardnerella", TRUE) + labs(y = NULL) +
  panel("Ureaplasma", TRUE) + labs(y = NULL) +
  plot_annotation(
    title    = "Average gestational frequency by birth outcome",
    subtitle = "Stanford: 39 women (9 preterm) | UAB: 96 women (41 preterm)",
    caption  = paste("One point per woman. One-sided Wilcoxon rank-sum test within each cohort",
                     "(Lactobacillus: lower in preterm; Gardnerella, Ureaplasma: higher in preterm).",
                     "\nLog panels: women with zero frequency are drawn at the bottom (labelled 0)."))
save_fig(fig1, "genus_by_outcome", 10, 6.5)

# Odds ratios at the published thresholds (Lactobacillus < 70%, Gardnerella > 0.1%)
odds_it <- function(d, taxon, thresh, risk = c("greater", "less")) {
  risk <- match.arg(risk)
  exposed <- if (risk == "greater") d[[taxon]] > thresh else d[[taxon]] < thresh
  tab <- table(exposed, ordered(d$Outcome, levels = c("Term", "Preterm")))
  r <- oddsratio(tab)
  c(OR = r$measure["TRUE", "estimate"], lower = r$measure["TRUE", "lower"],
    upper = r$measure["TRUE", "upper"], p_fisher = r$p.value["TRUE", "fisher.exact"],
    n_exposed_PTB = tab["TRUE", "Preterm"], n_exposed_term = tab["TRUE", "Term"])
}
ors <- as.data.frame(rbind(
  "Lactobacillus < 70% | Stanford" = odds_it(subdf[subdf$Cohort == "Stanford", ], "Lactobacillus", 0.7, "less"),
  "Lactobacillus < 70% | UAB"      = odds_it(subdf[subdf$Cohort == "UAB", ], "Lactobacillus", 0.7, "less"),
  "Gardnerella > 0.1% | Stanford"  = odds_it(subdf[subdf$Cohort == "Stanford", ], "Gardnerella", 0.001, "greater"),
  "Gardnerella > 0.1% | UAB"       = odds_it(subdf[subdf$Cohort == "UAB", ], "Gardnerella", 0.001, "greater")))
ors$paper_OR    <- c(5.812704, 1.296377, 5.123754, 1.57168)
ors$paper_lower <- c(1.122489, 0.5681911, 1.051368, 0.5012872)
ors$paper_upper <- c(33.74815, 2.999572, 31.12366, 5.575016)
write.csv(cbind(test = rownames(ors), ors), file.path(out, "tables/odds_ratios.csv"), row.names = FALSE)

fp <- data.frame(test   = rep(rownames(ors), 2),
                 source = rep(c("Ours", "Paper"), each = nrow(ors)),
                 OR     = c(ors$OR, ors$paper_OR),
                 lower  = c(ors$lower, ors$paper_lower),
                 upper  = c(ors$upper, ors$paper_upper))
fp$taxon  <- sub(" \\|.*", "", fp$test)
fp$cohort <- sub(".*\\| ", "", fp$test)
fig_or <- ggplot(fp, aes(x = OR, y = cohort, colour = source)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 2.3, position = position_dodge(width = 0.5)) +
  facet_wrap(~ taxon, ncol = 1) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 30)) +
  scale_colour_manual(values = c(Ours = "#2C3E50", Paper = "#E67E22")) +
  labs(title = "Odds of preterm birth in the high-risk group",
       x = "Odds ratio (95% CI, log scale)", y = NULL, colour = NULL,
       caption = "High risk: Lactobacillus < 70% or Gardnerella > 0.1% (per-woman average). epitools::oddsratio, median-unbiased estimate.") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top", strip.text = element_text(face = "italic"))
save_fig(fig_or, "odds_ratios", 7, 4.5)

# Per-variant tests, ours and on the deposited table; Fisher-combined p, BH FDR < 0.1
n_per_woman <- as.vector(table(subj)[subdf$SubjectID])
sv_mean <- rowsum(as.matrix(ft), subj)[subdf$SubjectID, ] / n_per_woman
p_stan <- wil_cols(sv_mean, "Stanford"); p_stan[is.na(p_stan)] <- 1
p_uab  <- wil_cols(sv_mean, "UAB");      p_uab[is.na(p_uab)]   <- 1

taxdf <- data.frame(pStanford = p_stan, pUAB = p_uab, Genus = taxtab$Genus_v123, Gvariant = taxtab$Gvariant,
                    freqStanford = colMeans(ft[meta$Cohort == "Stanford", ]),
                    freqUAB      = colMeans(ft[meta$Cohort == "UAB", ]),
                    row.names = colnames(ft))
taxdf$Joint <- pchisq(-2 * log(taxdf$pStanford * taxdf$pUAB), df = 4, lower.tail = FALSE)

t_subj <- as.character(A$df$SubjectID)
t_mean <- rowsum(A$ft, t_subj)[subdf$SubjectID, ] / as.vector(table(t_subj)[subdf$SubjectID])
pt_stan <- wil_cols(t_mean, "Stanford")
pt_uab  <- wil_cols(t_mean, "UAB")

sig_ours_S   <- rownames(taxdf)[p.adjust(taxdf$pStanford, "BH") < 0.1]
sig_ours_U   <- rownames(taxdf)[p.adjust(taxdf$pUAB, "BH") < 0.1]
sig_theirs_S <- names(pt_stan)[p.adjust(pt_stan, "BH") < 0.1]
sig_theirs_U <- names(pt_uab)[p.adjust(pt_uab, "BH") < 0.1]

full <- data.frame(variant = rownames(taxdf),
                   taxdf[, c("Genus", "Gvariant", "pStanford", "pUAB", "Joint", "freqStanford", "freqUAB")],
                   BH_Stanford = p.adjust(taxdf$pStanford, "BH"),
                   BH_UAB      = p.adjust(taxdf$pUAB, "BH"), row.names = NULL)
write.csv(full, file.path(out, "tables/variant_tests.csv"), row.names = FALSE)

hits <- union(union(sig_ours_S, sig_theirs_S), union(sig_ours_U, sig_theirs_U))
cmp83 <- data.frame(Genus = taxdf[hits, "Genus"],
                    p_ours_S = signif(taxdf[hits, "pStanford"], 3), p_theirs_S = signif(pt_stan[hits], 3),
                    sig_ours_S = hits %in% sig_ours_S, sig_theirs_S = hits %in% sig_theirs_S,
                    p_ours_U = signif(taxdf[hits, "pUAB"], 3), p_theirs_U = signif(pt_uab[hits], 3),
                    sig_ours_U = hits %in% sig_ours_U, sig_theirs_U = hits %in% sig_theirs_U,
                    freqStanford = signif(taxdf[hits, "freqStanford"], 3),
                    variant = hits, row.names = NULL)
cmp83 <- cmp83[order(cmp83$p_ours_S), ]
write.csv(cmp83, file.path(out, "tables/variant_significant_ours_vs_deposit.csv"), row.names = FALSE)

summ83 <- data.frame(source   = c("Paper (printed, R 3.4)", "Authors' deposit, re-run in R 4.x", "Our pipeline"),
                     Stanford = c(43, length(sig_theirs_S), length(sig_ours_S)),
                     UAB      = c(1, length(sig_theirs_U), length(sig_ours_U)))
write.csv(summ83, file.path(out, "tables/variant_significant_counts.csv"), row.names = FALSE)

# Lactobacillus species and Gardnerella variants
for (gv in c("G1", "G2", "G3"))
  subdf[[paste("Gardnerella", gv)]] <-
    tapply(rowSums(ft[, taxtab$Gvariant %in% gv, drop = FALSE]), subj, mean)[subdf$SubjectID]
taxa85 <- c(lactos, paste("Gardnerella", c("G1", "G2", "G3")))
labs85 <- sub("Lactobacillus", "L.", taxa85)

long85 <- do.call(rbind, lapply(taxa85, function(tx)
  data.frame(Cohort = subdf$Cohort, Outcome = subdf$Outcome, Taxon = tx, Freq = subdf[[tx]])))
long85$Taxon <- factor(long85$Taxon, levels = taxa85, labels = labs85)
long85$Plot  <- pmax(long85$Freq, 1e-6)

pv85 <- expand.grid(Taxon = taxa85, Cohort = c("Stanford", "UAB"), stringsAsFactors = FALSE)
pv85$p <- mapply(function(tx, co) wil(tx, co, if (grepl("^Lactobacillus", tx)) "less" else "greater"),
                 pv85$Taxon, pv85$Cohort)
pv85$lab   <- paste0("p = ", signif(pv85$p, 2))
pv85$Taxon <- factor(pv85$Taxon, levels = taxa85, labels = labs85)

fig2a <- ggplot(long85, aes(Outcome, Plot, colour = Outcome)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, height = 0, size = 1.3, alpha = 0.7) +
  geom_text(data = pv85, aes(x = 1.5, y = 3, label = lab), inherit.aes = FALSE, size = 3) +
  facet_grid(Cohort ~ Taxon) +
  scale_y_log10(breaks = c(1e-6, 1e-4, 1e-2, 1), labels = c("0", expression(10^-4), expression(10^-2), "1"),
                limits = c(1e-6, 6)) +
  scale_colour_manual(values = cols) +
  labs(x = NULL, y = "Mean gestational frequency (per woman)",
       title = "Lactobacillus species and Gardnerella variants by birth outcome",
       subtitle = "One-sided Wilcoxon: Lactobacillus 'less' in preterm, Gardnerella 'greater' in preterm") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", strip.text.x = element_text(face = "italic"),
        axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(fig2a, "lacto_species_gardnerella_variants", 13, 7)

# All Gardnerella variants (G1..G9 by abundance) and the genus without G2
gard_cols <- colnames(ft)[taxtab$Genus_v123 %in% "Gardnerella"]
gard_cols <- gard_cols[order(-colSums(counts[, gard_cols]))]
stopifnot(identical(taxtab[gard_cols[1:3], "Gvariant"], c("G1", "G2", "G3")))
gv_mean <- sv_mean[, gard_cols, drop = FALSE]
colnames(gv_mean) <- paste0("G", seq_along(gard_cols))
gv_mean <- cbind(gv_mean, `Gardnerella minus G2` = rowSums(sv_mean[, gard_cols[-2], drop = FALSE]))

lac_S <- grepl("^L\\.", pv85$Taxon) & pv85$Cohort == "Stanford"
lac_U <- grepl("^L\\.", pv85$Taxon) & pv85$Cohort == "UAB"
tab85 <- rbind(
  data.frame(Taxon = colnames(gv_mean), pStanford = wil_cols(gv_mean, "Stanford"),
             pUAB = wil_cols(gv_mean, "UAB"), women_with_it = colSums(gv_mean > 0), row.names = NULL),
  data.frame(Taxon = pv85$Taxon[lac_S], pStanford = pv85$p[lac_S], pUAB = pv85$p[lac_U], women_with_it = NA))
tab85$pStanford <- signif(tab85$pStanford, 3)
tab85$pUAB      <- signif(tab85$pUAB, 3)
write.csv(tab85, file.path(out, "tables/species_gardnerella_variant_tests.csv"), row.names = FALSE)

# Genus-level tests (genera + 4 Lactobacillus species), ours and on the deposited table
lactos_less <- c("Lactobacillus", lactos)
t_gen <- A$tax[colnames(A$ft), "Genus"]
t_spc <- A$tax[colnames(A$ft), "Species"]
t_genera <- c(unique(t_gen[!is.na(t_gen)]), lactos)
t_sub <- sapply(t_genera, function(g) {
  is_g <- t_gen %in% g | paste(t_gen, t_spc) %in% g
  tapply(rowSums(A$ft[, is_g, drop = FALSE]), t_subj, mean)[subdf$SubjectID]
})
o_sub <- as.matrix(subdf[, genera])

gtest <- function(m, cohort) sapply(colnames(m), function(g)
  wilcox_p(m[grp(cohort, TRUE), g], m[grp(cohort, FALSE), g], if (g %in% lactos_less) "less" else "greater"))
res84 <- lapply(list(ours = o_sub, theirs = t_sub), function(m) {
  pS <- gtest(m, "Stanford"); pU <- gtest(m, "UAB")
  list(pS = pS, pU = pU, sigS = names(pS)[p.adjust(pS, "BH") < 0.1], sigU = names(pU)[p.adjust(pU, "BH") < 0.1])
})

summ84 <- data.frame(source   = c("Paper (printed, R 3.4)", "Authors' deposit, re-run in R 4.x", "Our pipeline"),
                     n_tests  = c(349, ncol(t_sub), ncol(o_sub)),
                     Stanford = c(44, length(res84$theirs$sigS), length(res84$ours$sigS)),
                     UAB      = c(2, length(res84$theirs$sigU), length(res84$ours$sigU)))
write.csv(summ84, file.path(out, "tables/genus_significant_counts.csv"), row.names = FALSE)

hits84 <- union(union(res84$ours$sigS, res84$theirs$sigS), union(res84$ours$sigU, res84$theirs$sigU))
cmp84 <- data.frame(Genus = hits84,
                    p_ours_S = signif(res84$ours$pS[hits84], 3), p_theirs_S = signif(res84$theirs$pS[hits84], 3),
                    sig_ours_S = hits84 %in% res84$ours$sigS, sig_theirs_S = hits84 %in% res84$theirs$sigS,
                    p_ours_U = signif(res84$ours$pU[hits84], 3), p_theirs_U = signif(res84$theirs$pU[hits84], 3),
                    sig_ours_U = hits84 %in% res84$ours$sigU, sig_theirs_U = hits84 %in% res84$theirs$sigU,
                    row.names = NULL)
cmp84 <- cmp84[order(cmp84$p_ours_S), ]
write.csv(cmp84, file.path(out, "tables/genus_significant_ours_vs_deposit.csv"), row.names = FALSE)

genFreq <- function(g, coh_rows) mean(rowSums(ft[coh_rows, genus %in% g | paste(genus, taxtab$Species) %in% g, drop = FALSE]))
gendf <- data.frame(Genus = genera, pStanford = res84$ours$pS[genera], pUAB = res84$ours$pU[genera],
                    freqStanford = sapply(genera, genFreq, coh_rows = meta$Cohort == "Stanford"),
                    freqUAB      = sapply(genera, genFreq, coh_rows = meta$Cohort == "UAB"),
                    row.names = genera)
gendf$Joint <- pchisq(-2 * log(gendf$pStanford * gendf$pUAB), df = 4, lower.tail = FALSE)
gendf$Sig   <- p.adjust(gendf$Joint, "BH") < 0.1
gendf$Freq  <- (gendf$freqStanford + gendf$freqUAB) / 2
write.csv(gendf[, c("Genus", "pStanford", "pUAB", "Joint", "Sig", "freqStanford", "freqUAB", "Freq")],
          file.path(out, "tables/genus_tests.csv"), row.names = FALSE)

t_joint <- pchisq(-2 * log(res84$theirs$pS * res84$theirs$pU), df = 4, lower.tail = FALSE)
summ87 <- data.frame(source = c("Authors' deposit, re-run in R 4.x", "Our pipeline"),
                     joint_BH_sig = c(sum(p.adjust(t_joint, "BH") < 0.1), sum(gendf$Sig)))
write.csv(summ87, file.path(out, "tables/genus_joint_significant_counts.csv"), row.names = FALSE)

# Figure: every genus, Stanford vs UAB p-values
fig4a <- ggplot(gendf, aes(-log10(pStanford), -log10(pUAB), label = Genus, colour = Sig)) +
  geom_text(aes(size = sqrt(Freq)), show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "grey60") +
  scale_size_continuous(range = c(1.5, 7)) +
  scale_colour_manual(values = c(`FALSE` = "grey35", `TRUE` = "#C0392B")) +
  scale_x_continuous(expand = expansion(mult = c(0.25, 0.05))) +
  labs(x = expression(-log[10]~P~"(Stanford)"), y = expression(-log[10]~P~"(UAB)"),
       title = "Association of each genus with preterm birth, by cohort",
       subtitle = "One-sided Wilcoxon on per-woman mean frequency; red = Fisher-combined BH FDR < 0.1; text size ~ sqrt(study-wide frequency)") +
  coord_equal() + theme_bw(base_size = 11)
save_fig(fig4a, "genus_overview", 10, 9)

# BV-associated taxa (Table 1): BVAB1-3 by exact match inside the reference 16S, plus BV genera
bv <- sapply(c("BVAB1", "BVAB2", "BVAB3"), function(b) {
  hit <- colnames(ft)[vapply(colnames(ft), grepl, logical(1), x = refs$sequence[refs$name == b], fixed = TRUE)]
  stopifnot(length(hit) == 1)
  hit
})
bv_genera <- c("Gardnerella", "Prevotella", "Atopobium", "Mobiluncus", "Megasphaera", "Dialister",
               "Peptoniphilus", "Mycoplasma")
t1 <- rbind(
  data.frame(Taxon = names(bv), taxdf[bv, c("pStanford", "pUAB", "freqStanford", "freqUAB")], Resolution = "SV"),
  data.frame(Taxon = bv_genera, gendf[bv_genera, c("pStanford", "pUAB", "freqStanford", "freqUAB")], Resolution = "Genus"))
table1 <- data.frame(Taxon = t1$Taxon,
                     pS_ours = signif(t1$pStanford, 3),
                     pS_paper = c(0.907, 0.823, 0.729, 0.007, 0.008, 0.043, 0.034, 0.108, 0.012, 0.016, 0.042),
                     pU_ours = signif(t1$pUAB, 3),
                     pU_paper = c(0.339, 0.475, 0.280, 0.193, 0.024, 0.476, 0.619, 0.320, 0.140, 0.843, 0.283),
                     fS_ours = signif(t1$freqStanford, 2),
                     fS_paper = c(0.0017, 0.00049, 1.4e-06, 0.051, 0.023, 0.0072, 0.00017, 0.010, 0.0033, 0.0038, 0.00058),
                     fU_ours = signif(t1$freqUAB, 2),
                     fU_paper = c(0.034, 0.00081, 3.9e-05, 0.083, 0.062, 0.012, 0.00059, 0.054, 0.0069, 0.0089, 0.0058),
                     Resolution = t1$Resolution,
                     pS_deposit = signif(c(pt_stan[bv], res84$theirs$pS[bv_genera]), 3),
                     pU_deposit = signif(c(pt_uab[bv], res84$theirs$pU[bv_genera]), 3))
write.csv(table1, file.path(out, "tables/bv_taxa.csv"), row.names = FALSE)
