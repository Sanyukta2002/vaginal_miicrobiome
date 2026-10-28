# Compare VALENCIA CSTs with continuous frequencies of L. crispatus, L. iners and Gardnerella:
# (1) how much frequency variation the labels keep, (2) the CST III/IV boundary,
# (3) PTB associations detected by each representation.
# In:  data/processed (from 02), data/processed/valencia (from 05-06)
# Out: results/07_cst_vs_continuous/, data/processed/valencia/{swabs,women}_cst_and_frequencies.rds

library(here)
library(ggplot2)
library(vegan)
source(here("03_analysis/functions.R"))

d11 <- here("data/processed/valencia")
out <- here("03_analysis/results/07_cst_vs_continuous")
for (d in c("figures", "tables")) dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)
save_fig <- function(p, name, w, h) for (ext in c("pdf", "png"))
  ggsave(file.path(out, "figures", paste0(name, ".", ext)), p, width = w, height = h, dpi = 300)

counts <- readRDS(here("data/processed/counts_swabs_x_variants.rds"))
taxtab <- readRDS(here("data/processed/taxonomy_variants.rds"))
meta   <- readRDS(here("data/processed/sample_metadata.rds"))
vfull  <- read.csv(file.path(d11, "valencia_output.csv"), check.names = FALSE, colClasses = c(sampleID = "character"))
vin    <- read.csv(file.path(d11, "valencia_input.csv"), check.names = FALSE, colClasses = c(sampleID = "character"))
vfull  <- vfull[match(as.character(meta$SampleID), vfull$sampleID), ]
vin    <- vin[match(as.character(meta$SampleID), vin$sampleID), ]
stopifnot(identical(vfull$sampleID, as.character(meta$SampleID)), identical(vin$sampleID, vfull$sampleID))
vout <- vfull[, c("sampleID", "subCST", "CST", "score")]

cst_col <- c("I-A" = "#ff6868", "I-B" = "#ffd4da", "II" = "#b4ff68", "III-A" = "#ffbc6b", "III-B" = "#e4a67b",
             "IV-A" = "#c1adec", "IV-B" = "#91a8ed", "IV-C0" = "#989898", "IV-C1" = "#ffc0cb", "IV-C2" = "#a8e5e5",
             "IV-C3" = "#9acc9a", "IV-C4" = "#800080", "V" = "#ffff71")
grp34 <- function(x) ifelse(grepl("^III", x), "III", ifelse(grepl("^IV", x), "IV", "other"))

# Bray-Curtis MDS coloured by subCST
mds <- cmdscale(vegdist(counts, "bray"), k = 2, eig = TRUE)
rel <- round(100 * mds$eig[1:2] / sum(mds$eig[mds$eig > 0]), 1)
ordc <- data.frame(Axis.1 = mds$points[, 1], Axis.2 = mds$points[, 2], Cohort = meta$Cohort,
                   subCST = factor(vout$subCST, levels = names(cst_col)))
fig_ord <- ggplot(ordc, aes(Axis.1, Axis.2, fill = subCST)) +
  geom_point(shape = 21, size = 1.8, stroke = 0.15, colour = "grey25", alpha = 0.85) +
  scale_fill_manual(values = cst_col, drop = FALSE) +
  facet_grid(Cohort ~ .) + coord_equal() +
  labs(x = paste0("Axis 1 [", rel[1], "%]"), y = paste0("Axis 2 [", rel[2], "%]"),
       fill = "VALENCIA\nsubCST", title = "Bray-Curtis MDS coloured by VALENCIA community state type") +
  guides(fill = guide_legend(override.aes = list(size = 3))) + theme_bw()
save_fig(fig_ord, "mds_by_subcst", 7, 9)

# Per-swab composition ordered by subCST and its leading taxon
rab <- as.matrix(vin[, -(1:2)]) / vin$read_count
top <- c("Lactobacillus_crispatus", "Lactobacillus_iners", "Lactobacillus_jensenii", "Lactobacillus_gasseri",
         "Gardnerella_vaginalis", "BVAB1", "g_Atopobium", "g_Prevotella", "g_Sneathia", "g_Megasphaera", "g_Streptococcus")
rt <- cbind(rab[, top], Other = pmax(0, 1 - rowSums(rab[, top])))
subf <- factor(vout$subCST, levels = names(cst_col))
lead <- sapply(levels(subf), function(s)
  if (any(subf == s)) names(which.max(colMeans(rt[subf == s, , drop = FALSE]))) else NA)
lv <- rt[cbind(seq_len(nrow(rt)), match(lead[as.character(subf)], colnames(rt)))]
o  <- order(subf, -lv)

bar <- data.frame(pos = rep(seq_along(o), ncol(rt)), subCST = rep(subf[o], ncol(rt)),
                  taxon = factor(rep(colnames(rt), each = length(o)), levels = rev(colnames(rt))),
                  value = as.vector(rt[o, ]))
bar$panel <- factor(ifelse(bar$subCST %in% c("IV-C1", "IV-C2", "IV-C3", "IV-C4"), "IV-C1–4", as.character(bar$subCST)),
                    levels = c("I-A", "I-B", "II", "III-A", "III-B", "IV-A", "IV-B", "IV-C0", "IV-C1–4", "V"))
tax_col <- c(Lactobacillus_crispatus = "#ff6868", Lactobacillus_iners = "#ffbc6b", Lactobacillus_jensenii = "#e6c200",
             Lactobacillus_gasseri = "#7fd13b", Gardnerella_vaginalis = "#4a6fe3", BVAB1 = "#9b59b6",
             g_Atopobium = "#1abc9c", g_Prevotella = "#6d4c41", g_Sneathia = "#d35400", g_Megasphaera = "#2c3e50",
             g_Streptococcus = "#ff9ff3", Other = "#d9d9d9")
tax_lab <- sub("_", " ", sub("^g_", "", names(tax_col)))
tax_lab <- sub("^Lactobacillus ", "L. ", sub("^Gardnerella vaginalis$", "Gardnerella", tax_lab))
fig_bar <- ggplot(bar, aes(pos, value, fill = taxon)) +
  geom_col(width = 1) +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = tax_col, breaks = names(tax_col), labels = tax_lab) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Swabs (n = 2,179), ordered by VALENCIA subCST", y = "Relative abundance", fill = NULL,
       title = "Per-swab composition within each VALENCIA community state type") +
  theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing = unit(1, "pt"),
                     strip.text.x = element_text(size = 8), legend.text = element_text(face = "italic"))
save_fig(fig_bar, "composition_by_subcst", 14, 6)

# Swab table: labels and continuous key-taxon frequencies
ft  <- counts / rowSums(counts)
spc <- paste(taxtab$Genus_v123, taxtab$Species)
sw  <- data.frame(vout, SubjectID = as.character(meta$SubjectID), Cohort = meta$Cohort,
                  preterm = meta$preterm, GDcoll = meta$GDcoll,
                  crispatus   = rowSums(ft[, spc %in% "Lactobacillus crispatus", drop = FALSE]),
                  iners       = rowSums(ft[, spc %in% "Lactobacillus iners", drop = FALSE]),
                  Gardnerella = rowSums(ft[, taxtab$Genus_v123 %in% "Gardnerella", drop = FALSE]))
sw$subCST <- factor(sw$subCST, levels = names(cst_col))

# (1) Information kept: frequencies within each label, and eta^2 (rank and raw scale)
q <- function(x) sprintf("%.3f [%.3f-%.3f]", median(x), quantile(x, .25), quantile(x, .75))
agree <- do.call(rbind, lapply(split(sw, sw$subCST, drop = TRUE), function(d)
  data.frame(subCST = d$subCST[1], n = nrow(d), crispatus = q(d$crispatus), iners = q(d$iners),
             Gardnerella = q(d$Gardnerella), median_score = round(median(d$score), 3))))
write.csv(agree, file.path(out, "tables/frequency_by_subcst.csv"), row.names = FALSE)

eta2_rank <- function(y, g) { r <- rank(y); 1 - sum((r - ave(r, g))^2) / sum((r - mean(r))^2) }
eta2_raw  <- function(y, g) 1 - sum((y - ave(y, g))^2) / sum((y - mean(y))^2)
eta <- function(f) sapply(c("crispatus", "iners", "Gardnerella"), function(v)
  c(subCST = f(sw[[v]], sw$subCST), CST = f(sw[[v]], sw$CST)))
eta_raw  <- eta(eta2_raw)
eta_both <- rbind(round(eta(eta2_rank), 3), round(eta_raw, 3))
rownames(eta_both) <- c("rank_subCST", "rank_CST", "raw_subCST", "raw_CST")
write.csv(eta_both, file.path(out, "tables/eta2_label_vs_frequency.csv"))

etaC <- data.frame(taxon = rep(c("L. crispatus", "L. iners", "Gardnerella"), each = 2),
                   label = rep(c("CST (7)", "subCST (13)"), 3),
                   eta2  = as.vector(eta_raw[c("CST", "subCST"), ]))
etaC$taxon <- factor(etaC$taxon, levels = c("L. crispatus", "L. iners", "Gardnerella"))
figC <- ggplot(etaC, aes(taxon, eta2, fill = label)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", eta2)), position = position_dodge(0.75), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("CST (7)" = "grey55", "subCST (13)" = "grey25")) +
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  labs(x = NULL, y = expression(eta^2~"(share of frequency variation explained by the label)"), fill = NULL,
       title = "Information the label keeps") +
  theme_bw() + theme(axis.text.x = element_text(face = "italic"), legend.position = "bottom")
save_fig(figC, "eta2_information_retained", 5, 5)

# (2) The III/IV boundary: Gardnerella share of (L. iners + Gardnerella), label margin, coexistence zone
sims <- as.matrix(vfull[, grep("_sim$", names(vfull))])
colnames(sims) <- sub("_sim$", "", colnames(sims))
o2 <- t(apply(sims, 1, order, decreasing = TRUE))[, 1:2]
sw$margin <- sims[cbind(seq_len(nrow(sims)), o2[, 1])] - sims[cbind(seq_len(nrow(sims)), o2[, 2])]
sw$second <- colnames(sims)[o2[, 2]]
sw$runner_up_crosses <- grp34(as.character(sw$subCST)) != grp34(sw$second) &
  grp34(as.character(sw$subCST)) %in% c("III", "IV") & grp34(sw$second) %in% c("III", "IV")
sw$pairsum <- sw$iners + sw$Gardnerella
sw$gshare  <- ifelse(sw$pairsum > 0, sw$Gardnerella / sw$pairsum, NA)
sw$zone    <- sw$pairsum > 0.5 & sw$gshare >= 1/11 & sw$gshare <= 10/11
sw$d34     <- apply(sims[, grepl("^III", colnames(sims))], 1, max) - apply(sims[, grepl("^IV", colnames(sims))], 1, max)

pz <- subset(sw, pairsum > 0.5 & subCST %in% c("III-A", "III-B", "IV-A", "IV-B"))
gq <- do.call(rbind, lapply(split(pz$gshare, droplevels(pz$subCST)), function(x)
  round(c(n = length(x), quantile(x, c(.05, .25, .5, .75, .95))), 3)))
write.csv(gq, file.path(out, "tables/gardnerella_share_by_label.csv"))

zones <- list("coexistence zone"         = sw$zone,
              "I-A / III-A outside zone" = !sw$zone & sw$subCST %in% c("I-A", "III-A"),
              "all other swabs"          = !sw$zone & !sw$subCST %in% c("I-A", "III-A"))
conf <- data.frame(where = names(zones), n = sapply(zones, sum),
                   median_margin = round(sapply(zones, function(z) median(sw$margin[z])), 3),
                   pct_margin_below_0.1 = round(100 * sapply(zones, function(z) mean(sw$margin[z] < 0.1)), 1),
                   pct_runner_up_across_III_IV = round(100 * sapply(zones, function(z) mean(sw$runner_up_crosses[z])), 1))
write.csv(conf, file.path(out, "tables/label_confidence_by_zone.csv"), row.names = FALSE)

pz34 <- subset(sw, pairsum > 0.5 & grp34(as.character(subCST)) %in% c("III", "IV"))
fig_d34 <- ggplot(pz34, aes(gshare, d34, fill = subCST)) +
  annotate("rect", xmin = gq["III-B", "95%"], xmax = gq["IV-B", "5%"], ymin = -Inf, ymax = Inf, alpha = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(shape = 21, size = 2, stroke = 0.2, colour = "grey25", alpha = 0.85) +
  facet_grid(. ~ Cohort) +
  scale_fill_manual(values = cst_col, drop = TRUE) +
  labs(x = "Gardnerella share of (L. iners + Gardnerella)", y = "VALENCIA similarity: best CST III − best CST IV",
       fill = "subCST", title = "The CST III/IV label is the sign of a continuous score",
       subtitle = "Swabs where L. iners + Gardnerella > 50% of the community; shaded = G-share band where the label switches") +
  theme_bw()
save_fig(fig_d34, "iii_iv_score_vs_gardnerella_share", 12, 5.5)

# III/IV label flips between consecutive swabs (4-10 days) vs actual change in Gardnerella and community
o <- order(sw$SubjectID, sw$GDcoll)
i <- o[-length(o)]; j <- o[-1]
keep <- sw$SubjectID[i] == sw$SubjectID[j] & (sw$GDcoll[j] - sw$GDcoll[i]) %in% 4:10 &
  grp34(as.character(sw$subCST[i])) %in% c("III", "IV") & grp34(as.character(sw$subCST[j])) %in% c("III", "IV")
i <- i[keep]; j <- j[keep]
pr <- data.frame(SubjectID = sw$SubjectID[i], Cohort = sw$Cohort[i],
                 type = ifelse(grp34(as.character(sw$subCST[i])) == grp34(as.character(sw$subCST[j])),
                               "no flip", "III<->IV flip"),
                 dG = abs(sw$Gardnerella[j] - sw$Gardnerella[i]), d34_from = sw$d34[i], d34_to = sw$d34[j],
                 BC = 0.5 * rowSums(abs(ft[i, , drop = FALSE] - ft[j, , drop = FALSE])))

flips <- do.call(rbind, lapply(split(pr, list(pr$type, pr$Cohort), drop = TRUE), function(d) data.frame(
  Cohort = d$Cohort[1], type = d$type[1], n_pairs = nrow(d), women = length(unique(d$SubjectID)),
  median_abs_dG = round(median(d$dG), 3),
  pct_dG_below_0.05 = round(100 * mean(d$dG < 0.05), 1),
  pct_dG_above_0.10 = round(100 * mean(d$dG > 0.10), 1),
  pct_both_near_boundary = round(100 * mean(abs(d$d34_from) < 0.1 & abs(d$d34_to) < 0.1), 1),
  median_BC = round(median(d$BC), 3),
  pct_BC_below_0.2 = round(100 * mean(d$BC < 0.2), 1),
  pct_BC_above_0.5 = round(100 * mean(d$BC > 0.5), 1))))
write.csv(flips[order(flips$Cohort, flips$type), ], file.path(out, "tables/iii_iv_flips_vs_change.csv"), row.names = FALSE)
saveRDS(sw, file.path(d11, "swabs_cst_and_frequencies.rds"))

# (3) PTB association: continuous frequency vs CST-based exposure, per woman, per cohort
sw$G2 <- rowSums(ft[, taxtab$Gvariant %in% "G2", drop = FALSE])
grpC  <- ifelse(grepl("^IV", sw$CST), "IV", sw$CST)
wom <- unique(sw[, c("SubjectID", "Cohort", "preterm")])
rownames(wom) <- wom$SubjectID
pw  <- function(x) tapply(x, sw$SubjectID, mean)[wom$SubjectID]
for (v in c("crispatus", "iners", "Gardnerella", "G2")) wom[[v]] <- pw(sw[[v]])
wom$frac_I   <- pw(grpC == "I")
wom$frac_III <- pw(grpC == "III")
wom$frac_IVB <- pw(sw$CST == "IV-B")
wom$frac_IV  <- pw(grpC == "IV")
wom$ever_IV  <- wom$frac_IV > 0
wom$modal_IV <- tapply(grpC, sw$SubjectID, function(x) names(which.max(table(x))) == "IV")[wom$SubjectID]
wom$Outcome  <- factor(ifelse(wom$preterm, "Preterm", "Term"), levels = c("Preterm", "Term"))
saveRDS(wom, file.path(d11, "women_cst_and_frequencies.rds"))

wil1 <- function(x, co, alt) wilcox_p(x[wom$Cohort == co & wom$preterm], x[wom$Cohort == co & !wom$preterm], alt)
fis1 <- function(x, co) {
  d <- wom$Cohort == co
  fisher.test(table(factor(x[d], c(FALSE, TRUE)), factor(wom$preterm[d], c(FALSE, TRUE))), alternative = "greater")$p.value
}
eff <- function(x, co, bin) {
  a <- x[wom$Cohort == co & wom$preterm]; b <- x[wom$Cohort == co & !wom$preterm]
  if (bin) sprintf("%d/%d vs %d/%d", sum(a), length(a), sum(b), length(b)) else sprintf("%.3f vs %.3f", median(a), median(b))
}
spec <- data.frame(
  taxon   = c("L. crispatus", "L. crispatus", "L. iners", "L. iners", "Gardnerella", "Gardnerella",
              "Gardnerella / BV-like", "BV-like", "BV-like", "Gardnerella G2"),
  measure = c("crispatus", "frac_I", "iners", "frac_III", "Gardnerella", "frac_IVB", "frac_IV", "ever_IV", "modal_IV", "G2"),
  kind    = c("continuous", "CST", "continuous", "CST", "continuous", "CST", "CST", "CST (yes/no)", "CST (yes/no)", "continuous"),
  alt     = c("less", "less", "less", "less", "greater", "greater", "greater", NA, NA, "greater"))
cmp3 <- do.call(rbind, lapply(seq_len(nrow(spec)), function(k) {
  s <- spec[k, ]; bin <- is.na(s$alt); x <- wom[[s$measure]]
  data.frame(s[, c("taxon", "measure", "kind")],
             Stanford_PTB_vs_term = eff(x, "Stanford", bin),
             p_Stanford = signif(if (bin) fis1(x, "Stanford") else wil1(x, "Stanford", s$alt), 3),
             UAB_PTB_vs_term = eff(x, "UAB", bin),
             p_UAB = signif(if (bin) fis1(x, "UAB") else wil1(x, "UAB", s$alt), 3))
}))
write.csv(cmp3, file.path(out, "tables/ptb_continuous_vs_cst.csv"), row.names = FALSE)

# L. crispatus in women never labelled CST I; Gardnerella in women never labelled CST IV
neverI <- wom$frac_I == 0
mech <- do.call(rbind, lapply(c("Stanford", "UAB"), function(co) {
  d <- wom$Cohort == co & neverI
  data.frame(Cohort = co, women_never_CST_I = sum(d), PTB = sum(d & wom$preterm), term = sum(d & !wom$preterm),
             mean_crispatus_PTB  = signif(mean(wom$crispatus[d & wom$preterm]), 3),
             mean_crispatus_term = signif(mean(wom$crispatus[d & !wom$preterm]), 3),
             term_women_with_crispatus_over_1pct = sum(d & !wom$preterm & wom$crispatus > 0.01),
             PTB_women_with_crispatus_over_1pct  = sum(d & wom$preterm & wom$crispatus > 0.01),
             p_less = signif(wilcox_p(wom$crispatus[d & wom$preterm], wom$crispatus[d & !wom$preterm], "less"), 3))
}))
write.csv(mech, file.path(out, "tables/crispatus_in_never_cst_i_women.csv"), row.names = FALSE)

hid <- do.call(rbind, lapply(c("Stanford", "UAB"), function(co) {
  d <- wom$Cohort == co & wom$frac_IV == 0
  data.frame(Cohort = co, never_IV_PTB = sum(d & wom$preterm), never_IV_term = sum(d & !wom$preterm),
             median_G_PTB = signif(median(wom$Gardnerella[d & wom$preterm]), 2),
             median_G_term = signif(median(wom$Gardnerella[d & !wom$preterm]), 2),
             range_G_never_IV = paste(signif(range(wom$Gardnerella[d]), 2), collapse = " - "),
             p_greater = signif(wilcox_p(wom$Gardnerella[d & wom$preterm], wom$Gardnerella[d & !wom$preterm], "greater"), 3))
}))
write.csv(hid, file.path(out, "tables/gardnerella_in_never_cst_iv_women.csv"), row.names = FALSE)

# Figures: per-woman Gardnerella vs CST IV, and PTB tests by exposure measure
nlab <- with(wom, sapply(c("Stanford", "UAB"), function(co)
  sprintf("%s (n = %d; %d preterm)", co, sum(Cohort == co), sum(Cohort == co & preterm))))
wom$cohort_lab <- factor(nlab[wom$Cohort], levels = nlab)
figB <- ggplot(wom, aes(pmax(Gardnerella, 1e-5), frac_IV, colour = Outcome)) +
  geom_point(size = 2.2, alpha = 0.8, position = position_jitter(width = 0, height = 0.012, seed = 1)) +
  facet_grid(. ~ cohort_lab) +
  scale_x_log10(breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1), labels = c("0", "0.0001", "0.001", "0.01", "0.1", "1")) +
  scale_colour_manual(values = c(Preterm = "#C0392B", Term = "#2E86C1")) +
  labs(x = "Mean Gardnerella relative abundance", y = "Proportion of swabs classified CST IV",
       colour = NULL, title = "Gardnerella exposure per woman: continuous vs CST-based") +
  theme_bw() + theme(legend.position = "bottom")
save_fig(figB, "women_gardnerella_vs_cst_iv", 11, 5.5)

nice <- c(crispatus = "Mean abundance", frac_I = "Proportion CST I", iners = "Mean abundance",
          frac_III = "Proportion CST III", Gardnerella = "Mean abundance", frac_IVB = "Proportion CST IV-B",
          frac_IV = "Proportion CST IV", ever_IV = "Any CST IV", modal_IV = "Majority CST IV", G2 = "Mean abundance ")
grpD <- c(crispatus = "L. crispatus", frac_I = "L. crispatus", iners = "L. iners", frac_III = "L. iners",
          Gardnerella = "Gardnerella", frac_IVB = "Gardnerella", frac_IV = "Gardnerella", ever_IV = "Gardnerella",
          modal_IV = "Gardnerella", G2 = "Gardnerella G2")
cl <- data.frame(measure = rep(cmp3$measure, 2), kind = rep(cmp3$kind, 2),
                 cohort = rep(c("Stanford", "UAB"), each = nrow(cmp3)), p = c(cmp3$p_Stanford, cmp3$p_UAB))
cl$label      <- factor(cl$measure, levels = rev(names(nice)), labels = rev(nice))
cl$cohort_lab <- factor(nlab[cl$cohort], levels = nlab)
cl$method     <- factor(ifelse(cl$kind == "continuous", "Continuous relative abundance", "CST-based (VALENCIA)"),
                        levels = c("Continuous relative abundance", "CST-based (VALENCIA)"))
cl$group      <- factor(grpD[cl$measure], levels = unique(grpD))
figD <- ggplot(cl, aes(-log10(p), label, colour = method, shape = method)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_point(size = 3) +
  facet_grid(group ~ cohort_lab, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c("Continuous relative abundance" = "#C0392B", "CST-based (VALENCIA)" = "grey35")) +
  scale_shape_manual(values = c("Continuous relative abundance" = 16, "CST-based (VALENCIA)" = 17)) +
  labs(x = expression(-log[10]~italic(P)), y = NULL, colour = NULL, shape = NULL,
       title = "Association with preterm birth by exposure measure") +
  theme_bw() + theme(legend.position = "bottom", strip.text.y = element_text(angle = 0, face = "italic"))
save_fig(figD, "ptb_tests_continuous_vs_cst", 10, 7)

# Key-taxon abundance within the CSTs defined by L. crispatus, L. iners and Gardnerella
k5 <- subset(sw, subCST %in% c("I-A", "I-B", "III-A", "III-B", "IV-B"))
k5$subCST <- droplevels(k5$subCST)
long5 <- rbind(data.frame(k5[, c("subCST", "Cohort")], taxon = "L. crispatus", abund = k5$crispatus),
               data.frame(k5[, c("subCST", "Cohort")], taxon = "L. iners", abund = k5$iners),
               data.frame(k5[, c("subCST", "Cohort")], taxon = "Gardnerella", abund = k5$Gardnerella))
long5$taxon <- factor(long5$taxon, levels = c("L. crispatus", "L. iners", "Gardnerella"))
fig_within <- ggplot(long5, aes(subCST, abund, colour = subCST)) +
  geom_jitter(width = 0.25, height = 0, size = 0.7, alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "black", width = 0.5) +
  facet_grid(taxon ~ Cohort) +
  scale_colour_manual(values = c("I-A" = "#1F66D9", "I-B" = "#9DBBEF", "III-A" = "#F5A800",
                                 "III-B" = "#F8D58A", "IV-B" = "#D91A26"), guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "VALENCIA subCST", y = "Relative abundance in the swab", title = "Key-taxon abundance within each CST label") +
  theme_bw() + theme(strip.text.y = element_text(face = "italic"))
save_fig(fig_within, "abundance_within_cst_labels", 11, 7)
write.csv(aggregate(abund ~ taxon + subCST + Cohort, data = long5,
                    FUN = function(x) round(quantile(x, c(.05, .25, .5, .75, .95)), 3)),
          file.path(out, "tables/abundance_quantiles_within_cst.csv"), row.names = FALSE)
