# =============================================================================
# 02_mantel.R
# -----------------------------------------------------------------------------
# Core statistical analysis:
#   1. Combine the FLEXDIST blank + core distance matrices into one lithic distance (range-normalised mean) and check blank/core agreement.
#   2. Simple Mantel tests, full (n = 6) and Paglicci-excluded (n = 5).
#   3. Leave-one-out leverage across all predictors.
#   4. Partial Mantel tests at n = 5.
#   -> Figure 6 (tanglegram + heatmap + PCoA), Figure 7 (leverage + dumbbell),
#      and the Mantel/partial tables.
#
# UPSTREAM  : 00_flexdist.R  -> data/distance_matrix_*.csv
#             01_distances.R -> output/rds/01_distances.rds
# WRITES    : output/rds/02_mantel.rds  (results the Rmd needs)
#             output/figures/Figure_6.*, Figure_7.*, F_*.tiff
#
# Permutations: at n = 6 (720) and n = 5 (120) vegan enumerates the COMPLETE
# permutation set, so results are exact and seed-independent; permutations =
# 9999 simply requests "all of them". The seed only affects sampled steps.
# =============================================================================

# ---- Packages ---------------------------------------------------------------
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(dendextend)
library(viridis)
library(flextable)

set.seed(123)

# ---- Load upstream objects --------------------------------------------------
d01 <- readRDS("output/rds/01_distances.rds")
group_order    <- d01$group_order
predictor_mats <- d01$predictor_mats
pretty_pred    <- d01$pretty_pred

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/rds",     recursive = TRUE, showWarnings = FALSE)

# ---- Shared aesthetics ------------------------------------------------------
CLUST_COLS <- c("#1b9e77", "#d95f02", "#7570b3")     # k = 3 cluster colours
PAL <- list(grey_dark = "#2C3E50", red = "#E74C3C", green = "#16A085",
            purple = "#9B59B6", blue = "#2980B9", tan = "#D4A574")

pa <- "Pa_24A"                                        # Paglicci label

make_ft <- function(df) {
  flextable(df) |>
    theme_vanilla() |>
    flextable::align(j = 2:ncol(df), align = "center", part = "all") |>
    set_table_properties(layout = "autofit") |>
    autofit()
}

# p-value cannot resolve below 1/(nperm + 1); helper for reporting.
fmt_p <- function(p, nperm) {
  floor_p <- 1 / (nperm + 1)
  if (p <= floor_p) sprintf("< %.4f", floor_p) else sprintf("= %.4f", p)
}


# =============================================================================
# 1. COMBINE BLANK + CORE DISTANCES  (n = 6)
# =============================================================================
md_blanks <- as.matrix(read.csv("data/distance_matrix_blanks.csv",
                                row.names = 1, check.names = FALSE))
md_cores  <- as.matrix(read.csv("data/distance_matrix_cores.csv",
                                row.names = 1, check.names = FALSE))

stopifnot(setequal(group_order, rownames(md_blanks)),
          setequal(group_order, rownames(md_cores)))
md_blanks <- md_blanks[group_order, group_order]
md_cores  <- md_cores [group_order, group_order]

# Range-normalise each, then average -> equal weight to blanks and cores.
md_norm <- (md_blanks / max(md_blanks) + md_cores / max(md_cores)) / 2


# =============================================================================
# 2. CLUSTERING + BLANK/CORE AGREEMENT
# =============================================================================
hc       <- hclust(as.dist(md_norm), method = "ward.D2")
ord      <- hc$labels[hc$order]
clust_k3 <- cutree(hc, k = 3)

hc_blanks <- hclust(as.dist(md_blanks), method = "ward.D2")
hc_cores  <- hclust(as.dist(md_cores),  method = "ward.D2")

# Do blanks and cores split into the same three clusters?
cat("\n--- k = 3 cluster cross-tabulation (blanks vs cores) ---\n")
print(table(blanks = cutree(hc_blanks, 3)[group_order],
            cores  = cutree(hc_cores,  3)[group_order]))

dend_blanks <- as.dendrogram(hc_blanks) |>
  color_branches(k = 3, col = CLUST_COLS) |> set("branches_lwd", 2)
dend_cores  <- as.dendrogram(hc_cores) |>
  color_branches(k = 3, col = CLUST_COLS) |> set("branches_lwd", 2)
dl <- dendlist(dend_blanks, dend_cores) |> untangle(method = "step2side")

# Agreement metrics (reported in the Results text)
ent       <- round(entanglement(dl), 3)
cor_coph  <- round(cor_cophenetic(dl[[1]], dl[[2]]), 3)
cor_bg    <- round(cor_bakers_gamma(dl[[1]], dl[[2]]), 3)
mantel_bc <- mantel(as.dist(md_blanks), as.dist(md_cores),
                    method = "spearman", permutations = 9999)

cat("\n--- Blank vs core agreement ---\n")
cat(sprintf("Mantel r (Spearman):    %.3f (p %s)\n",
            mantel_bc$statistic, fmt_p(mantel_bc$signif, 9999)))
cat(sprintf("Cophenetic correlation: %.3f\n", cor_coph))
cat(sprintf("Baker's gamma:          %.3f\n", cor_bg))
cat(sprintf("Entanglement:           %.3f\n", ent))


# =============================================================================
# 3. FIGURE 6  (a tanglegram, b heatmap, c PCoA)
# =============================================================================
plot_half_dend <- function(dend, side = c("left", "right"), title = "") {
  side <- match.arg(side)
  gg   <- as.ggdend(dend); segs <- gg$segments; labs <- gg$labels
  segs$col[is.na(segs$col)] <- "grey45"
  if (side == "left") {
    p <- ggplot() +
      geom_segment(data = segs, aes(x = -y, xend = -yend, y = x, yend = xend),
                   colour = segs$col, linewidth = 0.8, lineend = "round") +
      geom_text(data = labs, aes(x = 0, y = x, label = label),
                hjust = 0, nudge_x = 0.015, size = 3.7, colour = "grey15") +
      scale_x_continuous(expand = expansion(mult = c(0.05, 0.30)))
  } else {
    p <- ggplot() +
      geom_segment(data = segs, aes(x = y, xend = yend, y = x, yend = xend),
                   colour = segs$col, linewidth = 0.8, lineend = "round") +
      geom_text(data = labs, aes(x = 0, y = x, label = label),
                hjust = 1, nudge_x = -0.015, size = 3.7, colour = "grey15") +
      scale_x_continuous(expand = expansion(mult = c(0.30, 0.05)))
  }
  p + scale_y_continuous(expand = expansion(mult = 0.04)) +
    labs(title = title) + theme_void(base_size = 11) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold", size = 12,
                                     margin = margin(b = 6)),
          plot.margin = margin(6, 4, 6, 4))
}

p_tangle <- (plot_half_dend(dl[[1]], "left",  "Blanks") +
             plot_half_dend(dl[[2]], "right", "Cores")) +
  plot_layout(widths = c(1, 1)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))
p_tangle <- wrap_elements(p_tangle) + labs(tag = "a") +
  theme(plot.tag = element_text(face = "bold", size = 16),
        plot.tag.position = c(0.01, 0.98))

heatmap_df <- expand.grid(site1 = ord, site2 = ord)
heatmap_df$dist <- md_norm[cbind(as.character(heatmap_df$site1),
                                 as.character(heatmap_df$site2))]
p_heatmap <- ggplot(heatmap_df, aes(site1, site2, fill = dist)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", dist)),
            colour = ifelse(heatmap_df$dist > 0.6, "white", "grey20"), size = 3) +
  scale_fill_viridis_c(option = "mako", direction = -1,
                       name = "Distance", limits = c(0, 1)) +
  scale_y_discrete(limits = rev) + coord_equal() +
  labs(x = NULL, y = NULL, tag = "b") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid = element_blank(),
        plot.tag = element_text(face = "bold", size = 16),
        plot.tag.position = c(0.01, 0.98))

pcoa    <- cmdscale(as.dist(md_norm), k = 2, eig = TRUE)
pcoa_df <- data.frame(site = rownames(pcoa$points),
                      PCo1 = pcoa$points[, 1], PCo2 = pcoa$points[, 2],
                      cluster = factor(clust_k3[rownames(pcoa$points)]))
ve1 <- round(pcoa$eig[1] / sum(abs(pcoa$eig)) * 100, 1)
ve2 <- round(pcoa$eig[2] / sum(abs(pcoa$eig)) * 100, 1)
p_pcoa <- ggplot(pcoa_df, aes(PCo1, PCo2, label = site)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(aes(colour = cluster), size = 4.8) +
  geom_text_repel(size = 3.8, box.padding = 0.5, segment.colour = "grey60") +
  scale_colour_manual(values = setNames(CLUST_COLS, levels(pcoa_df$cluster)),
                      name = "Cluster") +
  labs(tag = "c", x = paste0("PCo1 (", ve1, "%)"), y = paste0("PCo2 (", ve2, "%)")) +
  theme_minimal(base_size = 11) +
  theme(plot.tag = element_text(face = "bold", size = 16),
        plot.tag.position = c(0.01, 0.98), legend.position = "right")

Figure_6 <- p_tangle / (p_heatmap + p_pcoa + plot_layout(widths = c(1, 1.4))) +
  plot_layout(heights = c(0.9, 1.1))
ggsave("output/figures/Figure_6.tiff", Figure_6, width = 10.5, height = 12,
       dpi = 300, bg = "white", compression = "lzw")
ggsave("output/figures/Figure_6.pdf",  Figure_6, width = 10.5, height = 12,
       bg = "white")


# =============================================================================
# 4. PCA VALIDATION: subset of n = 6 vs independent rerun on n = 5
# =============================================================================
keep_no_pa <- setdiff(group_order, pa)

md.blanks.noPA <- as.matrix(read.csv("data/distance_matrix_blanks_noPA24A.csv",
                                     row.names = 1, check.names = FALSE))[keep_no_pa, keep_no_pa]
md.cores.noPA  <- as.matrix(read.csv("data/distance_matrix_cores_noPA24A.csv",
                                     row.names = 1, check.names = FALSE))[keep_no_pa, keep_no_pa]
md_norm_noPA   <- (md.blanks.noPA / max(md.blanks.noPA) +
                   md.cores.noPA  / max(md.cores.noPA)) / 2
md_norm_subset <- md_norm[keep_no_pa, keep_no_pa]

lt <- function(m) m[lower.tri(m)]
validation <- c(
  pearson  = cor(lt(md_norm_subset), lt(md_norm_noPA)),
  spearman = cor(lt(md_norm_subset), lt(md_norm_noPA), method = "spearman"))
cat("\n=== PCA validation: subset of n=6 vs rerun on n=5 ===\n")
cat(sprintf("Pearson r:  %.3f\nSpearman r: %.3f\n",
            validation["pearson"], validation["spearman"]))


# =============================================================================
# 5. SIMPLE MANTEL: full (n = 6) and Paglicci-excluded (n = 5)
# =============================================================================
# Alignment guard
stopifnot(all(vapply(predictor_mats,
                     function(m) setequal(rownames(m), group_order), logical(1))))
predictor_mats <- lapply(predictor_mats, function(m) m[group_order, group_order])

simple_mantel <- function(keep) {
  tech <- as.dist(md_norm[keep, keep])
  do.call(rbind, lapply(names(predictor_mats), function(nm) {
    pd <- as.dist(predictor_mats[[nm]][keep, keep])
    sp <- mantel(tech, pd, method = "spearman", permutations = 9999)
    pe <- mantel(tech, pd, method = "pearson",  permutations = 9999)
    data.frame(predictor = nm,
               r_spearman = sp$statistic, p_spearman = sp$signif,
               r_pearson  = pe$statistic, p_pearson  = pe$signif)
  }))
}

mantel_n6 <- simple_mantel(group_order) |> arrange(desc(r_spearman))
mantel_n5 <- simple_mantel(keep_no_pa)  |> arrange(desc(r_spearman))
cat("\n=== Simple Mantel (n = 6) ===\n"); print(mantel_n6, row.names = FALSE, digits = 3)
cat("\n=== Simple Mantel (n = 5) ===\n"); print(mantel_n5, row.names = FALSE, digits = 3)

# Table 6: full vs Paglicci-excluded, Spearman, with the shift.
Table_6 <- mantel_n6 |>
  transmute(predictor, r6 = r_spearman, p6 = p_spearman) |>
  left_join(transmute(mantel_n5, predictor, r5 = r_spearman, p5 = p_spearman),
            by = "predictor") |>
  mutate(shift = r5 - r6) |> arrange(desc(r6)) |>
  transmute(Predictor = pretty_pred[predictor],
            `r (n=6)` = sprintf("%.3f", r6), `p (n=6)` = sprintf("%.3f", p6),
            `r (n=5)` = sprintf("%.3f", r5), `p (n=5)` = sprintf("%.3f", p5),
            `dr (n5 - n6)` = sprintf("%+.3f", shift)) |>
  make_ft() |> flextable::bold(part = "header")


# =============================================================================
# 6. LEAVE-ONE-OUT LEVERAGE  (computed once, reused)
# =============================================================================
loo_r <- sapply(group_order, function(drop) {
  keep <- setdiff(group_order, drop)
  tech <- as.dist(md_norm[keep, keep])
  vapply(predictor_mats, function(m)
    mantel(tech, as.dist(m[keep, keep]),
           method = "spearman", permutations = 9999)$statistic, numeric(1))
})
r_full <- setNames(mantel_n6$r_spearman, mantel_n6$predictor)[rownames(loo_r)]

# Largest r change from dropping any NON-Paglicci assemblage
other_sites   <- setdiff(colnames(loo_r), pa)
other_changes <- abs(loo_r[, other_sites] - r_full)
worst_per_pred <- data.frame(
  predictor  = rownames(loo_r),
  max_change = round(apply(other_changes, 1, max), 3),
  worst_site = other_sites[apply(other_changes, 1, which.max)])
cat("\n=== Worst non-Paglicci leverage per predictor ===\n")
print(worst_per_pred, row.names = FALSE)
cat(sprintf("Largest r change dropping any non-Paglicci assemblage: %.3f\n",
            max(other_changes)))
cat("\n=== Change when dropping Paglicci (for contrast) ===\n")
print(data.frame(predictor = rownames(loo_r),
                 pa_change = round(abs(loo_r[, pa] - r_full), 3)), row.names = FALSE)


# =============================================================================
# 7. FIGURE 7  (a leverage heatmap, b dumbbell n=6 -> n=5)
# =============================================================================
loo_heatmap <- data.frame(
  predictor = rep(rownames(loo_r), ncol(loo_r) + 1),
  condition = rep(c(colnames(loo_r), "(full n = 6)"), each = nrow(loo_r)),
  r         = c(as.vector(loo_r), r_full))
loo_heatmap$condition <- factor(loo_heatmap$condition,
                                levels = c(group_order, "(full n = 6)"))
loo_heatmap$predictor <- factor(pretty_pred[loo_heatmap$predictor],
                                levels = rev(pretty_pred[rownames(loo_r)]))

F_leverage <- ggplot(loo_heatmap, aes(condition, predictor, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r)),
            colour = ifelse(abs(loo_heatmap$r) > 0.5, "white", "grey15"), size = 3.2) +
  geom_vline(xintercept = length(group_order) + 0.5, linewidth = 0.6, colour = "grey20") +
  scale_fill_gradient2(low = PAL$purple, mid = "white", high = PAL$green,
                       midpoint = 0, name = "Mantel r",
                       limits = c(-0.4, 0.85), oob = scales::squish) +
  labs(x = "Assemblage dropped", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())

dumbbell <- mantel_n6 |>
  transmute(predictor, r_n6 = r_spearman) |>
  left_join(transmute(mantel_n5, predictor, r_n5 = r_spearman), by = "predictor") |>
  mutate(shift_positive = r_n5 > r_n6) |> arrange(r_n6)
dumbbell$predictor <- factor(pretty_pred[dumbbell$predictor],
                             levels = pretty_pred[dumbbell$predictor])

F_dumbbell <- ggplot(dumbbell) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_segment(aes(x = r_n6, xend = r_n5, y = predictor, yend = predictor,
                   colour = shift_positive), linewidth = 0.7,
               arrow = arrow(length = unit(0.09, "inches"), type = "closed")) +
  geom_point(aes(r_n6, predictor, fill = "Full dataset (n = 6)"),
             shape = 21, size = 4.2, colour = "white") +
  geom_point(aes(r_n5, predictor, fill = "Paglicci excluded (n = 5)"),
             shape = 21, size = 4.2, colour = "white", alpha = 0.6) +
  scale_fill_manual(values = c("Full dataset (n = 6)" = PAL$grey_dark,
                               "Paglicci excluded (n = 5)" = PAL$red), name = NULL) +
  scale_colour_manual(values = c("TRUE" = PAL$green, "FALSE" = "grey55"), guide = "none") +
  scale_x_continuous(name = "Mantel r (Spearman)", limits = c(-0.3, 0.85),
                     breaks = seq(-0.2, 0.8, 0.2)) +
  labs(y = NULL) + theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "grey95"))

Figure_7 <- (F_leverage / F_dumbbell) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14))
ggsave("output/figures/Figure_7.tiff", Figure_7, width = 10, height = 11,
       dpi = 300, bg = "white", compression = "lzw")
ggsave("output/figures/Figure_7.pdf",  Figure_7, width = 10, height = 11, bg = "white")


# =============================================================================
# 8. SCATTER: combined lithic distance vs each predictor (n = 6)  [Figure S3]
# =============================================================================
mat_to_pairs <- function(mat, value_name) {
  lt2 <- which(lower.tri(mat), arr.ind = TRUE)
  out <- data.frame(site1 = rownames(mat)[lt2[, 1]],
                    site2 = rownames(mat)[lt2[, 2]],
                    val = mat[lt2], stringsAsFactors = FALSE)
  out$pair <- paste0(pmin(out$site1, out$site2), "_", pmax(out$site1, out$site2))
  out$involves_paglicci <- out$site1 == pa | out$site2 == pa
  names(out)[names(out) == "val"] <- value_name
  out[, c("pair", value_name, "involves_paglicci")]
}
predictors_to_plot <- c("TAP", "TAP_ds", "MAT", "MAT_ds", "chrono", "costpath")

pred_long <- do.call(rbind, lapply(predictors_to_plot, function(p) {
  df <- mat_to_pairs(predictor_mats[[p]], "predictor_value"); df$predictor <- p
  df[, c("pair", "predictor", "predictor_value")]
}))
resp_long <- mat_to_pairs(md_norm, "lithic")
all_pairs <- pred_long |>
  left_join(resp_long[, c("pair", "lithic", "involves_paglicci")], by = "pair")
all_pairs$predictor <- factor(all_pairs$predictor, levels = predictors_to_plot)

mantel_grid <- mantel_n6 |>
  filter(predictor %in% predictors_to_plot) |>
  transmute(predictor = factor(predictor, levels = predictors_to_plot),
            label = sprintf("r = %.2f, p = %.3f", r_spearman, p_spearman))

F_scatter <- ggplot(all_pairs, aes(predictor_value, lithic)) +
  geom_smooth(method = "lm", se = TRUE, colour = PAL$grey_dark,
              fill = "grey88", alpha = 0.4, linewidth = 0.6) +
  geom_point(aes(colour = involves_paglicci), size = 2.6, alpha = 0.85) +
  geom_text(data = mantel_grid, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.08, vjust = 1.4, size = 3.3, colour = "grey30",
            inherit.aes = FALSE) +
  facet_wrap(~ predictor, scales = "free_x", ncol = 3,
             labeller = labeller(predictor = pretty_pred)) +
  scale_colour_manual(values = c("FALSE" = PAL$grey_dark, "TRUE" = PAL$red),
                      labels = c("Other pairs", "Includes Paglicci"), name = NULL) +
  labs(x = "Predictor distance", y = "Combined lithic distance") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom", panel.grid.minor = element_blank())
ggsave("output/figures/F_combined_vs_predictors.tiff", F_scatter,
       width = 10, height = 7, dpi = 300, bg = "white", compression = "lzw")


# =============================================================================
# 9. PARTIAL MANTELS at n = 5
# =============================================================================
tech5 <- as.dist(md_norm[keep_no_pa, keep_no_pa])
d5    <- function(nm) as.dist(predictor_mats[[nm]][keep_no_pa, keep_no_pa])

pm <- function(x, z) {
  m <- mantel.partial(tech5, d5(x), d5(z), method = "spearman", permutations = 9999)
  data.frame(test = sprintf("%s | %s", pretty_pred[x], pretty_pred[z]),
             r = round(m$statistic, 3), p = round(m$signif, 3))
}
partials <- rbind(pm("TAP", "costpath"), pm("costpath", "TAP"),
                  pm("TAP", "chrono"),   pm("chrono", "TAP"))
cat("\n=== Partial Mantels (n = 5) ===\n"); print(partials, row.names = FALSE)

# Supplementary: full sweep (control one predictor, test each of the rest).
partial_sweep <- function(control_name) {
  others <- setdiff(names(predictor_mats), control_name)
  do.call(rbind, lapply(others, function(p) {
    m <- mantel.partial(tech5, d5(p), d5(control_name),
                        method = "spearman", permutations = 9999)
    data.frame(control = pretty_pred[control_name], predictor = pretty_pred[p],
               r_partial = round(m$statistic, 3), p_partial = round(m$signif, 3))
  }))
}
supp_partials <- rbind(partial_sweep("costpath"),
                       partial_sweep("chrono"),
                       partial_sweep("TAP"))
cat("\n=== Supplementary partial Mantel sweep (n = 5) ===\n")
print(supp_partials, row.names = FALSE)


# =============================================================================
# SAVE for the supplementary Rmd
# =============================================================================
saveRDS(
  list(
    md_norm = md_norm, md_blanks = md_blanks, md_cores = md_cores,
    mantel_n6 = mantel_n6, mantel_n5 = mantel_n5,
    Table_6 = Table_6, partials = partials, supp_partials = supp_partials,
    loo_r = loo_r, validation = validation,
    agreement = list(entanglement = ent, cophenetic = cor_coph,
                     bakers_gamma = cor_bg, mantel_bc = mantel_bc),
    pa = pa, keep_no_pa = keep_no_pa,
    predictor_mats = predictor_mats, pretty_pred = pretty_pred,
    PAL = PAL, make_ft = make_ft
  ),
  "output/rds/02_mantel.rds"
)
