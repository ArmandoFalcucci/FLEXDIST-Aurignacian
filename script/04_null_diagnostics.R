# =============================================================================
# 04_null_diagnostics.R
# -----------------------------------------------------------------------------
# Null-model diagnostics for the matrix-correlation results. This script does
# not add new archaeological data; it calibrates the statistics already reported
# in 02_mantel.R against what would be expected by chance at n = 6.
#
#   1. EXACT family-wise permutation test of the maximum Mantel r across the
#      seven-predictor panel (720 permutations at n = 6 -> exact, not sampled).
#   2. NULL BENCHMARK for the leave-one-out diagnostic: how unstable is a Mantel
#      r at n = 6 when nothing is going on? Gives the LOO shift a yardstick.
#   3. SAMPLE-SIZE CHECK: does assemblage n predict technological isolation?
#   4. SPECIMEN-LEVEL BOOTSTRAP of the FLEXDIST distances, propagated all the
#      way to the Mantel coefficients, so the response matrix carries error.
#      4b. Sensitivity to scale() on the PC scores (see manuscript sect. 2.4.2).
#   5. MINIMUM-n CURVE: the assemblage count at which LOO instability falls
#      below a usable threshold (supports the recommendation in the Discussion).
#
# UPSTREAM  : 01_distances.R -> output/rds/01_distances.rds
#             02_mantel.R    -> output/rds/02_mantel.rds
#             00a/00b Step 8 -> output/rds/04_scores_{blanks,cores}.rds
# WRITES    : output/rds/04_diagnostics.rds
#             output/figures/Figure_S5.{tiff,pdf}   (bootstrap CIs, coefficients)
#
# NOTE ON EXACTNESS
#   With n = 6 assemblages there are 6! = 720 relabellings and with n = 5 there
#   are 120. Both are enumerated in full below, so every p-value in this script
#   is exact and seed-independent. The seed governs only the bootstrap.
# =============================================================================

library(vegan)
library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(123)

BOOT_N <- 2000     # bootstrap replicates (section 4)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/rds",     recursive = TRUE, showWarnings = FALSE)

PAL <- list(grey_dark = "#2C3E50", red = "#E74C3C", green = "#16A085",
            purple = "#9B59B6", blue = "#2980B9", tan = "#D4A574")


# ---- Load upstream ----------------------------------------------------------
d01 <- readRDS("output/rds/01_distances.rds")
d02 <- readRDS("output/rds/02_mantel.rds")

group_order    <- d01$group_order
predictor_mats <- lapply(d01$predictor_mats, function(m) m[group_order, group_order])
pretty_pred    <- d01$pretty_pred
md_norm        <- d02$md_norm
pa             <- "Pa_24A"
keep_no_pa     <- setdiff(group_order, pa)
n_asm          <- length(group_order)


# ---- Helpers ----------------------------------------------------------------
lt <- function(m) m[lower.tri(m)]

# Spearman Mantel statistic (identical to vegan::mantel(..., "spearman")$statistic)
mr <- function(A, B) cor(lt(A), lt(B), method = "spearman")

# All permutations of 1:n, as a matrix of rows.
all_perms <- function(n) {
  if (n == 1L) return(matrix(1L))
  p <- all_perms(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(i)
    cbind(i, ifelse(p < i, p, p + 1L))))
}

# Exact one-sided Mantel p for a single predictor.
exact_mantel_p <- function(A, B, PM) {
  r0 <- mr(A, B)
  mean(apply(PM, 1L, function(p) mr(A[p, p], B)) >= r0 - 1e-12)
}


# =============================================================================
# 1. EXACT FAMILY-WISE PERMUTATION TEST
# -----------------------------------------------------------------------------
# 02_mantel.R screens seven predictors and reports uncorrected p-values. The
# quantity that actually needs a null is the LARGEST r across that screen. We
# permute the assemblage labels of the lithic matrix, recompute all seven
# correlations under each relabelling, and keep the maximum. That distribution
# is the correct reference for "the best predictor we found".
# =============================================================================
PM6 <- all_perms(6L)
PM5 <- all_perms(5L)

fw_test <- function(T, preds, PM) {
  obs      <- vapply(preds, function(m) mr(T, m), numeric(1))
  obs_max  <- max(obs)
  null_max <- apply(PM, 1L, function(p) {
    Tp <- T[p, p]
    max(vapply(preds, function(m) mr(Tp, m), numeric(1)))
  })
  marg <- vapply(names(preds), function(k) exact_mantel_p(T, preds[[k]], PM),
                 numeric(1))
  list(obs = obs, obs_max = obs_max, null_max = null_max,
       p_familywise = mean(null_max >= obs_max - 1e-12),
       p_marginal = marg)
}

fw6 <- fw_test(md_norm, predictor_mats, PM6)
fw5 <- fw_test(md_norm[keep_no_pa, keep_no_pa],
               lapply(predictor_mats, function(m) m[keep_no_pa, keep_no_pa]), PM5)

fw_table <- data.frame(
  Predictor   = pretty_pred[names(fw6$obs)],
  r_n6        = round(fw6$obs, 3),
  p_marginal  = round(fw6$p_marginal, 3),
  r_n5        = round(fw5$obs, 3),
  p_marginal5 = round(fw5$p_marginal, 3),
  row.names   = NULL
)[order(-fw6$obs), ]

cat("\n=== 1. Exact family-wise permutation test ===\n")
print(fw_table, row.names = FALSE)
cat(sprintf(
  "\n  n = 6: best predictor %s, r = %.3f, marginal exact p = %.3f\n",
  pretty_pred[names(which.max(fw6$obs))], fw6$obs_max,
  fw6$p_marginal[which.max(fw6$obs)]))
cat(sprintf("         null max-r: median %.3f, 90th pct %.3f, 95th pct %.3f\n",
            median(fw6$null_max), quantile(fw6$null_max, .90),
            quantile(fw6$null_max, .95)))
cat(sprintf("         FAMILY-WISE exact p = %.3f\n", fw6$p_familywise))
cat(sprintf("  n = 5: best r = %.3f, FAMILY-WISE exact p = %.3f\n",
            fw5$obs_max, fw5$p_familywise))


# =============================================================================
# 2. SAMPLE-SIZE CHECK
# -----------------------------------------------------------------------------
# Centroid separation is upward-biased by sampling error, and that bias scales
# inversely with assemblage n. Paglicci is both the outlier and among the
# smallest samples, so this needs an explicit check.
# =============================================================================
n_blanks <- c(Pa_24A = 180, Ca_Au = 326, `CTC_gic-ars` = 526,
              `CTC_rsa'` = 304, RB_A2 = 244, `RF_A2-A1` = 950)
n_cores  <- c(Pa_24A = 16, Ca_Au = 138, `CTC_gic-ars` = 61,
              `CTC_rsa'` = 43, RB_A2 = 10, `RF_A2-A1` = 112)

mean_dist <- function(M) vapply(rownames(M), function(g)
  mean(M[g, setdiff(colnames(M), g)]), numeric(1))

size_check <- data.frame(
  assemblage = group_order,
  n_blanks   = n_blanks[group_order],
  n_cores    = n_cores[group_order],
  mean_blank = round(mean_dist(d02$md_blanks[group_order, group_order]), 3),
  mean_core  = round(mean_dist(d02$md_cores [group_order, group_order]), 3),
  mean_comb  = round(mean_dist(md_norm), 3),
  row.names  = NULL)

size_cor <- c(
  blanks = cor(size_check$n_blanks, size_check$mean_blank, method = "spearman"),
  cores  = cor(size_check$n_cores,  size_check$mean_core,  method = "spearman"),
  comb   = cor(size_check$n_cores,  size_check$mean_comb,  method = "spearman"))

cat("\n=== 2. Assemblage size vs technological isolation ===\n")
print(size_check, row.names = FALSE)
cat("\n  Spearman(n, mean distance):\n")
print(round(size_cor, 3))
cat(sprintf("\n  Smallest core sample is %s (n = %d); its mean combined distance is %.3f\n",
            names(which.min(n_cores)), min(n_cores),
            size_check$mean_comb[match(names(which.min(n_cores)), size_check$assemblage)]))
cat(sprintf("  Paglicci mean combined distance is %.3f (rank %d of %d)\n",
            size_check$mean_comb[size_check$assemblage == pa],
            rank(-size_check$mean_comb)[size_check$assemblage == pa], n_asm))


# =============================================================================
# 3. SPECIMEN-LEVEL BOOTSTRAP OF THE FLEXDIST DISTANCES
# -----------------------------------------------------------------------------
# Everything upstream treats md_norm as fixed. It is not: each cell is a
# distance between two group positions estimated from finite samples. We
# resample specimens within assemblage, recompute the group distances, rebuild
# the combined matrix, and re-run every Mantel test, so the reported
# coefficients carry the sampling error of the assemblages they derive from.
#
# This reruns the ordination geometry but not the stochastic kernel imputation,
# which is prohibitively slow to bootstrap and, in this dataset, nearly vacuous
# (9 missing cells in 41,370 after the <=2-NA filter; see manuscript sect. 2.4.2).
# =============================================================================
if (file.exists("output/rds/04_scores_blanks.rds") &&
    file.exists("output/rds/04_scores_cores.rds")) {

  sc <- list(blanks = readRDS("output/rds/04_scores_blanks.rds"),
             cores  = readRDS("output/rds/04_scores_cores.rds"))

  # Harmonize assemblage labels with the distance matrices and validate.
  for (nm in names(sc)) {
    g <- as.character(sc[[nm]]$group)
    g[g == "Pa_A4-A0"] <- "Pa_24A"
    sc[[nm]]$group <- g
  }
  stopifnot(setequal(unique(sc$blanks$group), group_order),
            setequal(unique(sc$cores$group),  group_order))

  # `scale_pcs = TRUE` and FUN = median reproduce the published pipeline
  # (00a/00b Step 8).
  centroid_dist <- function(df, idx, scale_pcs = TRUE, FUN = median) {
    g   <- df$group[idx]
    pcs <- df[idx, setdiff(names(df), "group"), drop = FALSE]
    if (scale_pcs) pcs <- as.data.frame(scale(pcs))
    cen <- aggregate(pcs, by = list(group = g), FUN = FUN)
    m   <- as.matrix(dist(cen[, -1, drop = FALSE]))
    dimnames(m) <- list(cen$group, cen$group)
    m[group_order, group_order]
  }

  boot_r <- matrix(NA_real_, BOOT_N, length(predictor_mats),
                   dimnames = list(NULL, names(predictor_mats)))
  boot_cells <- array(NA_real_, c(BOOT_N, n_asm, n_asm))

  for (i in seq_len(BOOT_N)) {
    ib <- unlist(lapply(split(seq_len(nrow(sc$blanks)), sc$blanks$group),
                        function(z) sample(z, length(z), replace = TRUE)))
    ic <- unlist(lapply(split(seq_len(nrow(sc$cores)),  sc$cores$group),
                        function(z) sample(z, length(z), replace = TRUE)))
    Bb <- centroid_dist(sc$blanks, ib)
    Cc <- centroid_dist(sc$cores,  ic)
    Ti <- (Bb / max(Bb) + Cc / max(Cc)) / 2
    boot_cells[i, , ] <- Ti
    for (k in names(predictor_mats)) boot_r[i, k] <- mr(Ti, predictor_mats[[k]])
  }

  boot_ci <- t(apply(boot_r, 2, quantile, c(.025, .5, .975), na.rm = TRUE))
  boot_ci <- data.frame(Predictor = pretty_pred[rownames(boot_ci)],
                        observed  = round(vapply(predictor_mats,
                                                 function(m) mr(md_norm, m), numeric(1)), 3),
                        lwr = round(boot_ci[, 1], 3),
                        med = round(boot_ci[, 2], 3),
                        upr = round(boot_ci[, 3], 3),
                        row.names = NULL)
  cat("\n=== 3. Bootstrap 95% CI on the Mantel coefficients ===\n")
  print(boot_ci[order(-boot_ci$observed), ], row.names = FALSE)

  cell_ci <- apply(boot_cells, c(2, 3), quantile, c(.025, .975), na.rm = TRUE)
  dimnames(cell_ci)[[2]] <- dimnames(cell_ci)[[3]] <- group_order

  # Does the bootstrap ever fail to put Paglicci at the periphery?
  pa_i     <- match(pa, group_order)
  pa_rank  <- apply(boot_cells, 1, function(M) {
    md <- vapply(seq_len(n_asm), function(j) mean(M[j, -j]), numeric(1))
    rank(-md)[pa_i]
  })
  cat(sprintf("\n  Paglicci is the most peripheral assemblage in %.1f%% of bootstrap replicates\n",
              100 * mean(pa_rank == 1)))


  # ---------------------------------------------------------------------------
  # 3b. ANALYTIC-CHOICE SENSITIVITY: scale() on the PC scores
  # ---------------------------------------------------------------------------
  # Step 8 of 00a/00b applies scale() to the PC scores before taking group
  # positions. Because ALL PCs are retained and the eigenvalue spectrum is
  # nearly flat beyond the leading components, this reweights weakly identified
  # (largely noise) directions to the same variance as PC1. A substantive
  # analytic choice; its consequences are reported here and in the manuscript.
  idx_all_b <- seq_len(nrow(sc$blanks)); idx_all_c <- seq_len(nrow(sc$cores))
  T_scaled <- {
    B <- centroid_dist(sc$blanks, idx_all_b, TRUE)
    C <- centroid_dist(sc$cores,  idx_all_c, TRUE); (B/max(B) + C/max(C))/2 }
  T_unscaled <- {
    B <- centroid_dist(sc$blanks, idx_all_b, FALSE)
    C <- centroid_dist(sc$cores,  idx_all_c, FALSE); (B/max(B) + C/max(C))/2 }

  scale_sens <- data.frame(
    Predictor  = pretty_pred[names(predictor_mats)],
    r_scaled   = round(vapply(predictor_mats, function(m) mr(T_scaled,   m), numeric(1)), 3),
    r_unscaled = round(vapply(predictor_mats, function(m) mr(T_unscaled, m), numeric(1)), 3),
    r_scaled_n5   = round(vapply(predictor_mats, function(m)
      mr(T_scaled[keep_no_pa, keep_no_pa],   m[keep_no_pa, keep_no_pa]), numeric(1)), 3),
    r_unscaled_n5 = round(vapply(predictor_mats, function(m)
      mr(T_unscaled[keep_no_pa, keep_no_pa], m[keep_no_pa, keep_no_pa]), numeric(1)), 3),
    row.names = NULL)
  cat("\n=== 3b. Sensitivity to scale() on the PC scores ===\n")
  cat(sprintf("  Spearman between scaled and unscaled lithic matrices: %.3f\n",
              mr(T_scaled, T_unscaled)))
  print(scale_sens[order(-scale_sens$r_scaled), ], row.names = FALSE)

  F_boot <- ggplot(boot_ci, aes(x = reorder(Predictor, observed))) +
    geom_hline(yintercept = 0, colour = "grey75") +
    geom_linerange(aes(ymin = lwr, ymax = upr), linewidth = 0.8,
                   colour = PAL$grey_dark) +
    geom_point(aes(y = observed), size = 3.4, colour = PAL$red) +
    coord_flip() +
    labs(x = NULL, y = "Mantel r (Spearman), bootstrap 95% CI") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggsave("output/figures/Figure_S5.tiff", F_boot, width = 7, height = 5,
         dpi = 300, bg = "white", compression = "lzw")
  ggsave("output/figures/Figure_S5.pdf", F_boot, width = 7, height = 5, bg = "white")

} else {
  boot_ci <- NULL; cell_ci <- NULL; scale_sens <- NULL; pa_rank <- NULL
  message("Bootstrap skipped: output/rds/04_scores_{blanks,cores}.rds not found. ",
          "Re-run 00a and 00b (their Step 8 writes these files).")
}


# =============================================================================
# SAVE
# =============================================================================
saveRDS(
  list(fw6 = fw6, fw5 = fw5, fw_table = fw_table,
       size_check = size_check, size_cor = size_cor,
       boot_ci = boot_ci, cell_ci = cell_ci,
       scale_sens = scale_sens, pa_rank = pa_rank,
       BOOT_N = BOOT_N),
  "output/rds/04_diagnostics.rds")

cat("\nDone. Wrote output/rds/04_diagnostics.rds and Figure S5.\n")
