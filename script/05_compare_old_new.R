# =============================================================================
# 05_compare_old_new.R
# -----------------------------------------------------------------------------
# Old vs new results after correcting the Castelcivita coordinates and the
# Paglicci 24A age. Base R only; changes nothing in the analysis outputs.
#
# BEFORE rerunning anything, keep a copy of the current caches:
#   dir.create("output/rds_old", showWarnings = FALSE)
#   file.copy(list.files("output/rds", full.names = TRUE), "output/rds_old")
# Then edit 01 (CTC coordinates, Pa_24A date), run 01 -> 02 -> 04, and source
# this file from the project root.
#
# WRITES : output/check/compare_old_new.txt  (everything printed below)
#          output/check/compare_*.csv        (the same tables as CSV)
# =============================================================================

options(width = 200)

old_dir <- "output/rds_old"
new_dir <- "output/rds"
out_dir <- "output/check"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rd <- function(dir, f) readRDS(file.path(dir, f))
o1 <- rd(old_dir, "01_distances.rds"); n1 <- rd(new_dir, "01_distances.rds")
o2 <- rd(old_dir, "02_mantel.rds");    n2 <- rd(new_dir, "02_mantel.rds")
has04 <- all(file.exists(file.path(c(old_dir, new_dir), "04_diagnostics.rds")))
if (has04) {
  o4 <- rd(old_dir, "04_diagnostics.rds"); n4 <- rd(new_dir, "04_diagnostics.rds")
}

g    <- n1$group_order
pp   <- n1$pretty_pred
lt   <- function(m) m[lower.tri(m)]
hdr  <- function(x) cat("\n\n==== ", x, " ====\n", sep = "")
save <- function(df, name) write.csv(df, file.path(out_dir, paste0("compare_", name, ".csv")),
                                     row.names = FALSE)
flip <- function(p_old, p_new) (p_old < 0.05) != (p_new < 0.05)

sink(file.path(out_dir, "compare_old_new.txt"), split = TRUE)
cat("Old:", normalizePath(old_dir), "\nNew:", normalizePath(new_dir), "\n")


# ---- 1. Inputs --------------------------------------------------------------
hdr("1. Inputs (coordinates, median ages, climate slice in ka)")
inp <- data.frame(
  assemblage = g,
  lat_old   = o1$coordinates[g, "latitude"],  lat_new = n1$coordinates[g, "latitude"],
  lon_old   = o1$coordinates[g, "longitude"], lon_new = n1$coordinates[g, "longitude"],
  age_old   = unname(o1$site_dates[g]),       age_new = unname(n1$site_dates[g]),
  slice_old = round(unname(o1$site_dates[g]) / 1000),
  slice_new = round(unname(n1$site_dates[g]) / 1000),
  row.names = NULL)
print(inp, row.names = FALSE)
save(inp, "inputs")


# ---- 2. Climate values (Table 5) -------------------------------------------
hdr("2. Climate values per assemblage (Table 5): rows that changed")
clim_cmp <- function(o, n, vars, tag) do.call(rbind, lapply(vars, function(v)
  data.frame(assemblage = g, variable = paste0(v, tag),
             old = round(o[[v]], 2), new = round(n[[v]], 2),
             change = round(n[[v]] - o[[v]], 2), row.names = NULL)))
clim <- rbind(
  clim_cmp(o1$clim_native,     n1$clim_native,     c("bio01", "bio12", "npp"), ""),
  clim_cmp(o1$clim_downscaled, n1$clim_downscaled, c("bio01", "bio12"),        "_ds"))
clim_changed <- clim[clim$change != 0, ]
if (nrow(clim_changed)) print(clim_changed, row.names = FALSE) else cat("No change.\n")
cat("bio01 = MAT (C), bio12 = TAP (mm), npp = NPP; _ds = downscaled 10 arc-min\n")
save(clim, "climate")


# ---- 3. Predictor matrices --------------------------------------------------
hdr("3. Predictor matrices: how much did each change?")
pm_cmp <- do.call(rbind, lapply(names(n1$predictor_mats), function(p) {
  a <- lt(o1$predictor_mats[[p]][g, g]); b <- lt(n1$predictor_mats[[p]][g, g])
  data.frame(predictor         = unname(pp[p]),
             pairs_changed     = sum(abs(a - b) > 1e-9),
             spearman_old_new  = round(suppressWarnings(cor(a, b, method = "spearman")), 3),
             same_rank_order   = identical(rank(a), rank(b)),
             max_abs_change    = signif(max(abs(a - b)), 3),
             row.names = NULL)
}))
print(pm_cmp, row.names = FALSE)
cat("same_rank_order = TRUE means Spearman Mantel results for that predictor cannot change.\n")
save(pm_cmp, "predictor_matrices")

hdr("3b. Least-cost distances from Castelcivita (old -> new)")
to  <- setdiff(g, c("CTC_gic-ars", "CTC_rsa'"))
lcp <- data.frame(to = to,
                  old = round(o1$costpath_dist["CTC_gic-ars", to]),
                  new = round(n1$costpath_dist["CTC_gic-ars", to]),
                  row.names = NULL)
lcp$change_pct <- round(100 * (lcp$new / lcp$old - 1), 1)
print(lcp, row.names = FALSE)
save(lcp, "costpath_ctc")


# ---- 4-5. Simple Mantel -----------------------------------------------------
mantel_cmp <- function(o, n, label, file) {
  m <- merge(o, n, by = "predictor", suffixes = c("_old", "_new"))
  m <- data.frame(predictor      = unname(pp[m$predictor]),
                  r_sp_old       = round(m$r_spearman_old, 3),
                  r_sp_new       = round(m$r_spearman_new, 3),
                  d_r_sp         = round(m$r_spearman_new - m$r_spearman_old, 3),
                  p_sp_old       = round(m$p_spearman_old, 3),
                  p_sp_new       = round(m$p_spearman_new, 3),
                  r_pe_old       = round(m$r_pearson_old, 3),
                  r_pe_new       = round(m$r_pearson_new, 3),
                  p_pe_old       = round(m$p_pearson_old, 3),
                  p_pe_new       = round(m$p_pearson_new, 3),
                  crosses_0.05   = flip(m$p_spearman_old, m$p_spearman_new) |
                                   flip(m$p_pearson_old,  m$p_pearson_new))
  m <- m[order(-m$r_sp_new), ]
  hdr(label)
  print(m, row.names = FALSE)
  cat("Order old:", paste(pp[o$predictor[order(-o$r_spearman)]], collapse = " > "), "\n")
  cat("Order new:", paste(pp[n$predictor[order(-n$r_spearman)]], collapse = " > "), "\n")
  save(m, file)
  m
}
m6 <- mantel_cmp(o2$mantel_n6, n2$mantel_n6, "4. Simple Mantel, n = 6 (Table 6, Table S2)", "mantel_n6")
m5 <- mantel_cmp(o2$mantel_n5, n2$mantel_n5, "5. Simple Mantel, n = 5 (Table 6, Table S2)", "mantel_n5")


# ---- 6-7. Partial Mantel ----------------------------------------------------
hdr("6. Partial Mantel, n = 5 (main text)")
pt <- merge(o2$partials, n2$partials, by = "test", suffixes = c("_old", "_new"))
pt$d_r <- round(pt$r_new - pt$r_old, 3)
pt$crosses_0.05 <- flip(pt$p_old, pt$p_new)
print(pt, row.names = FALSE)
save(pt, "partials")

hdr("7. Partial Mantel sweep, n = 5 (Table S3), sorted by size of change")
sp <- merge(o2$supp_partials, n2$supp_partials, by = c("control", "predictor"),
            suffixes = c("_old", "_new"))
sp$d_r <- round(sp$r_partial_new - sp$r_partial_old, 3)
sp$crosses_0.05 <- flip(sp$p_partial_old, sp$p_partial_new)
sp <- sp[order(-abs(sp$d_r)), ]
print(sp, row.names = FALSE)
save(sp, "partials_sweep")


# ---- 8. Leave-one-out (Figure 7a) -------------------------------------------
hdr("8. Leave-one-out r (Figure 7a): new minus old")
d_loo <- n2$loo_r[rownames(o2$loo_r), colnames(o2$loo_r)] - o2$loo_r
rownames(d_loo) <- pp[rownames(d_loo)]
print(round(d_loo, 3))
save(data.frame(predictor = rownames(d_loo), round(d_loo, 3), check.names = FALSE), "loo_change")


# ---- 9. Lithic side: must be identical --------------------------------------
hdr("9. Lithic side (should be identical)")
cat("Combined lithic matrix identical:", isTRUE(all.equal(o2$md_norm, n2$md_norm)), "\n")
cat("PCA validation Spearman:", round(o2$validation[["spearman"]], 3), "->",
    round(n2$validation[["spearman"]], 3), "\n")
cat("Blank-core Mantel r:    ", round(o2$agreement$mantel_bc$statistic, 3), "->",
    round(n2$agreement$mantel_bc$statistic, 3), "\n")


# ---- 10-12. Null diagnostics (04) -------------------------------------------
if (has04) {
  hdr("10. Exact family-wise test (Table S4 caption)")
  fw <- data.frame(
    n         = c(6, 5),
    best_old  = unname(pp[c(names(which.max(o4$fw6$obs)), names(which.max(o4$fw5$obs)))]),
    best_new  = unname(pp[c(names(which.max(n4$fw6$obs)), names(which.max(n4$fw5$obs)))]),
    max_r_old = round(c(o4$fw6$obs_max, o4$fw5$obs_max), 3),
    max_r_new = round(c(n4$fw6$obs_max, n4$fw5$obs_max), 3),
    p_fw_old  = round(c(o4$fw6$p_familywise, o4$fw5$p_familywise), 3),
    p_fw_new  = round(c(n4$fw6$p_familywise, n4$fw5$p_familywise), 3))
  print(fw, row.names = FALSE)
  save(fw, "familywise")

  hdr("11. Exact marginal p (Table S4)")
  ft <- merge(o4$fw_table, n4$fw_table, by = "Predictor", suffixes = c("_old", "_new"))
  print(ft, row.names = FALSE)
  save(ft, "exact_marginal")

  if (!is.null(o4$boot_ci) && !is.null(n4$boot_ci)) {
    hdr("12. Bootstrap 95% CI on Mantel r (Figure S5)")
    bc <- merge(o4$boot_ci, n4$boot_ci, by = "Predictor", suffixes = c("_old", "_new"))
    bc <- bc[order(-bc$observed_new), ]
    print(bc, row.names = FALSE)
    save(bc, "bootstrap_ci")
    cat("\nPa_24A most peripheral (% of replicates):",
        round(100 * mean(o4$pa_rank == 1), 1), "->",
        round(100 * mean(n4$pa_rank == 1), 1), "(should be identical)\n")
  }
}


# ---- Summary ----------------------------------------------------------------
hdr("SUMMARY")
cat(sprintf("Largest change in Spearman r, n = 6: %+.3f (%s)\n",
            m6$d_r_sp[which.max(abs(m6$d_r_sp))], m6$predictor[which.max(abs(m6$d_r_sp))]))
cat(sprintf("Largest change in Spearman r, n = 5: %+.3f (%s)\n",
            m5$d_r_sp[which.max(abs(m5$d_r_sp))], m5$predictor[which.max(abs(m5$d_r_sp))]))
cat("Predictor ranking changed (n = 6):",
    !identical(o2$mantel_n6$predictor[order(-o2$mantel_n6$r_spearman)],
               n2$mantel_n6$predictor[order(-n2$mantel_n6$r_spearman)]), "\n")
cat("Predictor ranking changed (n = 5):",
    !identical(o2$mantel_n5$predictor[order(-o2$mantel_n5$r_spearman)],
               n2$mantel_n5$predictor[order(-n2$mantel_n5$r_spearman)]), "\n")
cat("Tests crossing p = 0.05: simple n = 6:", sum(m6$crosses_0.05),
    "| simple n = 5:", sum(m5$crosses_0.05),
    "| partials:", sum(pt$crosses_0.05), "| sweep:", sum(sp$crosses_0.05), "\n")
if (has04) cat("Family-wise p, n = 6:", fw$p_fw_old[1], "->", fw$p_fw_new[1],
               "| n = 5:", fw$p_fw_old[2], "->", fw$p_fw_new[2], "\n")

sink()
cat("\nWrote", file.path(out_dir, "compare_old_new.txt"), "and compare_*.csv\n")
