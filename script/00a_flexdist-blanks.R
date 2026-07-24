# =============================================================================
# R workflow for analyzing lithic flake variation with FLEXDIST
# -----------------------------------------------------------------------------
# Author:  H. Rathmann (updated by A. Falcucci; cleaned for submission)
# Version date:  September 2026
#
# This workflow uses functions from FLEXDIST, introduced in:
#   Rathmann et al. 2023. Estimating inter-individual Mahalanobis distances
#   from mixed incomplete high-dimensional data: Application to human skeletal
#   remains from 3rd to 1st millennia BC Southwest Germany.
#   Journal of Archaeological Science 156: 105802.
#   https://doi.org/10.1016/j.jas.2023.105802
#
# DATA FORMAT REQUIREMENTS
#   - A .csv matrix with individuals in rows and variables in columns.
#   - Column 1: unique individual IDs (e.g., Indiv_1, Indiv_2, ...).
#   - Column 2: grouping variable (e.g., Group_1, Group_2, ...).
#   - Remaining columns: nominal, ordinal, or continuous variables (any mix).
#       * Header "N" = nominal   -> non-negative integers differing by 1 (0/1, 1/2)
#       * Header "O" = ordinal   -> ordered non-negative integers (0,1,2,3; 1..6)
#       * Header "C" = continuous-> measurements on any scale
#   - Missing values must be coded as NA.
#
# REPRODUCIBILITY NOTES
#   - Parts of this workflow are stochastic (kernel-based imputation, jitter,
#     sina plots). A global seed is set below so results are reproducible.
#   - All paths are relative to the project root. Run from the project folder,
#     or open the accompanying .Rproj. (here::here() is a good alternative.)
# =============================================================================


### Step 1: Load packages -----------------------------------------------------
library(ade4)        # dudi.mix (mixed-data PCA), scalewt
library(ggplot2)     # plotting
library(dplyr)       # data wrangling / pipe
library(RColorBrewer)# palettes
library(ggrepel)     # non-overlapping plot labels
library(vegan)       # betadisper (group dispersions)
library(ggforce)     # sina plots, geom_circle
library(psych)       # mixedCor (correlations among mixed variables)
library(ggcorrplot)  # correlation heatmap
library(patchwork)   # combine plots into one figure
library(grid)        # arrow(), unit()

# Global seed for the stochastic steps (imputation, jitter, sina).
set.seed(123)


### Step 2: Load data and set parameters --------------------------------------

# Path to the data file, relative to the project root (edit if your file differs).
file_path <- "data/dataset_blanks_FLEXDIST.csv"
dat1 <- read.csv(file = file_path, header = FALSE)

# Remove the "site", "layer", and "blank" variables (original columns 2 and 4).
# The remaining "site+layer" column is used as the grouping variable.
dat1 <- dat1[, -c(2, 4)]

# Remove individuals with more than `max_na` missing values across the
# variable columns (columns 3 onward). Default = 2.
max_na <- 2
var_cols <- 3:ncol(dat1)
dat1 <- dat1[rowSums(is.na(dat1[var_cols])) <= max_na, ]

# Variable-type header row (N / O / C) and counts.
ivar  <- dat1[1, -c(1:2)]        # variable types
nvar  <- length(ivar)            # number of variables
ID    <- dat1[-1, 1]             # individual IDs
nind  <- length(ID)              # number of individuals
site  <- dat1[-1, 2]             # grouping variable (site + layer)
nsite <- length(unique(site))    # number of groups

# Relabel one assemblage for clearer figure legends.
# Display only: `site` is used solely for grouping/plotting (p1, p3),
# not in the imputation or distance estimation, so this does not affect results.
site[site == "Pa_A4-A0"] <- "Pa_24A"

# Extract the data block and coerce to numeric.
dat2 <- dat1[-1, -c(1:2)]
dat2[1:nvar] <- lapply(dat2[1:nvar], as.numeric)

# Sanity checks.
str(dat2)
summary(dat2)


### Step 3: Generate multiple differentially imputed datasets -----------------
# Number of resampling runs (increase for final publication-quality estimates).
runs <- 150

# Two imputed versions:
#   cent.imp = single "central" dataset (NAs -> median of kernel draws)
#   mult.imp = `runs` stacked replicates (NAs -> individual kernel draws)
cent.imp <- dat2
mult.imp <- do.call("rbind", replicate(runs, dat2, simplify = FALSE))

for (i in 1:nvar) {
  d1 <- cent.imp[, i]
  d2 <- mult.imp[, i]
  na <- sum(is.na(d2))
  mi <- min(d1, na.rm = TRUE)
  ma <- max(d1, na.rm = TRUE)
  dens <- density(d1, n = 1e5, from = mi, to = ma, na.rm = TRUE)
  kern.samp <- sample(dens$x, na, replace = TRUE, prob = dens$y)
  if (ivar[i] == "N") {
    d1[is.na(d1)] <- median(round(kern.samp, 0))
    d2[is.na(d2)] <- round(kern.samp, 0)
    cent.imp[, i] <- as.factor(d1)
    mult.imp[, i] <- as.factor(d2)
  } else if (ivar[i] == "O") {
    d1[is.na(d1)] <- median(round(kern.samp, 0))
    d2[is.na(d2)] <- round(kern.samp, 0)
    cent.imp[, i] <- as.ordered(d1)
    mult.imp[, i] <- as.ordered(d2)
  } else if (ivar[i] == "C") {
    d1[is.na(d1)] <- median(kern.samp)
    d2[is.na(d2)] <- kern.samp
    cent.imp[, i] <- as.numeric(d1)
    mult.imp[, i] <- as.numeric(d2)
  }
}

str(cent.imp)
summary(cent.imp, maxsum = 10)


### Step 4 (optional): size-correct (Q-mode) the continuous data --------------
# Applied here to the shape-related continuous variables:
#   1 = length, 2 = width, 3 = thickness, 4 = platform_width, 5 = platform_thickness
# To disable size correction, comment out this whole block.
cvar <- c(1, 2, 3, 4, 5)
cvar.cent <- cent.imp[, cvar]
cvar.mult <- mult.imp[, cvar]

# Geometric mean per individual, then divide each variable by it.
geo.mean.cent <- apply(cvar.cent, 1, function(x) exp(mean(log(x))))
geo.mean.mult <- apply(cvar.mult, 1, function(x) exp(mean(log(x))))
qmode.cent <- cvar.cent / geo.mean.cent
qmode.mult <- cvar.mult / geo.mean.mult

cent.imp[, cvar] <- qmode.cent
mult.imp[, cvar] <- qmode.mult

str(cent.imp)
summary(cent.imp, maxsum = 10)


### Step 5: Estimate Mahalanobis distances via mixed-data PCA -----------------
# Iterated over the imputed datasets. Axes retained to reach ~90% cumulative
# variance; PC scores scaled; Euclidean distance on scaled scores.
sta <- 1
sto <- nind
dist.mult <- array(0, dim = c(nind, nind, runs))
for (i in 1:runs) {
  temp.dat <- mult.imp[sta:sto, ]
  pca <- dudi.mix(temp.dat, scannf = FALSE, nf = 1e7)
  eigval  <- pca$eig
  propvar <- (eigval / sum(eigval)) * 100
  cumvar  <- cumsum(propvar)
  index   <- sum(cumvar < 90) + 1
  pcs <- pca$li
  if (length(pca$cw) >= nind) {
    pcs <- pcs[, 1:index]
  }
  pcs <- scale(pcs)
  dist.mult[, , i] <- as.matrix(dist(pcs))
  sta <- sta + nind
  sto <- sto + nind
}


### Step 6: Summary statistic across the imputed distances --------------------
# Median distance matrix (used downstream). The 2.5% / 97.5% percentiles below
# give a credible interval if needed for reporting.
FLEXDIST_median <- apply(dist.mult, 1:2, quantile, 0.5, names = FALSE)
# FLEXDIST_lower <- apply(dist.mult, 1:2, quantile, 0.025, names = FALSE)
# FLEXDIST_upper <- apply(dist.mult, 1:2, quantile, 0.975, names = FALSE)

# Save all key objects so the plotting steps can be re-run without recomputing.
save(
  dat1, dat2, ivar, nvar, ID, nind, site, nsite,
  cent.imp, mult.imp, qmode.cent, qmode.mult,
  dist.mult, FLEXDIST_median,
  file = "FLEXDIST_results.150runs.RData"
)


### Step 7: Plots -------------------------------------------------------------

# Short variable labels (order matches the data columns).
lab <- c("LEN", "WID", "THI", "PLW", "PLT", "PLTY", "BULB", "LIP",
         "CX", "SCARP", "CS", "CURV", "TWIST", "BLSH", "DIST")

# To skip Steps 2-6 and jump straight to plotting, uncomment the next line:
# load("FLEXDIST_results.150runs.RData")

# Directory for figure output.
if (!dir.exists("output")) dir.create("output")


## Step 7.1: scatterplot of principal components ------------------------------

# Mixed-data PCA on the central imputed dataset.
colnames(cent.imp) <- lab
pca     <- dudi.mix(cent.imp, scannf = FALSE, nf = 1e7)
eigval  <- pca$eig
propvar <- (eigval / sum(eigval)) * 100

# Project all imputed replicates onto the PCA space of the central dataset.
Xori <- cent.imp
Xsup <- mult.imp
res  <- matrix(0, nrow(Xsup), 1)
for (i in 1:ncol(Xsup)) {
  if (pca$index[i] == "q") {
    scale.Xori <- scalewt(df = Xori[, i], wt = pca$lw, center = TRUE, scale = TRUE)
    center  <- attr(scale.Xori, "scaled:center")
    norm    <- attr(scale.Xori, "scaled:scale")
    var.tmp <- scale(Xsup[, i], scale = norm, center = center)
    res <- cbind(res, var.tmp)
  } else if (pca$index[i] == "f") {
    fac <- Xsup[, i]
    n   <- length(fac)
    mat <- matrix(0, n, nlevels(fac))
    mat[(1:n) + n * (unclass(fac) - 1)] <- 1
    var.tmp <- mat
    col.w   <- pca$cw[pca$assign == i]
    var.tmp <- t(t(var.tmp) / col.w) - 1
    res <- cbind(res, var.tmp)
  } else if (pca$index[i] == "o") {
    tmp.ori <- as.numeric(Xori[, i])
    n       <- length(tmp.ori)
    polys   <- poly(tmp.ori, 2)
    tmp.sup <- as.numeric(Xsup[, i])
    var.tmp <- sqrt(n - 1) * predict(polys, tmp.sup)
    res <- cbind(res, var.tmp)
  }
}
res2 <- res[, -1]
Xsup.coo <- res2 %*% (as.matrix(pca$c1) * pca$cw)
colnames(Xsup.coo) <- names(pca$li)

# Combine the first two PC scores with ID and site labels.
plot.1 <- as.data.frame(Xsup.coo[, 1:2])
colnames(plot.1) <- c("V1", "V2")
plot.1$ID   <- rep(ID, runs)
plot.1$site <- rep(site, runs)

# Small jitter so ellipses can be drawn even for complete individuals.
set.seed(123)
plot.1$V1 <- jitter(plot.1$V1, amount = 0.0001)
plot.1$V2 <- jitter(plot.1$V2, amount = 0.0001)

# Preserve original ID order (prevents alphabetical reordering in aggregate()).
plot.1$ID <- factor(plot.1$ID, levels = unique(plot.1$ID))

# Individual centroids, then site centroids.
plot.2 <- aggregate(cbind(V1, V2) ~ ID, data = plot.1, FUN = median)
plot.2$group <- site
plot.3 <- aggregate(cbind(V1, V2) ~ site, data = plot.2, FUN = median)

p1 <- ggplot() +
  stat_ellipse(data = plot.1, aes(x = V1, y = V2, fill = site, group = ID),
               geom = "polygon", alpha = 0.05, level = 0.95) +
  geom_point(data = plot.2, aes(x = V1, y = V2, color = site, group = ID), size = 2) +
  stat_ellipse(data = plot.2, aes(x = V1, y = V2, group = group),
               color = "white", linewidth = 3) +
  stat_ellipse(data = plot.2, aes(x = V1, y = V2, color = group), linewidth = 1.5) +
  geom_point(data = plot.3, aes(x = V1, y = V2), color = "white", shape = 15, size = 5) +
  geom_point(data = plot.3, aes(x = V1, y = V2, color = site), shape = 15, size = 3.5) +
  xlab(paste0("PC1 - ", round(propvar[1], 2), "%")) +
  ylab(paste0("PC2 - ", round(propvar[2], 2), "%")) +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line   = element_line(colour = "black"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x  = element_text(size = 13),
        axis.text.y  = element_text(size = 13),
        legend.position = "bottom")

p1
# ggsave("output/p1.tiff", p1, device = "tiff", width = 10, height = 8, units = "in", dpi = 300)


## Step 7.2: PCA correlation circle -------------------------------------------

data.circle <- as.data.frame(pca$co[, 1:2])
colnames(data.circle) <- c("PC1", "PC2")

set.seed(123)
p2 <- ggplot(data.circle, aes(PC1, PC2)) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.03, "npc")),
               color = "steelblue", linewidth = 1, alpha = 0.8) +
  geom_circle(aes(x0 = 0, y0 = 0, r = 1),
              color = "darkgrey", linetype = "dashed", linewidth = 0.8) +
  geom_label_repel(aes(label = rownames(data.circle)),
                   box.padding = 0.3, point.padding = 0.2,
                   segment.color = "grey50", segment.size = 0.5,
                   size = 4, fill = "white", color = "black") +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line  = element_line(colour = "black", linewidth = 0.8),
        axis.title = element_text(size = 16),
        axis.text  = element_text(size = 14),
        legend.position = "none",
        plot.margin = margin(15, 15, 15, 15)) +
  labs(x = paste0("PC1 (", round(propvar[1], 2), "%)"),
       y = paste0("PC2 (", round(propvar[2], 2), "%)"))

p2
# ggsave("output/p2.tiff", p2, width = 8, height = 6, dpi = 300)


## Step 7.3: boxplot (sina) of group dispersions ------------------------------

mod <- betadisper(as.dist(FLEXDIST_median), site)
data.box <- data.frame(distance = mod$distance, Assemblage = mod$group)

# Median with 95% interval, for the summary overlay.
my_ci <- function(x) data.frame(
  y    = quantile(x, 0.5),
  ymin = quantile(x, 0.025),
  ymax = quantile(x, 0.975)
)

set.seed(123)
p3 <- ggplot(data.box, aes(x = Assemblage, y = distance, color = Assemblage)) +
  geom_sina(shape = 19, size = 2) +
  stat_summary(fun.data = "my_ci", colour = "white", linewidth = 3) +
  stat_summary(fun.data = "my_ci", colour = "black", linewidth = 1.5) +
  scale_y_continuous(name = "Distance to group centroid") +
  scale_color_brewer(palette = "Paired") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line   = element_line(colour = "black"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x  = element_text(size = 13, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 13),
        legend.position = "none")

p3
# ggsave("output/p3.tiff", p3, device = "tiff", width = 10, height = 8, units = "in", dpi = 300)


## Step 7.4: correlations among mixed variables -------------------------------

dat4 <- cent.imp
dat4[1:nvar] <- lapply(dat4[1:nvar], as.numeric)  # mixedCor requires numeric coding
r <- mixedCor(data = dat4, smooth = FALSE)
cormat <- r$rho

# Order variables from distal to length for the heatmap.
desired_order <- c("DIST", "BLSH", "TWIST", "CURV", "CS", "SCARP", "CX",
                   "LIP", "BULB", "PLTY", "PLT", "PLW", "THI", "WID", "LEN")
cormat <- cormat[desired_order, desired_order]

col1 <- colorRampPalette(brewer.pal(9, "BrBG"))(100)

p4 <- ggcorrplot(cormat,
                 hc.order = FALSE,
                 type = "lower",
                 lab = FALSE,
                 colors = c(col1[1], "white", col1[100]),
                 outline.col = "gray",
                 legend.title = "Correlation",
                 ggtheme = ggplot2::theme_void() +
                   theme(legend.position = c(0.15, 0.78),
                         legend.background = element_rect(fill = "white", color = "white"),
                         axis.text.x = element_text(angle = 45, hjust = 1)))

p4


## Step 7.5: combine plots for the paper --------------------------------------

combined_plot <- (p4 + p2) / (p1) +
  plot_annotation(tag_levels = "a")

tiff("output/combined_figure.tif", width = 12, height = 10, units = "in", res = 300)
print(combined_plot)
dev.off()

# Vector (PDF) — scalable, preferred for print/production
pdf("output/combined_figure.pdf", width = 12, height = 10)
print(combined_plot)
dev.off()


### Step 8: Distances among group centroids -----------------------------------
# PCA on the central imputed dataset, aggregate PC scores to group medians,
# then Euclidean distances among group centroids. Written out as a CSV.
pca     <- dudi.mix(cent.imp, scannf = FALSE, nf = 1e7)
eigval  <- pca$eig
propvar <- (eigval / sum(eigval)) * 100
cumvar  <- cumsum(propvar)
index   <- sum(cumvar < 90) + 1

pcs <- pca$li
if (length(pca$cw) >= nind) {
  pcs <- pcs[, 1:index]
}
pcs <- as.data.frame(scale(pcs))

tmp <- aggregate(x = pcs, by = list(site), FUN = median)
md  <- as.matrix(dist(tmp[, -1]))
rownames(md) <- tmp[, 1]
colnames(md) <- tmp[, 1]
md

write.csv(md, file = "output/distance_matrix_blanks.csv", row.names = TRUE)