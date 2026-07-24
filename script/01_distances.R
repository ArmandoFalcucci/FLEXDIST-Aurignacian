# =============================================================================
# 01_distances.R
# -----------------------------------------------------------------------------
# Builds every PREDICTOR distance matrix used in the Mantel analyses:
#   - chronological distance (from Bayesian median ages)
#   - environmental distances (native + downscaled Krapp2021, via pastclim)
#   - least-cost-path distance (Tobler hiking function over ETOPO)
#
# It also (optionally) regenerates the downscaled palaeoclimate model. That step
# is SLOW and needs large downloads, so it is OFF by default: the repo ships the
# precomputed result (data/Europe_bioclims_10min.nc) and this script simply
# loads it. Set RUN_DOWNSCALING <- TRUE only if you want to rebuild it.
#
# UPSTREAM  : 00_flexdist.R presents the steps done to obtain the distance matrices for both blanks and cores. 
# WRITES    : output/rds/01_distances.rds   (all objects 02 and the Rmd need)
# RUN ORDER : 00 -> 01 -> 02 -> 03 ; Rmd last.
# =============================================================================

# ---- Packages ---------------------------------------------------------------
library(pastclim)
library(terra)
library(gdistance)
library(sp)
library(raster)
library(geosphere)
library(ggplot2)
library(patchwork)


# ---- Config -----------------------------------------------------------------
RUN_DOWNSCALING <- FALSE   # TRUE = rebuild the 10-min model from scratch (slow)


if (!file.exists(model_path)) {
  stop("Downscaled model not found at '", model_path,
       "'. Place Europe_bioclims_10min.nc in data/ (from the Zenodo release), ",
       "or set RUN_DOWNSCALING <- TRUE to rebuild it.")
}

# Precomputed downscaled model shipped with the repo:
model_path   <- "data/Europe_bioclims_10min.nc"
n_europe_ext <- c(-10, 25, 35, 70)          # xmin, xmax, ymin, ymax
time_steps   <- seq(-50000, -30000, 1000)   # slices used downstream

dir.create("output/rds", recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# pastclim setup (REQUIRED before any climate call)
# -----------------------------------------------------------------------------

# Optional: point pastclim at a persistent data folder (edit for your machine)
# set_data_path("~/pastclim_data")

# set_data_path(on_CRAN = FALSE)      # persistent user data dir
# # set_data_path("~/pastclim_data")  # ...or a folder of your choice
# 
# # Krapp2021 is needed on EVERY run (for native_data + extractions), not only
# # when RUN_DOWNSCALING = TRUE. Safe to re-run: already-present data is skipped.
# download_dataset(dataset = "Krapp2021",
#                  bio_variables = c("bio01", "bio04", "bio12", "bio15", "npp"))


# =============================================================================
# (OPTIONAL) DOWNSCALING — rebuild data/Europe_bioclims_10min.nc
# =============================================================================
# Delta-change downscaling of Krapp2021 to 10-min using WorldClim as the
# high-resolution reference. Only runs if RUN_DOWNSCALING = TRUE.
if (RUN_DOWNSCALING) {

  tavg_vars <- c(paste0("temperature_0", 1:9),   paste0("temperature_", 10:12))
  prec_vars <- c(paste0("precipitation_0", 1:9), paste0("precipitation_", 10:12))

  download_dataset("Krapp2021",          bio_variables = c(tavg_vars, prec_vars))
  download_dataset("WorldClim_2.1_10m",  bio_variables = c(tavg_vars, prec_vars))

  # Low-resolution model series (Krapp2021)
  tavg_series <- region_series(tavg_vars, time_bp = time_steps,
                               dataset = "Krapp2021", ext = n_europe_ext)
  prec_series <- region_series(prec_vars, time_bp = time_steps,
                               dataset = "Krapp2021", ext = n_europe_ext)

  # High-resolution observations (WorldClim, present day)
  tavg_obs <- terra::crop(region_series(tavg_vars, time_ce = 1985,
                          dataset = "WorldClim_2.1_10m", ext = n_europe_ext),
                          n_europe_ext)
  prec_obs <- region_series(prec_vars, time_ce = 1985,
                            dataset = "WorldClim_2.1_10m", ext = n_europe_ext)
  tavg_range <- range(unlist(lapply(tavg_obs, minmax, compute = TRUE)))
  prec_range <- range(unlist(lapply(prec_obs, minmax, compute = TRUE)))

  # High-resolution land mask (ETOPO relief), minus ice and internal seas
  download_etopo()
  relief_rast <- terra::resample(load_etopo(), tavg_obs[[1]])
  land_high   <- make_land_mask(relief_rast = relief_rast, time_bp = 0)

  ice_low  <- get_ice_mask(time_bp = time_steps, dataset = "Krapp2021")
  ice_high <- downscale_ice_mask(ice_mask_low_res = ice_low,
                                 land_mask_high_res = land_high)
  land_high <- mask(land_high, ice_high, inverse = TRUE)

  internal_seas <- readRDS(system.file("extdata/internal_seas.RDS",
                                       package = "pastclim"))
  crs(internal_seas) <- crs(land_high)
  land_high <- mask(land_high, internal_seas, inverse = TRUE)

  downscale_all <- function(series, obs, rng) {
    out <- vector("list", 12)
    for (i in 1:12) {
      d <- delta_compute(x = series[[i]], ref_time = 0, obs = obs[[i]])
      out[[i]] <- delta_downscale(x = series[[i]], delta_rast = d,
                                  x_landmask_high = land_high, range_limits = rng)
    }
    terra::sds(out)
  }

  tavg_downscaled <- downscale_all(tavg_series, tavg_obs, tavg_range)
  prec_downscaled <- downscale_all(prec_series, prec_obs, prec_range)

  bioclim_downscaled <- bioclim_vars(tavg = tavg_downscaled, prec = prec_downscaled)

  terra::writeCDF(bioclim_downscaled, model_path, overwrite = TRUE)
  message("Downscaled model written to ", model_path)
}


# =============================================================================
# LOAD CLIMATE MODELS (native + downscaled)
# =============================================================================
downscaled_data <- region_series(
  bio_variables = c("bio01", "bio04", "bio12", "bio15"),
  dataset = "custom", path_to_nc = model_path, time_bp = time_steps)

native_data <- region_series(
  bio_variables = c("bio01", "bio04", "bio12", "bio15"),
  dataset = "Krapp2021", time_bp = time_steps, ext = n_europe_ext)


# --- Precompute Figure S4 inputs (rasters don't survive .rds; cache data frames)
pa_lat  <- coordinates["Pa_24A", "latitude"]
pa_lon  <- coordinates["Pa_24A", "longitude"]
lyr_s4  <- which(time_steps == -34000)
zoom_s4 <- ext(pa_lon - 2.5, pa_lon + 2.5, pa_lat - 2, pa_lat + 2)

.tap_ds  <- crop(downscaled_data[["bio12"]][[lyr_s4]], zoom_s4)
.tap_nat <- crop(native_data[["bio12"]][[lyr_s4]],     zoom_s4)

s4_df_ds  <- as.data.frame(.tap_ds,  xy = TRUE); names(s4_df_ds)[3]  <- "TAP"
s4_df_nat <- as.data.frame(.tap_nat, xy = TRUE); names(s4_df_nat)[3] <- "TAP"
s4_zlim   <- range(c(s4_df_ds$TAP, s4_df_nat$TAP), na.rm = TRUE)

# neighbour comparison (centre cell vs its 8 neighbours), computed while live
.cell <- cellFromXY(.tap_ds, cbind(pa_lon, pa_lat))
.adj9 <- as.vector(adjacent(.tap_ds, .cell, directions = 8, include = TRUE))
s4_neigh <- data.frame(
  position   = c("centre (Pa 24A)", rep("neighbour", length(.adj9) - 1)),
  downscaled = round(as.numeric(.tap_ds[.adj9][, 1])),
  native     = round(as.numeric(.tap_nat[.adj9][, 1])))
s4_point <- data.frame(x = pa_lon, y = pa_lat)


# =============================================================================
# SITE COORDINATES AND DATES
# =============================================================================
# Rows are the six analysed assemblages, in display order.
# Bombrini uses its original coordinates; pastclim's
# nearest-neighbour interpolation handles the ocean-flagged cell downstream.
coordinates <- rbind(
  Pa_24A        = c(41.654115,  15.6149763),  # Paglicci 24A
  Ca_Au         = c(40.000983,  15.381251),   # Grotta della Cala
  `CTC_gic-ars` = c(40.5419159, 15.3093058),  # Castelcivita gic-ars
  `CTC_rsa'`    = c(40.5419159, 15.3093058),  # Castelcivita rsa'
  RB_A2         = c(43.7838361, 7.5341963),   # Riparo Bombrini A2
  `RF_A2-A1`    = c(45.5918287, 10.9051824)   # Riparo Fumane A2-A1
)
colnames(coordinates) <- c("latitude", "longitude")

# Bayesian model MEDIAN ages (years BP) — the best chronological estimate.
# (See Table 1 / SI for provenance of each age.)
site_dates <- c(
  Pa_24A        = 34437,
  Ca_Au         = 38028,
  `CTC_gic-ars` = 39720,
  `CTC_rsa'`    = 40883,
  RB_A2         = 39919,
  `RF_A2-A1`    = 40532
)

group_order <- names(site_dates)
stopifnot(identical(rownames(coordinates), group_order))


# =============================================================================
# PREDICTOR DISTANCE MATRICES
# =============================================================================

## ---- Chronological distance -------------------------------------------------
chrono_dist <- as.matrix(dist(site_dates))
dimnames(chrono_dist) <- list(group_order, group_order)

## ---- Environmental extraction (native + downscaled) -------------------------
# Each assemblage queried at the climate slice nearest its median age.
coords_df <- as.data.frame(coordinates)
coords_df$time_bp <- round(site_dates / 1000) * -1000

clim_native <- location_slice(
  x = coords_df, coords = c("longitude", "latitude"),
  bio_variables = c("bio01", "bio04", "bio12", "bio15", "npp"),
  dataset = "Krapp2021")

clim_downscaled <- location_slice(
  x = coords_df, coords = c("longitude", "latitude"),
  bio_variables = c("bio01", "bio04", "bio12", "bio15"),
  dataset = "custom", path_to_nc = model_path)

# location_slice() prepends id/coord columns; keep only the bio* / npp columns.
clim_native     <- clim_native[,     c("bio01", "bio04", "bio12", "bio15", "npp")]
clim_downscaled <- clim_downscaled[, c("bio01", "bio04", "bio12", "bio15")]


# =============================================================================
# TABLE 5 — environmental conditions at each site (native + downscaled)
# =============================================================================
table5 <- data.frame(
  Assemblage = group_order,
  MAT        = round(clim_native$bio01,     1),
  MAT_ds     = round(clim_downscaled$bio01, 1),
  TAP        = round(clim_native$bio12),
  TAP_ds     = round(clim_downscaled$bio12),
  NPP        = round(clim_native$npp),
  row.names  = NULL
)
cat("\n=== Table 5: environmental conditions per assemblage ===\n")
print(table5, row.names = FALSE)


## ---- Environmental distance matrices ---------------------------------------
env_dist <- function(v) {
  m <- as.matrix(dist(v)); dimnames(m) <- list(group_order, group_order); m
}
MAT    <- env_dist(clim_native$bio01)      # mean annual temperature
TS     <- env_dist(clim_native$bio04)      # temperature seasonality
TAP    <- env_dist(clim_native$bio12)      # total annual precipitation
PS     <- env_dist(clim_native$bio15)      # precipitation seasonality
npp    <- env_dist(clim_native$npp)        # net primary productivity
MAT_ds <- env_dist(clim_downscaled$bio01)  # MAT (downscaled)
TS_ds  <- env_dist(clim_downscaled$bio04)  # TS  (downscaled)
TAP_ds <- env_dist(clim_downscaled$bio12)  # TAP (downscaled)
PS_ds  <- env_dist(clim_downscaled$bio15)  # PS  (downscaled)

## ---- Least-cost-path (Tobler) distance -------------------------------------
relief_rast  <- load_etopo()            # download_etopo() first if not present
Italy_extent <- c(5, 20, 36, 47)
Italy        <- crop(relief_rast, Italy_extent)

sites_vect <- vect(cbind(coordinates[, "longitude"], coordinates[, "latitude"]),
                   crs = crs(relief_rast))

# Palaeo land mask for the target window, then project to metres for gdistance.
land_mask    <- make_land_mask(relief_rast = Italy,
                               time_bp = seq(-38000, -42000, -1000))
Italy_masked <- mask(Italy, land_mask, inverse = FALSE)
Italy_proj   <- project(Italy_masked[[3]],
  "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
Italy_raster <- raster(Italy_proj)      # gdistance needs Raster* objects

sites_sp      <- as(sites_vect, "Spatial")
sites_sp_proj <- spTransform(sites_sp, CRS(proj4string(Italy_raster)))

inside_idx <- !is.na(cellFromXY(Italy_raster, coordinates(sites_sp_proj)))
stopifnot(sum(inside_idx) == length(group_order))   # all sites must land on-grid

sites_snapped <- SpatialPoints(
  xyFromCell(Italy_raster,
             cellFromXY(Italy_raster, coordinates(sites_sp_proj)[inside_idx, ])),
  proj4string = CRS(proj4string(Italy_raster)))

# Conductance = Tobler's hiking function over slope (altitude difference).
altDiff     <- function(x) x[2] - x[1]
hd          <- transition(Italy_raster, altDiff, 8, symm = TRUE)
slope       <- geoCorrection(hd, type = "c")
adj         <- adjacencyFromTransition(slope)
speed       <- slope
speed[adj]  <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))   # Tobler
conductance <- geoCorrection(speed, type = "c")

costpath_dist <- as.matrix(dist(costDistance(conductance, sites_snapped)))
dimnames(costpath_dist) <- list(group_order, group_order)


# =============================================================================
# PREDICTOR PANEL (reduced set used in the headline analyses)
# =============================================================================
# Seasonality terms (TS, PS, and their _ds) are dropped from the headline
# Mantel panel; they remain available above if you want the full variant.
predictor_mats <- list(
  chrono   = chrono_dist,
  costpath = costpath_dist,
  MAT      = MAT,  MAT_ds = MAT_ds,
  TAP      = TAP,  TAP_ds = TAP_ds,
  npp      = npp
)

pretty_pred <- c(chrono = "Chronology", costpath = "Costpath",
                 MAT = "MAT", MAT_ds = "MAT (ds)",
                 TAP = "TAP", TAP_ds = "TAP (ds)", npp = "NPP")


# =============================================================================
# SAVE for downstream scripts
# =============================================================================
saveRDS(
  list(
    coordinates     = coordinates,
    site_dates      = site_dates,
    group_order     = group_order,
    clim_native     = clim_native,
    clim_downscaled = clim_downscaled,
    predictor_mats  = predictor_mats,
    pretty_pred     = pretty_pred,
    # individual matrices kept for convenience / supplementary use
    chrono_dist = chrono_dist, costpath_dist = costpath_dist,
    MAT = MAT, MAT_ds = MAT_ds, TAP = TAP, TAP_ds = TAP_ds,
    TS = TS, TS_ds = TS_ds, PS = PS, PS_ds = PS_ds, npp = npp,
    # climate rasters for Figure S4 (single-cell artefact check)
    downscaled_data = downscaled_data, native_data = native_data,
    model_path = model_path, time_steps = time_steps,
    s4_df_ds = s4_df_ds, s4_df_nat = s4_df_nat, s4_zlim = s4_zlim,
    s4_neigh = s4_neigh, s4_point = s4_point
  ),
  "output/rds/01_distances.rds"
)


# =============================================================================
# (OPTIONAL) FIGURE: native vs downscaled climate time series per site
# =============================================================================
# Not a numbered main-text figure; useful diagnostic. Set to TRUE to draw.
if (FALSE) {
  make_ts <- function(dataset, path_to_nc = NULL, res_label) {
    args <- list(x = as.data.frame(coordinates), time_bp = time_steps,
                 bio_variables = c("bio01", "bio04", "bio12", "bio15"),
                 dataset = dataset)
    if (!is.null(path_to_nc)) args$path_to_nc <- path_to_nc
    out <- do.call(location_series, args)
    out$resolution <- res_label
    out$site <- factor(out$name, levels = seq_along(group_order),
                       labels = group_order)
    out
  }
  ts_both <- rbind(
    make_ts("Krapp2021", NULL,        "native (0.5 deg)"),
    make_ts("custom",    model_path,  "downscaled (10-min)")
  )

  panel <- function(var, title) {
    ggplot(ts_both, aes(time_bp, .data[[var]],
                        group = interaction(site, resolution),
                        colour = site, linetype = resolution)) +
      geom_line(linewidth = 0.8) +
      labs(title = title, x = "Time (years BP)", y = var,
           colour = "Site", linetype = "Resolution") +
      theme_bw()
  }
  ((panel("bio01", "Mean annual temperature") |
    panel("bio12", "Total annual precipitation")) /
   (panel("bio04", "Temperature seasonality") |
    panel("bio15", "Precipitation seasonality"))) +
    plot_layout(guides = "collect") & theme(legend.position = "right")
}
