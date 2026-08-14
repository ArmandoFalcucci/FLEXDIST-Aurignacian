# =============================================================================
# Paglicci climate at three candidate time slices (-35, -34, -33 ka)
# Mirrors the extraction in 01_distances.R exactly: location_slice() at the
# site coordinates, native from Krapp2021 and downscaled from the local .nc.
# Run from the project root.
# =============================================================================
library(pastclim)

model_path <- "data/Europe_bioclims_10min.nc"
pa_lon <- 15.615
pa_lat <- 41.654
slices <- c(-35000, -34000, -33000)

# location_slice() needs one row per time point
pa_df <- data.frame(longitude = rep(pa_lon, length(slices)),
                    latitude  = rep(pa_lat, length(slices)),
                    time_bp   = slices)

nat <- location_slice(
  x = pa_df, coords = c("longitude", "latitude"),
  bio_variables = c("bio01", "bio04", "bio12", "bio15", "npp"),
  dataset = "Krapp2021")

ds <- location_slice(
  x = pa_df, coords = c("longitude", "latitude"),
  bio_variables = c("bio01", "bio04", "bio12", "bio15"),
  dataset = "custom", path_to_nc = model_path)

out <- data.frame(
  slice  = slices,
  MAT    = round(nat$bio01, 2),
  MAT_ds = round(ds$bio01,  2),
  TAP    = round(nat$bio12),
  TAP_ds = round(ds$bio12),
  NPP    = round(nat$npp))

print(out, row.names = FALSE)

cat("\nTable 5 reports, at -34000: MAT 8.8 | MAT_ds 10.3 | TAP 447 | TAP_ds 278 | NPP 235\n")
cat("The -34000 row above must match those values.\n")