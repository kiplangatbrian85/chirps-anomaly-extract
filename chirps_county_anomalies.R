# =============================================================================
# chirps_county_anomalies.R
#
# Seasonal rainfall anomalies at sub-national level (ADM1) for East Africa,
# using CHIRPS v2.0 monthly rainfall, terra, sf and exactextractr.
#
# Steps:
#   1. Download and cache CHIRPS monthly GeoTIFFs for a target season plus a
#      baseline climatology period.
#   2. Build seasonal rainfall totals per year and compare the target season
#      against the baseline (anomaly in mm, % of normal, standardized anomaly).
#   3. Extract area-weighted means per ADM1 unit with exactextractr, plus the
#      share of each unit's area in "dry" conditions (standardized anomaly < -1).
#   4. Run quality checks and write tidy CSV outputs, a GeoPackage, and a map.
#
# The script handles several countries, so the same code can run across a
# multi-country portfolio (e.g. KEN, SOM, ETH, SSD, UGA, TZA).
#
# Author: Brian Kiplangat Kilele
# =============================================================================

# ---- 0. Configuration --------------------------------------------------------

config <- list(
  countries      = c("KEN"),            # ISO3 codes, e.g. c("KEN","SOM","ETH")
  season_name    = "OND",               # label used in outputs
  season_months  = 10:12,               # Oct-Nov-Dec short rains
  target_year    = 2025,
  baseline_years = 1991:2020,           # WMO standard normal period
  dry_threshold  = -1,                  # standardized anomaly for "dry" area share
  chirps_url     = "https://data.chc.ucsb.edu/products/CHIRPS-2.0/africa_monthly/tifs",
  cache_dir      = "data/chirps_monthly",
  out_dir        = "outputs"
)

# ---- 1. Packages -------------------------------------------------------------

pkgs <- c("terra", "sf", "exactextractr", "dplyr", "tidyr", "readr",
          "ggplot2", "jsonlite", "R.utils")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) install.packages(missing, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(exactextractr)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

dir.create(config$cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$out_dir,   recursive = TRUE, showWarnings = FALSE)
options(timeout = 600)

log_msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", ...)

# ---- 2. Helpers --------------------------------------------------------------

# Download one CHIRPS monthly file (gzipped GeoTIFF) and cache the unzipped tif.
get_chirps_month <- function(year, month, cfg = config) {
  fname <- sprintf("chirps-v2.0.%d.%02d.tif", year, month)
  tif   <- file.path(cfg$cache_dir, fname)
  if (file.exists(tif)) return(tif)

  gz  <- paste0(tif, ".gz")
  url <- paste0(cfg$chirps_url, "/", fname, ".gz")
  ok  <- tryCatch({
    download.file(url, gz, mode = "wb", quiet = TRUE)
    R.utils::gunzip(gz, destname = tif, remove = TRUE, overwrite = TRUE)
    TRUE
  }, error = function(e) {
    warning("Failed to fetch ", url, ": ", conditionMessage(e))
    FALSE
  })
  if (ok) tif else NA_character_
}

# Seasonal total for one year, cropped to the area of interest.
seasonal_total <- function(year, aoi, cfg = config) {
  files <- vapply(cfg$season_months, function(m) get_chirps_month(year, m, cfg),
                  character(1))
  if (anyNA(files)) {
    warning("Incomplete season for ", year, "; skipped")
    return(NULL)
  }
  r <- rast(files)
  r <- crop(r, aoi, snap = "out")
  r <- classify(r, cbind(-Inf, 0, NA), right = FALSE)  # nodata (-9999) -> NA
  s <- app(r, sum, na.rm = FALSE)                       # a missing month -> NA
  names(s) <- as.character(year)
  s
}

# ADM1 boundaries from the geoBoundaries API (open licence, 47 counties for KEN).
get_adm1 <- function(iso3) {
  api  <- sprintf("https://www.geoboundaries.org/api/current/gbOpen/%s/ADM1/", iso3)
  meta <- jsonlite::fromJSON(api)
  b    <- st_read(meta$gjDownloadURL, quiet = TRUE)
  b    |>
    st_make_valid() |>
    st_transform(4326) |>
    transmute(iso3 = iso3, adm1_name = shapeName, adm1_id = shapeID)
}

# ---- 3. Main loop over countries ---------------------------------------------

results <- list()

for (iso3 in config$countries) {
  log_msg("Processing ", iso3)

  adm1 <- get_adm1(iso3)
  aoi  <- ext(vect(adm1)) + 0.1   # small buffer so edge pixels are kept

  # 3a. Baseline seasonal totals and climatology
  log_msg("  Building baseline ", min(config$baseline_years), "-",
          max(config$baseline_years))
  base_list <- lapply(config$baseline_years, seasonal_total, aoi = aoi)
  base_list <- Filter(Negate(is.null), base_list)
  base      <- rast(base_list)

  n_base <- nlyr(base)
  if (n_base < 0.9 * length(config$baseline_years)) {
    stop(iso3, ": only ", n_base, " baseline years available")
  }

  clim_mean <- app(base, mean, na.rm = TRUE)
  clim_sd   <- app(base, sd,   na.rm = TRUE)

  # 3b. Target season and anomaly layers
  log_msg("  Target season ", config$season_name, " ", config$target_year)
  target <- seasonal_total(config$target_year, aoi)
  if (is.null(target)) stop(iso3, ": target season data unavailable")

  anom_mm  <- target - clim_mean
  pct_norm <- ifel(clim_mean > 1, 100 * target / clim_mean, NA)  # avoid ~0 means
  std_anom <- ifel(clim_sd > 0, anom_mm / clim_sd, NA)
  is_dry   <- std_anom < config$dry_threshold

  stack <- c(target, clim_mean, anom_mm, pct_norm, std_anom, is_dry)
  names(stack) <- c("rain_mm", "normal_mm", "anom_mm", "pct_normal",
                    "std_anom", "dry_frac")

  # 3c. Zonal statistics: area-weighted mean per ADM1 unit.
  # exactextractr weights each pixel by the fraction covered by the polygon,
  # so small counties and irregular boundaries are handled properly.
  # The mean of the 0/1 dry layer gives the share of area classed as dry.
  zs <- exact_extract(stack, adm1, fun = "mean", progress = FALSE)
  names(zs) <- sub("^mean\\.", "", names(zs))

  county_tbl <- adm1 |>
    st_drop_geometry() |>
    bind_cols(zs) |>
    mutate(
      season       = config$season_name,
      year         = config$target_year,
      baseline     = paste0(min(config$baseline_years), "-",
                            max(config$baseline_years)),
      dry_frac     = round(dry_frac, 3),
      anomaly_class = cut(
        std_anom,
        breaks = c(-Inf, -1.5, -1, -0.5, 0.5, 1, 1.5, Inf),
        labels = c("Very dry", "Dry", "Moderately dry", "Near normal",
                   "Moderately wet", "Wet", "Very wet")
      )
    ) |>
    mutate(across(c(rain_mm, normal_mm, anom_mm, pct_normal, std_anom),
                  ~ round(.x, 2)))

  # ---- 4. Quality checks -----------------------------------------------------
  qa <- tibble(
    iso3               = iso3,
    n_units            = nrow(adm1),
    n_units_with_value = sum(!is.na(county_tbl$std_anom)),
    baseline_years     = n_base,
    min_rain_mm        = min(county_tbl$rain_mm, na.rm = TRUE),
    max_rain_mm        = max(county_tbl$rain_mm, na.rm = TRUE),
    negative_rain      = any(county_tbl$rain_mm < 0, na.rm = TRUE),
    extreme_std_anom   = sum(abs(county_tbl$std_anom) > 4, na.rm = TRUE)
  )
  if (qa$negative_rain) warning(iso3, ": negative rainfall found, check nodata handling")
  if (qa$n_units_with_value < qa$n_units) {
    warning(iso3, ": ", qa$n_units - qa$n_units_with_value,
            " unit(s) have no value (check small or coastal units)")
  }

  results[[iso3]] <- list(table = county_tbl, adm1 = adm1, qa = qa,
                          std_anom = std_anom)
}

# ---- 5. Outputs --------------------------------------------------------------

tag <- sprintf("%s_%d", config$season_name, config$target_year)

all_tbl <- bind_rows(lapply(results, `[[`, "table"))
all_qa  <- bind_rows(lapply(results, `[[`, "qa"))

write_csv(all_tbl, file.path(config$out_dir, paste0("adm1_rainfall_anomalies_", tag, ".csv")))
write_csv(all_qa,  file.path(config$out_dir, paste0("qa_summary_", tag, ".csv")))

all_sf <- bind_rows(lapply(results, `[[`, "adm1")) |>
  left_join(all_tbl, by = c("iso3", "adm1_name", "adm1_id"))
st_write(all_sf, file.path(config$out_dir, paste0("adm1_rainfall_anomalies_", tag, ".gpkg")),
         delete_dsn = TRUE, quiet = TRUE)

# Choropleth of standardized anomaly with a diverging brown-to-green palette
p <- ggplot(all_sf) +
  geom_sf(aes(fill = std_anom), colour = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#8c510a", mid = "#f5f5f5", high = "#01665e",
                       midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish,
                       name = "Standardized\nanomaly") +
  labs(
    title    = sprintf("%s %d rainfall anomaly by ADM1 unit", config$season_name,
                       config$target_year),
    subtitle = sprintf("CHIRPS v2.0, baseline %d-%d, area-weighted (exactextractr)",
                       min(config$baseline_years), max(config$baseline_years)),
    caption  = "Boundaries: geoBoundaries gbOpen ADM1"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text = element_blank())

ggsave(file.path(config$out_dir, paste0("map_std_anomaly_", tag, ".png")),
       p, width = 8, height = 8, dpi = 200, bg = "white")

log_msg("Done. Outputs written to ", normalizePath(config$out_dir))
print(all_qa)
print(all_tbl |> arrange(std_anom) |> select(iso3, adm1_name, rain_mm,
                                             pct_normal, std_anom, dry_frac) |>
        head(10))
