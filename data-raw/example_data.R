## Code to prepare the example data returned by windscape_example().
##
## All three data sets derive from CFSR 10 m hourly wind, originally downloaded with cfsr_dl()
## before NCAR's archive moved (Sept 2025), which broke those download URLs:
##   wind  <- cfsr_dl(years = 2000, months = 1:12, days = c(1, 15), xlim = c(-120, -90) + 360,  # approx.; exact call not recorded
##                    ylim = c(30, 50)) %>% shift(dx = -360)   # saved as legacy wind.tif
##   field <- cfsr_dl(years = 2005, months = 8, days = 28, hlim = c(23, 23),
##                    xlim = c(-99, -78) + 360, ylim = c(17, 35)) %>% shift(dx = -360)  # katrina
## Below, `legacy` holds those original rasters.

library(terra)
devtools::load_all()
legacy <- "/tmp/legacy"
out <- "inst/extdata"
gdal_int <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9")
gdal_flt <- c("COMPRESS=DEFLATE", "PREDICTOR=3", "ZLEVEL=9")

# hourly u and v, 1st and 15th of each month in 2000 (576 hours), western/central US
w <- rast(file.path(legacy, "wind.tif"))
n <- nlyr(w) / 2
crs(w) <- "EPSG:4326"

# example wind_series: every 6th hour (96 time steps), stored to 0.1 m/s
keep <- seq(1, n, 6)
ws <- c(w[[keep]], w[[n + keep]])
writeRaster(round(ws, 1), file.path(out, "wind_usa.tif"), datatype = "INT2S", scale = 0.1,
            gdal = gdal_int, overwrite = TRUE)

# example wind_rose: built from all 576 hours at full precision, so it is a better estimate of
# the wind regime than a rose built from the 96-step example series
rose <- wind_rose(wind_series(w), trans = 1)
writeRaster(rose, file.path(out, "rose_usa.tif"), datatype = "FLT4S", gdal = gdal_flt,
            overwrite = TRUE)

# example wind_field: Hurricane Katrina, 2005-08-28 23:00 UTC
load(file.path(legacy, "katrina.rda"))
k <- unwrap(katrina)
crs(k) <- "EPSG:4326"
writeRaster(round(k, 1), file.path(out, "katrina.tif"), datatype = "INT2S", scale = 0.1,
            gdal = gdal_int, overwrite = TRUE)
