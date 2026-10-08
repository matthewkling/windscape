## Code to build the pre-made global CFSR wind roses hosted for download.
##
## Products (all trans = 1, 10 m wind, CFSR 1979-2010):
##   month          one rose per year-month (384)
##   year           one rose per calendar year (32)
##   month_of_year  one rose per calendar month, across all years (12)
##   all            the full period (1)
## Aggregates are hour-weighted combinations of the monthly roses (combine_roses()), so each equals
## the rose built from all of its hours at once.
##
## Monthly roses were built from hourly NCAR d093001 files (wnd10m.gdas.YYYYMM.grb2) with
##   wind_rose(wind_series(f, order = "uvuv"), trans = 1)
## 1979-1989 and 1991-2010 in May 2024; 1990 rebuilt Oct 2026 from fresh downloads (several
## original 1990 files were corrupt). Months that were fine in both batches matched exactly.
## Each monthly file holds one time step per hour (checked: n_steps = hours in month).
##
## Output: cloud-optimized GeoTIFFs (pixel-interleaved 256 x 256 tiles, DEFLATE with the
## floating-point predictor, no overviews), longitudes -180 to 180, with windscape_* GDAL
## metadata (read with terra::describe(f, meta = TRUE)); plus catalog.csv. Requires terra >= 1.8-42,
## which writes metadata tags (terra::metags()) into the file itself.

library(terra)
devtools::load_all()

src_dir <- "~/data/CFSR/monthly_roses"
out_dir <- "~/data/CFSR/windscape_roses"
years <- 1979:2010
source <- "cfsr"
level <- "10m"
format_version <- 1   # bump if the meaning of rose values or the metadata ever changes

stopifnot(packageVersion("terra") >= "1.8-42")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# monthly inputs --------------------------------------------------------------------------------

ym <- expand.grid(month = 1:12, year = years)
ym$file <- file.path(src_dir, sprintf("wnd10m_%d%02d.tif", ym$year, ym$month))
stopifnot(all(file.exists(path.expand(ym$file))))
first_day <- as.Date(sprintf("%d-%02d-01", ym$year, ym$month))
next_month <- as.Date(sprintf("%d-%02d-01", ym$year + (ym$month == 12), ym$month %% 12 + 1))
ym$n_steps <- as.integer(next_month - first_day) * 24L

monthly <- Map(function(f, n) wind_rose(f, trans = 1, n_steps = n), ym$file, ym$n_steps)


# writing ---------------------------------------------------------------------------------------

catalog <- list()

# Write a rose as a COG with windscape metadata, and add it to the catalog. `tier`, `year`,
# `month` describe the product; `first` and `last` are the first and last months it covers.
write_rose <- function(rose, tier, year = NA, month = NA, first, last){
      name <- paste0(paste(c(source, level, tier,
                             if(!is.na(year)) year,
                             if(!is.na(month)) sprintf("%02d", month)), collapse = "_"), ".tif")
      r <- methods::as(rose, "SpatRaster")
      if(xmax(r) > 180 + xres(r)) r <- rotate(r)   # 0 to 360 -> -180 to 180
      stopifnot(abs(xmax(r) - xmin(r) - 360) < xres(r) / 2, xmax(r) <= 180 + xres(r))
      meta <- c(rose_format = format_version, source = source, level = level, tier = tier,
                trans = 1, n_steps = rose@n_steps, first = first, last = last,
                units = "1/hour")
      metags(r) <- setNames(as.character(meta), paste0("windscape_", names(meta)))
      out <- file.path(path.expand(out_dir), name)
      writeRaster(r, out, filetype = "COG", datatype = "FLT4S", overwrite = TRUE,
                  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=YES", "BLOCKSIZE=256", "OVERVIEWS=NONE"))

      catalog[[length(catalog) + 1]] <<- data.frame(
            source = source, level = level, tier = tier, year = year, month = month,
            first = first, last = last, n_steps = rose@n_steps, trans = 1, file = name,
            bytes = file.size(out), md5 = unname(tools::md5sum(out)))
      invisible(out)
}

ymstr <- function(y, m) sprintf("%d-%02d", y, m)
span <- c(ymstr(min(years), 1), ymstr(max(years), 12))

for(i in seq_len(nrow(ym))){
      write_rose(monthly[[i]], "month", ym$year[i], ym$month[i],
                 ymstr(ym$year[i], ym$month[i]), ymstr(ym$year[i], ym$month[i]))
}

annual <- lapply(years, function(y) combine_roses(monthly[ym$year == y]))
for(i in seq_along(years)){
      write_rose(annual[[i]], "year", years[i], NA, ymstr(years[i], 1), ymstr(years[i], 12))
}

moy <- lapply(1:12, function(m) combine_roses(monthly[ym$month == m]))
for(m in 1:12) write_rose(moy[[m]], "month_of_year", NA, m, span[1], span[2])

full <- combine_roses(annual)
write_rose(full, "all", NA, NA, span[1], span[2])

catalog <- do.call(rbind, catalog)
write.csv(catalog, file.path(path.expand(out_dir), "catalog.csv"), row.names = FALSE)


# checks ----------------------------------------------------------------------------------------

# the full period two ways: from annual roses and from month-of-year roses
alt <- combine_roses(moy)
stopifnot(full@n_steps == sum(ym$n_steps), alt@n_steps == full@n_steps)
rel_diff <- function(a, b) max(global(abs(a - b), "max", na.rm = TRUE)[, 1]) /
      max(global(abs(b), "max", na.rm = TRUE)[, 1])
stopifnot(rel_diff(full, alt) < 1e-6)

# a written file reads back with its values and metadata
f <- file.path(path.expand(out_dir), catalog$file[catalog$tier == "all"])
meta <- describe(f, meta = TRUE)
stopifnot(paste0("windscape_n_steps=", full@n_steps) %in% meta)
back <- rast(f)
stopifnot(identical(names(back), c("SW", "W", "NW", "N", "NE", "E", "SE", "S")))
r <- methods::as(full, "SpatRaster")
if(xmax(r) > 180 + xres(r)) r <- rotate(r)
stopifnot(compareGeom(back, r), rel_diff(back, r) < 1e-6)

table(catalog$tier)
sum(catalog$bytes) / 1e9   # GB
