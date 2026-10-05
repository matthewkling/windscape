# Download wind data from NCAR's Geoscience Data Exchange (GDEX) THREDDS server, using the
# NetCDF Subset Service (NCSS) to clip each request to a bounding box on the server. Each month
# of data is converted to a GeoTIFF that can be loaded with wind_series(), and cached.


#' Download hourly wind data from NCAR
#'
#' Downloads gridded hourly wind data from the NCAR Geoscience Data Exchange (GDEX) for a
#' bounding box and a set of months, saving one GeoTIFF file per month. No account is needed.
#' Data are clipped to the bounding box on the server, so only the requested region is
#' transferred.
#'
#' Each output file holds one month of data in `wind_series` layout: all u layers followed by
#' all v layers, named like `"u 2005-08-28 23:00:00"` (UTC), in m/s, oriented to true east and
#' north, on a longitude/latitude grid. Combine the files into one series with [wind_series()],
#' which keeps the data on disk, and summarize it with [wind_rose()].
#'
#' Downloads are cached: a month whose file already exists in `dir` is not downloaded again
#' unless `overwrite = TRUE`, so an interrupted download can be resumed by rerunning the same
#' call. Large requests (many years, or a large region) can take a long time and use
#' substantial disk space; thinning hourly data with `time_stride` reduces both.
#'
#' Requires the \pkg{ncdf4} package.
#'
#' @param source Data set to download:
#'   * `"era5"` (the default): ECMWF ERA5 reanalysis, 0.25 degree grid, 1940 to present
#'     (NCAR data set d633000).
#'   * `"cfsr"`: NCEP Climate Forecast System Reanalysis, ~0.31 degree Gaussian grid, 1979 to
#'     2010 (d093001).
#'   * `"cfsv2"`: NCEP Climate Forecast System version 2, the operational continuation of CFSR,
#'     ~0.2 degree Gaussian grid, 2011 to present (d094001).
#' @param level Height of the wind data. `"10m"` (the default) is wind 10 m above the ground,
#'   available for all sources. `"100m"` is available for ERA5 only. Pressure levels
#'   `"1000hPa"`, `"850hPa"`, `"700hPa"`, `"500hPa"`, and `"200hPa"` are available for CFSR and
#'   CFSv2 only.
#' @param xlim Numeric vector of length 2 giving the western and eastern limits of the bounding
#'   box, in degrees. Longitudes can be given either in the -180 to 180 range or the 0 to 360
#'   range, and output files use the same convention. Boxes crossing the prime meridian (in -180
#'   to 180 coordinates, e.g. `c(-10, 30)`) or the antimeridian (in 0 to 360 coordinates, e.g.
#'   `c(170, 200)`) are supported.
#' @param ylim Numeric vector of length 2 giving the southern and northern limits of the bounding
#'   box, in degrees, in the range -90 to 90.
#' @param years Integer vector of years to download.
#' @param months Integer vector of months (1-12) to download within each year. Defaults to all
#'   months.
#' @param time_stride Integer. Download every `time_stride`-th hourly time step, starting with
#'   each month's first time step. The default, 1, downloads all hours; for example, 3 downloads
#'   every third hour and reduces download time and file size about threefold.
#' @param dir Directory where monthly files are saved. Defaults to a temporary directory, which
#'   is deleted at the end of the R session; set it to a permanent location to keep the data and
#'   make use of caching across sessions.
#' @param overwrite Logical. If `FALSE` (the default), months already downloaded to `dir` are
#'   skipped.
#' @param quiet Logical. If `TRUE`, progress messages are suppressed.
#' @return A character vector of file paths, one per month, in chronological order (returned
#'   invisibly if all files were already cached).
#' @seealso [wind_series()] to load the files; [wind_rose()] to summarize them;
#'   [ncar_land()] to download a matching land-water layer.
#' @examples
#' \dontrun{
#' # every third hour of 10 m ERA5 wind for summer 2020, Pacific Northwest
#' files <- ncar_download("era5", xlim = c(-125, -115), ylim = c(42, 49),
#'                        years = 2020, months = 6:8, time_stride = 3,
#'                        dir = "~/wind_data")
#' ws <- wind_series(files)
#' rose <- wind_rose(files)
#' }
#' @export
ncar_download <- function(source = c("era5", "cfsr", "cfsv2"), level = "10m",
                          xlim, ylim, years, months = 1:12, time_stride = 1,
                          dir = tempdir(), overwrite = FALSE, quiet = FALSE){

      source <- match.arg(source)
      if(!requireNamespace("ncdf4", quietly = TRUE))
            stop("the ncdf4 package is required to download NCAR data", call. = FALSE)
      spec <- ncar_spec(source, level)
      bbox <- check_bbox(xlim, ylim)
      if(length(time_stride) != 1 || !is.finite(time_stride) || time_stride < 1 ||
         time_stride %% 1 != 0) stop("`time_stride` must be a positive integer", call. = FALSE)
      if(length(months) == 0 || any(!months %in% 1:12))
            stop("`months` must be integers between 1 and 12", call. = FALSE)
      if(length(years) == 0 || any(years %% 1 != 0))
            stop("`years` must be integers", call. = FALSE)
      ym <- expand.grid(month = sort(unique(months)), year = sort(unique(years)))
      check_years(ym, spec)
      dir.create(dir, showWarnings = FALSE, recursive = TRUE)

      files <- file.path(dir, sprintf("%s_%s_%04d%02d_%s%s.tif", source, level, ym$year, ym$month,
                                      bbox$tag, ifelse(time_stride > 1, paste0("_t", time_stride), "")))
      todo <- overwrite | !file.exists(files)
      for(i in which(todo)){
            if(!quiet) message(sprintf("downloading %s %s %04d-%02d (%d of %d)", source, level,
                                       ym$year[i], ym$month[i], sum(todo[seq_len(i)]), sum(todo)))
            ncar_month(spec, ym$year[i], ym$month[i], bbox, time_stride, files[i])
      }
      if(!quiet && any(!todo)) message(sum(!todo), " month(s) already downloaded; skipped")
      if(any(todo)) files else invisible(files)
}


#' Download a land-water layer from NCAR
#'
#' Downloads a land layer matching the grid of [ncar_download()] data for the same `source` and
#' bounding box, e.g. for weighting a wind rose to reduce conductance over water.
#'
#' @param source Data set: `"era5"` or `"cfsr"`. (Not yet available for `"cfsv2"`.)
#' @param xlim,ylim Bounding box; see [ncar_download()].
#' @return A single-layer `SpatRaster` named `"land"`. For ERA5 this is the land fraction of each
#'   cell, from 0 (all water) to 1 (all land); use e.g. `land >= 0.5` for a binary layer. For CFSR
#'   it is binary, 1 for land and 0 for water. (CFSR publishes no land mask, so it is derived from
#'   which cells have soil temperature data.)
#' @examples
#' \dontrun{
#' land <- ncar_land("era5", xlim = c(-125, -115), ylim = c(42, 49))
#' }
#' @export
ncar_land <- function(source = c("era5", "cfsr", "cfsv2"), xlim, ylim){
      source <- match.arg(source)
      if(!requireNamespace("ncdf4", quietly = TRUE))
            stop("the ncdf4 package is required to download NCAR data", call. = FALSE)
      bbox <- check_bbox(xlim, ylim)
      req <- switch(source,
                    era5 = list(path = paste0("files/g/d633000/e5.oper.invariant/197901/",
                                              "e5.oper.invariant.128_172_lsm.ll025sc.1979010100_1979010100.nc"),
                                vars = c(land = "LSM")),
                    cfsr = list(path = "files/g/d093001/1980/soilt1.gdas.198001.grb2",
                                vars = c(land = "Temperature_depth_below_surface_layer")),
                    cfsv2 = stop("ncar_land() is not yet available for CFSv2", call. = FALSE))
      g <- ncar_fetch_grid(req, bbox, time_stride = 1)
      x <- grid_to_rast(g, "land")[[1]]
      if(source == "cfsr") x <- terra::ifel(is.na(x), 0, 1)
      names(x) <- "land"
      x
}


# Data set definitions -------------------------------------------------------------------------

# Base URL of NCAR's THREDDS server; an option, so it can be changed if the server moves again.
ncar_thredds <- function() getOption("windscape.ncar_thredds", "https://thredds.rda.ucar.edu/thredds")

ncar_spec <- function(source, level){
      levels <- switch(source,
                       era5 = c("10m", "100m"),
                       cfsr = , cfsv2 = c("10m", "1000hPa", "850hPa", "700hPa", "500hPa", "200hPa"))
      if(length(level) != 1 || !level %in% levels)
            stop("`level` for ", source, " must be one of: ", paste0('"', levels, '"', collapse = ", "),
                 call. = FALSE)

      if(source == "era5"){
            code <- switch(level, "10m" = c(u = "128_165_10u", v = "128_166_10v"),
                           "100m" = c(u = "228_246_100u", v = "228_247_100v"))
            var <- switch(level, "10m" = c(u = "VAR_10U", v = "VAR_10V"),
                          "100m" = c(u = "VAR_100U", v = "VAR_100V"))
            requests <- function(year, month){
                  ym <- sprintf("%04d%02d", year, month)
                  path <- sprintf("files/g/d633000/e5.oper.an.sfc/%s/e5.oper.an.sfc.%s.ll025sc.%s0100_%s%02d23.nc",
                                  ym, code, ym, ym, days_in_month(year, month))
                  list(list(path = path[1], vars = var["u"]),
                       list(path = path[2], vars = var["v"]))
            }
            return(list(source = source, years = c(1940, current_year()), requests = requests))
      }

      # CFSR and CFSv2 share a file structure, with u and v in one file
      file <- if(level == "10m") "wnd10m" else paste0("wnd", sub("hPa", "", level))
      tag <- if(level == "10m") "height_above_ground" else "isobaric"
      var <- c(u = paste0("u-component_of_wind_", tag), v = paste0("v-component_of_wind_", tag))
      ds <- switch(source, cfsr = "d093001", cfsv2 = "d094001")
      stream <- switch(source, cfsr = "gdas", cfsv2 = "cdas1")
      requests <- function(year, month){
            list(list(path = sprintf("files/g/%s/%04d/%s.%s.%04d%02d.grb2", ds, year, file, stream, year, month),
                      vars = var))
      }
      yrs <- switch(source, cfsr = c(1979, 2010), cfsv2 = c(2011, current_year()))
      list(source = source, years = yrs, requests = requests)
}

check_years <- function(ym, spec){
      bad <- ym$year < spec$years[1] | ym$year > spec$years[2]
      if(any(bad)) stop(spec$source, " data are available for ", spec$years[1], "-",
                        if(spec$years[2] == current_year()) "present" else spec$years[2],
                        "; requested years outside this range: ",
                        paste(unique(ym$year[bad]), collapse = ", "), call. = FALSE)
}

current_year <- function() as.integer(format(Sys.Date(), "%Y"))

days_in_month <- function(year, month){
      first <- as.Date(sprintf("%04d-%02d-01", year, month))
      as.integer(format(seq(first, by = "month", length.out = 2)[2] - 1, "%d"))
}


# Bounding boxes and longitude conventions ---------------------------------------------------

# Validate a bounding box and work out how to request it from a 0-360 data set: as one
# longitude range, or as two if it crosses the 0/360 seam. `convention` records whether the
# user gave -180..180 or 0..360 longitudes, which the output will match.
check_bbox <- function(xlim, ylim){
      if(missing(xlim) || missing(ylim)) stop("`xlim` and `ylim` are required", call. = FALSE)
      if(length(xlim) != 2 || length(ylim) != 2 || any(!is.finite(c(xlim, ylim))))
            stop("`xlim` and `ylim` must be numeric vectors of length 2", call. = FALSE)
      x <- sort(xlim)
      y <- sort(ylim)
      if(x[1] < -180 || x[2] > 360) stop("`xlim` must be within -180 to 360", call. = FALSE)
      if(x[2] - x[1] > 360) stop("`xlim` spans more than 360 degrees", call. = FALSE)
      if(x[1] < 0 && x[2] > 180) stop("`xlim` mixes the -180 to 180 and 0 to 360 longitude ",
                                      "conventions; use one or the other", call. = FALSE)
      if(y[1] < -90 || y[2] > 90) stop("`ylim` must be within -90 to 90", call. = FALSE)
      if(x[1] == x[2] || y[1] == y[2]) stop("`xlim` and `ylim` must each span a nonzero range",
                                            call. = FALSE)
      convention <- if(x[2] > 180) "0_360" else "-180_180"

      if(x[2] - x[1] >= 360){
            pieces <- list(c(0, 360))
      }else{
            w <- x[1] %% 360
            e <- x[2] %% 360
            if(e == 0) e <- 360
            pieces <- if(w < e) list(c(w, e)) else list(c(w, 360), c(0, e))
      }
      tag <- gsub("-", "m", sprintf("x%s_%s_y%s_%s", fmt_num(x[1]), fmt_num(x[2]),
                                    fmt_num(y[1]), fmt_num(y[2])))
      list(xlim = x, ylim = y, pieces = pieces, convention = convention, tag = tag)
}

fmt_num <- function(x) format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)

# convert data-set longitudes to the requested convention
convert_lon <- function(lon, convention){
      if(convention == "0_360") lon %% 360 else ((lon + 180) %% 360) - 180
}


# Fetching and assembling -----------------------------------------------------------------------

# Download, assemble, and save one month
ncar_month <- function(spec, year, month, bbox, time_stride, file){
      grids <- lapply(spec$requests(year, month), ncar_fetch_grid, bbox = bbox, time_stride = time_stride)
      g <- do.call(merge_vars, grids)
      x <- grid_to_rast(g, c("u", "v"))
      tmp <- paste0(file, ".part.tif")
      on.exit(unlink(tmp), add = TRUE)
      terra::writeRaster(round(x, 2), tmp, datatype = "INT2S", scale = 0.01, NAflag = -32768,
                         gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"), overwrite = TRUE)
      file.rename(tmp, file) # only now does the month count as downloaded
      invisible(file)
}

# Fetch the variables in one request (a single server file), across one or two longitude
# pieces, and assemble them into one grid: list(lon, lat, time, data = list(var = array[lon, lat, time]))
ncar_fetch_grid <- function(req, bbox, time_stride){
      pieces <- lapply(bbox$pieces, function(p){
            nc <- tempfile(fileext = ".nc")
            on.exit(unlink(nc))
            ncss_get(req$path, req$vars, west = p[1], east = p[2], south = bbox$ylim[1],
                     north = bbox$ylim[2], time_stride = time_stride, dest = nc)
            read_ncss(nc, req$vars)
      })
      join_lon(pieces, bbox$convention)
}

# Combine pieces covering different longitudes into one grid in the output convention
join_lon <- function(pieces, convention){
      lat <- pieces[[1]]$lat
      time <- pieces[[1]]$time
      for(p in pieces[-1]){
            if(!isTRUE(all.equal(p$lat, lat)) || !identical(p$time, time))
                  stop("data pieces on either side of the 0/360 meridian do not match", call. = FALSE)
      }
      lon <- convert_lon(unlist(lapply(pieces, `[[`, "lon")), convention)
      keep <- !duplicated(round(lon, 6))
      o <- order(lon[keep])
      vars <- names(pieces[[1]]$data)
      data <- lapply(stats::setNames(vars, vars), function(v){
            a <- do.call(abind_lon, lapply(pieces, function(p) p$data[[v]]))
            a[which(keep)[o], , , drop = FALSE]
      })
      list(lon = lon[keep][o], lat = lat, time = time, data = data)
}

abind_lon <- function(...){
      a <- list(...)
      d <- dim(a[[1]])
      out <- array(NA_real_, c(sum(vapply(a, function(z) dim(z)[1], numeric(1))), d[2], d[3]))
      i <- 0
      for(z in a){
            out[i + seq_len(dim(z)[1]), , ] <- z
            i <- i + dim(z)[1]
      }
      out
}

# Combine grids holding different variables (e.g. ERA5's separate u and v files)
merge_vars <- function(...){
      g <- list(...)
      for(h in g[-1]){
            if(!isTRUE(all.equal(h$lon, g[[1]]$lon)) || !isTRUE(all.equal(h$lat, g[[1]]$lat)))
                  stop("u and v data are on different grids", call. = FALSE)
            if(!identical(h$time, g[[1]]$time))
                  stop("u and v data have different time steps", call. = FALSE)
      }
      g[[1]]$data <- do.call(c, lapply(g, `[[`, "data"))
      g[[1]]
}

# Convert a grid to a SpatRaster with layers for each variable in `vars`, then each time step.
# Latitudes needn't be evenly spaced (CFSR uses a Gaussian grid): rows are spread evenly between
# the first and last latitude, as in the data set's nominal grid.
grid_to_rast <- function(g, vars){
      nx <- length(g$lon)
      ny <- length(g$lat)
      if(nx < 2 || ny < 2) stop("the bounding box must include at least 2 grid points in each ",
                                "direction; try a larger box", call. = FALSE)
      dx <- (max(g$lon) - min(g$lon)) / (nx - 1)
      dy <- (max(g$lat) - min(g$lat)) / (ny - 1)
      e <- terra::ext(min(g$lon) - dx / 2, max(g$lon) + dx / 2,
                      min(g$lat) - dy / 2, max(g$lat) + dy / 2)
      rows <- order(g$lat, decreasing = TRUE)
      layers <- lapply(vars, function(v){
            a <- aperm(g$data[[v]][, rows, , drop = FALSE], c(2, 1, 3)) # [lat, lon, time]
            r <- terra::rast(a, extent = e, crs = "EPSG:4326")
            names(r) <- if(length(g$time) == dim(a)[3] && !all(is.na(g$time)))
                  paste(v, format(g$time, "%Y-%m-%d %H:%M:%S", tz = "UTC")) else
                        paste(v, seq_len(dim(a)[3]))
            r
      })
      do.call(c, layers)
}


# Server access ------------------------------------------------------------------------------------

# Request a spatial subset of a file through the NetCDF Subset Service, trying the compressed
# netCDF-4 format first and falling back to netCDF-3.
ncss_get <- function(path, vars, west, east, south, north, time_stride, dest){
      query <- c(paste0("var=", vars),
                 paste0("north=", fmt_num(north)), paste0("south=", fmt_num(south)),
                 paste0("west=", fmt_num(west)), paste0("east=", fmt_num(east)),
                 "temporal=all", if(time_stride > 1) paste0("timeStride=", time_stride))
      url <- paste0(ncar_thredds(), "/ncss/grid/", path, "?",
                    paste(utils::URLencode(query, reserved = FALSE), collapse = "&"))
      err <- NULL
      for(fmt in c("netcdf4", "netcdf")){
            ok <- tryCatch({
                  ncss_fetch(paste0(url, "&accept=", fmt), dest)
                  is_netcdf(dest)
            }, error = function(e){
                  err <<- conditionMessage(e)
                  FALSE
            })
            if(ok) return(invisible(dest))
      }
      stop("failed to download ", path, " from NCAR", if(!is.null(err)) paste0(": ", err),
           "\nRequest URL: ", url, call. = FALSE)
}

# Download one URL to a file, with retries. (Separate from ncss_get() so tests can mock it.)
ncss_fetch <- function(url, dest, tries = 3){
      op <- options(timeout = max(3600, getOption("timeout")))
      on.exit(options(op))
      for(i in seq_len(tries)){
            res <- tryCatch(utils::download.file(url, dest, mode = "wb", quiet = TRUE),
                            error = function(e) e, warning = function(w) w)
            if(identical(res, 0L)) return(invisible(dest))
            if(i < tries) Sys.sleep(2 ^ i)
      }
      unlink(dest)
      stop(if(inherits(res, "condition")) conditionMessage(res) else "download failed", call. = FALSE)
}

# A server error page can arrive with a success code, so check the file's signature
is_netcdf <- function(file){
      if(!file.exists(file) || file.size(file) < 4) return(FALSE)
      sig <- readBin(file, "raw", 4)
      identical(sig[1:3], charToRaw("CDF")) || identical(sig, as.raw(c(0x89, 0x48, 0x44, 0x46)))
}


# Reading NetCDF ----------------------------------------------------------------------------------

# Read variables from a NetCDF file as arrays ordered [lon, lat, time], identifying dimensions
# by their units and names rather than assuming a layout. Extra dimensions (such as a single
# height level) must have length 1.
read_ncss <- function(file, vars){
      nc <- ncdf4::nc_open(file)
      on.exit(ncdf4::nc_close(nc))
      out <- list(data = list())
      for(k in names(vars)){
            v <- nc$var[[vars[[k]]]]
            if(is.null(v)) stop("variable '", vars[[k]], "' not found in downloaded file; found: ",
                                paste(names(nc$var), collapse = ", "), call. = FALSE)
            dims <- v$dim
            role <- vapply(dims, dim_role, character(1))
            if(sum(role == "lon") != 1 || sum(role == "lat") != 1)
                  stop("could not identify longitude and latitude dimensions of '", vars[[k]], "'",
                       call. = FALSE)
            len <- vapply(dims, function(d) d$len, numeric(1))
            if(any(role == "other" & len > 1))
                  stop("'", vars[[k]], "' has an unexpected extra dimension", call. = FALSE)
            a <- ncdf4::ncvar_get(nc, v, collapse_degen = FALSE)
            a <- array(a, dim = len)
            if(!any(role == "time")){ # a single time, or no time dimension
                  a <- array(a, dim = c(len, 1))
                  role <- c(role, "time")
                  tvals <- NA
            }else{
                  td <- dims[[which(role == "time")]]
                  tvals <- parse_nc_time(td$vals, td$units)
            }
            keep <- match(c("lon", "lat", "time"), role)
            a <- aperm(a, c(keep, setdiff(seq_along(role), keep)))
            a <- array(a, dim = dim(a)[1:3])
            lon <- dims[[keep[1]]]$vals
            lat <- dims[[keep[2]]]$vals
            if(is.null(out$lon)){
                  out$lon <- lon
                  out$lat <- lat
                  out$time <- tvals
            }else if(!isTRUE(all.equal(out$lon, lon)) || !isTRUE(all.equal(out$lat, lat)) ||
                     !identical(out$time, tvals)){
                  stop("variables in downloaded file are on different grids", call. = FALSE)
            }
            out$data[[k]] <- a
      }
      # ncss_get doesn't guarantee coordinate order; standardize to increasing longitude
      o <- order(out$lon)
      out$lon <- out$lon[o]
      out$data <- lapply(out$data, function(a) a[o, , , drop = FALSE])
      out
}

dim_role <- function(d){
      nm <- tolower(d$name)
      un <- tolower(if(is.null(d$units)) "" else d$units)
      if(grepl("^degrees?_?e(ast)?$", un) || nm %in% c("lon", "longitude", "x")) return("lon")
      if(grepl("^degrees?_?n(orth)?$", un) || nm %in% c("lat", "latitude", "y")) return("lat")
      if(grepl(" since ", un)) return("time")
      "other"
}

# Convert NetCDF time values to POSIXct, from units like "hours since 1900-01-01 00:00:00" or
# "Hour since 2000-04-01T00:00:00Z"
parse_nc_time <- function(vals, units){
      parts <- strsplit(units, " since ", fixed = TRUE)[[1]]
      if(length(parts) != 2) stop("unrecognized time units: ", units, call. = FALSE)
      unit <- sub("s$", "", tolower(trimws(parts[1])))
      mult <- unname(c(second = 1, sec = 1, minute = 60, min = 60, hour = 3600, hr = 3600, h = 3600,
                day = 86400, d = 86400)[unit])
      if(is.na(mult)) stop("unrecognized time units: ", units, call. = FALSE)
      origin_str <- sub("Z$", "", sub("T", " ", trimws(parts[2])))
      origin_str <- sub(" UTC$", "", origin_str)
      origin <- NA
      for(f in c("%Y-%m-%d %H:%M:%OS", "%Y-%m-%d %H:%M", "%Y-%m-%d")){
            origin <- as.POSIXct(origin_str, format = f, tz = "UTC")
            if(!is.na(origin)) break
      }
      if(is.na(origin)) stop("unrecognized time units: ", units, call. = FALSE)
      # round to the second, avoiding floating-point drift in layer names
      as.POSIXct(round(as.numeric(origin) + as.numeric(vals) * mult), origin = "1970-01-01", tz = "UTC")
}
