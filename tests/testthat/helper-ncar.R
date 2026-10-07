# A fake NCAR NetCDF Subset Service for testing download_wind_data() offline. fake_ncss() returns a
# replacement for ncss_fetch(url, dest) that parses the request URL and writes a NetCDF file
# laid out like the real data set (ERA5, CFSR/CFSv2, or the ERA5 land mask), on a coarse 2 degree
# global grid with 8 time steps per month. Values come from truth_u() and truth_v(), so tests can
# check them at known locations and times.

truth_u <- function(lon, lat, time) 10 * sin(lon * pi / 180) + lat / 10 + as.numeric(format(time, "%H", tz = "UTC")) / 10
truth_v <- function(lon, lat, time) 5 * cos(lon * pi / 180) - lat / 20 + as.numeric(format(time, "%d", tz = "UTC")) / 10

parse_ncss_url <- function(url){
      path <- sub("\\?.*$", "", sub("^.*/ncss/grid/", "", url))
      q <- strsplit(sub("^[^?]*\\?", "", url), "&")[[1]]
      key <- sub("=.*$", "", q)
      val <- utils::URLdecode(sub("^[^=]*=", "", q))
      p <- split(val, key)
      list(path = path, var = p$var, west = as.numeric(p$west), east = as.numeric(p$east),
           south = as.numeric(p$south), north = as.numeric(p$north),
           stride = if(is.null(p$timeStride)) 1 else as.integer(p$timeStride),
           accept = p$accept, temporal = p$temporal)
}

fake_ncss <- function(log = new.env(), reject_netcdf4 = FALSE, html = FALSE){
      log$urls <- character()
      function(url, dest, tries = 3){
            log$urls <- c(log$urls, url)
            q <- parse_ncss_url(url)
            if(reject_netcdf4 && q$accept == "netcdf4") stop("HTTP status 400")
            if(html){
                  writeLines("<html>error</html>", dest)
                  return(invisible(dest))
            }
            stopifnot(identical(q$temporal, "all"))
            lon <- seq(0, 358, 2)
            lon <- lon[lon >= q$west & lon <= q$east]
            lat_all <- seq(90, -90, -2)
            lat <- lat_all[lat_all >= q$south & lat_all <= q$north]
            v4 <- q$accept == "netcdf4"

            if(grepl("invariant", q$path)){ # ERA5 land-sea mask
                  t0 <- as.POSIXct("1979-01-01", tz = "UTC")
                  dims <- list(ncdf4::ncdim_def("longitude", "degrees_east", lon),
                               ncdf4::ncdim_def("latitude", "degrees_north", lat),
                               ncdf4::ncdim_def("time", "hours since 1900-01-01 00:00:00",
                                                as.integer(difftime(t0, as.POSIXct("1900-01-01", tz = "UTC"), units = "hours"))))
                  vd <- ncdf4::ncvar_def("LSM", "1", dims, prec = "float")
                  nc <- ncdf4::nc_create(dest, vd, force_v4 = v4)
                  ncdf4::ncvar_put(nc, vd, outer(lon, lat, function(x, y) as.numeric(x > 180) * 0.75))
                  ncdf4::nc_close(nc)
                  return(invisible(dest))
            }

            ym <- regmatches(q$path, regexpr("\\d{6}(?=(0100_|\\.grb2))", q$path, perl = TRUE))
            month_start <- as.POSIXct(paste0(substr(ym, 1, 4), "-", substr(ym, 5, 6), "-01"), tz = "UTC")

            if(grepl("d633000", q$path)){ # ERA5: one variable per file, int hours since 1900
                  times <- month_start + 3600 * (0:7)
                  times <- times[seq(1, length(times), q$stride)]
                  tv <- as.integer(round(as.numeric(difftime(times, as.POSIXct("1900-01-01", tz = "UTC"), units = "hours"))))
                  dims <- list(ncdf4::ncdim_def("longitude", "degrees_east", lon),
                               ncdf4::ncdim_def("latitude", "degrees_north", lat),
                               ncdf4::ncdim_def("time", "hours since 1900-01-01 00:00:00", tv))
                  fun <- if(grepl("U$", q$var)) truth_u else truth_v
                  vd <- ncdf4::ncvar_def(q$var, "m s**-1", dims, prec = "float")
                  a <- array(NA_real_, c(length(lon), length(lat), length(times)))
                  for(k in seq_along(times)) a[, , k] <- outer(lon, lat, fun, time = times[k])
                  nc <- ncdf4::nc_create(dest, vd, force_v4 = v4)
                  ncdf4::ncvar_put(nc, vd, a)
                  ncdf4::nc_close(nc)
                  return(invisible(dest))
            }

            # CFSR/CFSv2: u and v in one file, a degenerate height dimension, ascending latitudes,
            # and times starting at hour 1 in units like "Hour since 2000-04-01T00:00:00Z"
            lat <- rev(lat)
            hv <- 1:8
            hv <- hv[seq(1, length(hv), q$stride)]
            times <- month_start + 3600 * hv
            tunits <- paste0("Hour since ", format(month_start, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
            dims <- list(ncdf4::ncdim_def("lon", "degrees_east", lon),
                         ncdf4::ncdim_def("lat", "degrees_north", lat),
                         ncdf4::ncdim_def("height_above_ground", "m", 10),
                         ncdf4::ncdim_def("time", tunits, hv))
            vds <- lapply(q$var, ncdf4::ncvar_def, units = "m/s", dim = dims, prec = "float")
            nc <- ncdf4::nc_create(dest, vds, force_v4 = v4)
            for(i in seq_along(q$var)){
                  fun <- if(grepl("^u", q$var[i])) truth_u else truth_v
                  a <- array(NA_real_, c(length(lon), length(lat), 1, length(times)))
                  for(k in seq_along(times)) a[, , 1, k] <- outer(lon, lat, fun, time = times[k])
                  ncdf4::ncvar_put(nc, vds[[i]], a)
            }
            ncdf4::nc_close(nc)
            invisible(dest)
      }
}
