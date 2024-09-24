

cfsr_dl_block <- function(variable = "wnd10m",
                        year = 1979, month = 1, day = 1, hmin = 0, hmax = 23,
                        xmin = 260, xmax = 270, ymin = 40, ymax = 50){

      if(! year %in% 1979:2010) stop("'year' must be between 1979 and 2010")
      if(! month %in% 1:12) stop("'month' must be between 1 and 12")
      if(min(xmin, xmax) < 0 | max(xmin, xmax) > 360) stop("longitude and latitude must be between 0 and 360")

      # open connection
      # example url: "https://thredds.rda.ucar.edu/thredds/dodsC/files/g/ds093.1/1990/wnd10m.gdas.199001.grb2"
      url <- paste0("https://thredds.rda.ucar.edu/thredds/dodsC/files/g/ds093.1/", year, "/", variable, ".gdas.",
                    year, stringr::str_pad(month, 2, "left", 0), ".grb2")
      ds <- ncdf4::nc_open(url)

      # dimensions
      lon <- ncdf4::ncvar_get(ds, "lon")
      lat <- ncdf4::ncvar_get(ds, "lat")
      time <- ncdf4::ncvar_get(ds, "time")

      # convert time
      t_units <- ncdf4::ncatt_get(ds, "time", "units")
      startdate <- gsub("T00:00:00Z", "", unlist(strsplit(t_units$value, " "))[3])
      timestamp <- lubridate::ymd(startdate) + lubridate::dhours(time-1)

      # bounds
      x <- c(xmin, xmax)
      y <- c(ymin, ymax)
      date <- paste0(year, "-",
                     stringr::str_pad(month, 2, "left", 0), "-",
                     stringr::str_pad(day, 2, "left", 0))
      t <- as.POSIXct(paste0(date, c(paste0(" ", stringr::str_pad(hmin, 2, "left", 0), ":00:00 UTC"),
                                     paste0(" ", stringr::str_pad(hmax, 2, "left", 0), ":00:00 UTC"))), tz = "UTC")

            # indices
      btw <- function(data, z) range(which(data <= max(z) & data >= min(z)))
      lon_i <- btw(lon, x)
      lat_i <- btw(lat, y)
      time_i <- btw(timestamp, t)
      lon_count <- lon_i[-1] - lon_i[1] + 1
      lat_count <- lat_i[-1] - lat_i[1] + 1
      time_count <- time_i[-1] - time_i[1] + 1
      start <- c(lon_i[1], lat_i[1], 1, time_i[1]) # x,y,z,t
      count <- c(lon_count, lat_count, 1, time_count)

      # get data
      message(paste0("... retrieving CFSR data for ", date, " ..."))
      tag <- ifelse(variable == "wnd10m", "height_above_ground", "isobaric")
      v <- ncdf4::ncvar_get(ds, paste0("v-component_of_wind_", tag),
                            start = start, count = count)
      u <- ncdf4::ncvar_get(ds, paste0("u-component_of_wind_", tag),
                            start = start, count = count)

      # convert to raster object
      xres <- base::diff(sort(x)) / (lon_count - 1)
      yres <- base::diff(sort(y)) / (lat_count - 1)
      extent <- c(x[1] - 0.5 * xres, x[2] + 0.5 * xres,
                  y[1] - 0.5 * yres, y[2] + 0.5 * yres)
      nd <- 1:length(dim(v)) # 1:2 if only a single hour, else 1:3
      vr <- terra::rast(aperm(v, c(2, 1, 3)[nd]), extent = extent)
      ur <- terra::rast(aperm(u, c(2, 1, 3)[nd]), extent = extent)
      names(vr) <- paste("v", timestamp[time_i[1]:time_i[2]])
      names(ur) <- paste("u", timestamp[time_i[1]:time_i[2]])
      return(list(u = ur, v = vr))
}



#' Download wind data from the Climate Forecast System Reanalysis (CFSR)
#'
#' This function downloads hourly CFSR wind data from the NCAR THREDDS server, for one or more historic dates, within a spatial bounding box.
#' Data for all factorial combinations of 'years', 'months', 'days', and 'hours' are downloaded and combined.
#'
#' @param variable Name of the CFSR wind variable to download. Options include "wnd10m" (the default, representing wind speed
#'    10 m above the ground), as well as "wnd1000", "wnd850", "wnd700", "wnd500", and "wnd200" (representing higher altitudes,
#'    with the numbers indicating pressure levels in hPa).
#' @param years Integer vector representing the historic year(s) for which data are to be accessed, in the range 1979-2010.
#' @param months Integer vector representing the month(s) for which data are to be accessed, in the range 1-12.
#' @param days Integer vector representing day(s) for which data is to be requested, in the range 0-31.
#' @param hlim Vector of length 2, containing integers between 0 and 23, giving the range of hours (in UTC time zone) to download.
#' @param xlim Vector of length 2, giving longitudinal limits of the data to download, in the range 0 to 360.
#' @param ylim Vector of length 2, giving latitudinal limits of the data to download, in the range -90 to 90.
#' @return A list containing u and v components of the wind field, each as a SpatRaster with a layer for each hour of each requested date.
#' @export
cfsr_dl <- function(variable = "wnd10m",
                    years = 1979, months = 1, days = 1,
                    hlim = c(0, 23), xlim = c(260, 270), ylim = c(40, 50)){
      q <- expand.grid(variable = variable, day = days, month = months, year = years, hmin = min(hlim), hmax = max(hlim),
                       xmin = min(xlim), xmax = max(xlim), ymin = min(ylim), ymax = max(ylim))
      w <- purrr::pmap(q, cfsr_dl_block)
      c(terra::rast(purrr::map(w, "u")),
        terra::rast(purrr::map(w, "v")))
}



#' Download land-water raster from the Climate Forecast System Reanalysis (CFSR)
#'
#' This function downloads a CFSR land-water raster layer from the NCAR THREDDS server, for an area specified by a spatial bounding box.
#'
#' (Actually, since CFSR does not publish a land layer per se, this function downloads a soil temperature layer for an arbitrary date and classifies areas with data as land and no data as water.)
#'
#' @param xlim Vector of length 2, giving longitudinal limits of the data to download, in the range 0 to 360.
#' @param ylim Vector of length 2, giving latitudinal limits of the data to download, in the range -90 to 90.
#' @return A SpatRaster layer with 1 indicating land and 0 indicating water.
#' @export
cfsr_dl_land <- function(xlim = c(260, 270), ylim = c(40, 50)){

      # open connection
      require(ncdf4)
      url <- paste0("https://thredds.rda.ucar.edu/thredds/dodsC/files/g/ds093.1/1980/soilt1.gdas.198001.grb2")
      ds <- nc_open(url)

      # dimensions
      lon <- ncdf4::ncvar_get(ds, "lon")
      lat <- ncdf4::ncvar_get(ds, "lat")

      # bounds
      x <- sort(xlim)
      y <- sort(ylim)

      # indices
      btw <- function(data, z) range(which(data <= max(z) & data >= min(z)))
      lon_i <- btw(lon, x)
      lat_i <- btw(lat, y)
      lon_count <- lon_i[-1] - lon_i[1] + 1
      lat_count <- lat_i[-1] - lat_i[1] + 1
      start <- c(lon_i[1], lat_i[1], 1, 1) # x,y,z,t
      count <- c(lon_count, lat_count, 1, 1)

      # get data (land temperature data used as land indicator)
      message("... retrieving CFSR land layer ...")
      v <- ncdf4::ncvar_get(ds, "Temperature_depth_below_surface_layer",
                            start = start, count = count)

      # convert to raster object
      xres <- base::diff(sort(x)) / (lon_count - 1)
      yres <- base::diff(sort(y)) / (lat_count - 1)
      extent <- c(x[1] - 0.5 * xres, x[2] + 0.5 * xres,
                  y[1] - 0.5 * yres, y[2] + 0.5 * yres)
      vr <- terra::rast(aperm(v, c(2, 1, 3)[1:2]), extent = extent)

      # classify
      vr[is.finite(vr)] <- 1
      vr[is.na(vr)] <- 0
      names(vr) <- "land"
      return(vr)
}
