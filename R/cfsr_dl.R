cfsr_dl_day <- function(year = 1979, month = 1, day = 1,
                        xmin = 260, xmax = 270, ymin = 40, ymax = 50){

      if(! year %in% 1979:2010) stop("'year' must be between 1979 and 2010")
      if(! month %in% 1:12) stop("'month' must be between 1 and 12")
      if(min(xmin, xmax) < 0 | max(xmin, xmax) > 360) stop("longitude and latitude must be between 0 and 360")

      require(ncdf4)
      require(lubridate)
      require(terra)

      # open connection
      # example url: "https://thredds.rda.ucar.edu/thredds/dodsC/files/g/ds093.1/1990/wnd10m.gdas.199001.grb2"
      url <- paste0("https://thredds.rda.ucar.edu/thredds/dodsC/files/g/ds093.1/", year, "/wnd10m.gdas.",
                    year, stringr::str_pad(month, 2, "left", 0), ".grb2")
      ds <- nc_open(url)

      # dimensions
      lon <- ncvar_get(ds, "lon")
      lat <- ncvar_get(ds, "lat")
      time <- ncvar_get(ds, "time")

      # convert time
      t_units <- ncatt_get(ds, "time", "units")
      startdate <- gsub("T00:00:00Z", "", unlist(strsplit(t_units$value, " "))[3])
      timestamp <- ymd(startdate) + dhours(time-1)

      # bounds
      x <- c(xmin, xmax)
      # x[x<0] <- x[x<0] + 360
      y <- c(ymin, ymax)
      date <- paste0(year, "-",
                     stringr::str_pad(month, 2, "left", 0), "-",
                     stringr::str_pad(day, 2, "left", 0))
      t <- as.POSIXct(paste0(date, c(" 00:00:00 UTC", " 23:00:00 UTC")), tz = "UTC")

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
      v <- ncvar_get(ds, "v-component_of_wind_height_above_ground", start = start, count = count)
      u <- ncvar_get(ds, "u-component_of_wind_height_above_ground", start = start, count = count)

      # convert to raster object
      xres <- base::diff(sort(x)) / (lon_count - 1)
      yres <- base::diff(sort(y)) / (lat_count - 1)
      extent <- c(x[1] - 0.5 * xres, x[2] + 0.5 * xres,
                  y[1] - 0.5 * yres, y[2] + 0.5 * yres)
      vr <- rast(aperm(v, c(2, 1, 3)), extent = extent)
      ur <- rast(aperm(u, c(2, 1, 3)), extent = extent)
      names(vr) <- paste("v", timestamp[time_i[1]:time_i[2]])
      names(ur) <- paste("u", timestamp[time_i[1]:time_i[2]])
      return(list(u = ur, v = vr))
}



#' Download wind data from the Climate Forecast System Reanalysis (CFSR)
#'
#' This function downloads hourly CFSR wind data from the NCAR THREDDS server, for one or more historic dates, within a spatial bounding box.
#' Data for all factorial combinations of 'years', 'months', and 'days' are downloaded and combined.
#'
#' @param years Integer vector representing the historic year(s) for which data are to be accessed, in the range 1979-2010.
#' @param months Integer vector representing the month(s) for which data are to be accessed, in the range 1-12.
#' @param days Integer vector representing day(s) for which data is to be requested, in the range 0-31.
#' @param xlim Vector of length 2, giving longitudinal limits of the data to download, in the range 0 to 360.
#' @param ylim Vector of length 2, giving latitudinal limits of the data to download, in the range -90 to 90.
#' @return A list containing u and v components of the wind field, each as a SpatRaster with a layer for each hour of each requested date.
#' @export
cfsr_dl <- function(years = 1979, months = 1, days = 1,
                    xlim = c(260, 270), ylim = c(40, 50)){
      require(terra)
      require(purrr)
      q <- expand.grid(day = days, month = months, year = years,
                       xmin = min(xlim), xmax = max(xlim), ymin = min(ylim), ymax = max(ylim))
      w <- pmap(q, cfsr_dl_day)
      c(rast(map(w, "u")), rast(map(w, "v")))
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

      require(ncdf4)
      require(terra)

      # open connection
      url <- paste0("https://thredds.rda.ucar.edu/thredds/dodsC/files/g/ds093.1/1980/soilt1.gdas.198001.grb2")
      ds <- nc_open(url)

      # dimensions
      lon <- ncvar_get(ds, "lon")
      lat <- ncvar_get(ds, "lat")

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
      v <- ncvar_get(ds, "Temperature_depth_below_surface_layer", start = start, count = count)

      # convert to raster object
      xres <- base::diff(sort(x)) / (lon_count - 1)
      yres <- base::diff(sort(y)) / (lat_count - 1)
      extent <- c(x[1] - 0.5 * xres, x[2] + 0.5 * xres,
                  y[1] - 0.5 * yres, y[2] + 0.5 * yres)
      vr <- rast(aperm(v, c(2, 1, 3)[1:2]), extent = extent)

      # classify
      vr[is.finite(vr)] <- 1
      vr[is.na(vr)] <- 0
      names(vr) <- "land"
      return(vr)
}
