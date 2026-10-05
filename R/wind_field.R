#' An S4 object class representing a wind field time series
#'
setClass("wind_series",
         contains = "SpatRaster",
         slots = c(n_steps = "numeric"))

#' Generate a wind field time series data set from a set of rasters.
#'
#' @param x Multi-layer \code{SpatRaster} with layers containing u and v wind components, or
#'    an object like a file path that can be converted to a \code{SpatRaster}. Data must be
#'    on a longitude/latitude grid with square cells (equal x and y resolution in degrees),
#'    with u and v components oriented to true east and north; see \link{wind_rose}.
#' @param order Either \code{"uuvv"}, the default, indicating `x` has all u components
#'    followed by all v components, or \code{"uvuv"}, indicating the u and v components
#'    of `x` are alternating.
#' @return A `wind_series` object, which is a particular form of \code{SpatRaster}.
#' @export
wind_series <- function(x, order = c("uuvv", "uvuv")){

      order <- match.arg(order)
      if(!inherits(x, "SpatRaster")) x <- rast(x)
      check_grid(x)
      if(terra::nlyr(x) %% 2 != 0) stop("if `v` is not specified, `u` must have an even number of layers.")

      # collate layers
      if(order == "uvuv"){
            even <- function(x) x %% 2 == 0
            x <- x[[c(which(!even(1:nlyr(x))),
                      which(even(1:nlyr(x))))]]
      }

      # create wind_field
      y <- as(as(x, "SpatRaster"), "wind_series")
      y@n_steps <- nlyr(x)/2
      y
}



#' An S4 object class representing a wind field
#'
setClass("wind_field",
         contains = "SpatRaster",
         slots = c(n_steps = "numeric"))


#' Generate a wind field data set from a pair of rasters.
#'
#' @param x A `SpatRaster` with two layers representing u and v wind components;
#'    note that these must be in lat-long coordinates.
#' @return A `wind_field` object, which is a particular form of `SpatRaster`.
#' @export
wind_field <- function(x){
      if(terra::nlyr(x) != 2) stop("`x` must have two layers.")
      xt <- ext(x)
      if(any(c(xt$xmin < -180, xt$xmax > 360, xt$ymin < -90, xt$ymax > 90))) stop("wind field rasters must be in lon-lat coordinates")
      as(as(x, "SpatRaster"), "wind_field") # via SpatRaster, so subclasses (e.g. wind_series layers) work
}


#' Load a wind_series from one or more raster files on disk
#'
#' Reads files in `wind_series` layout (all u layers followed by all v layers), such as the
#' monthly files saved by [ncar_download()], and combines them into a single `wind_series`.
#' Data are not loaded into memory until needed.
#'
#' @param x Character vector of file paths. Each file must hold an equal number of u and v
#'   layers, with u layers first, and all files must share the same grid. Time steps are combined
#'   in the order the files are given.
#' @return A `wind_series` object whose time steps are those of all the files combined: all
#'   files' u layers, followed by all files' v layers.
#' @seealso [ncar_download()]
#' @export
read_wind_series <- function(x){
      if(!is.character(x) || length(x) == 0) stop("`x` must be a character vector of file paths")
      missing <- !file.exists(x)
      if(any(missing)) stop("file(s) not found: ", paste(x[missing], collapse = ", "))
      r <- lapply(x, terra::rast)
      for(i in seq_along(r)){
            if(terra::nlyr(r[[i]]) %% 2 != 0) stop("file has an odd number of layers: ", x[i])
            if(i > 1 && !terra::compareGeom(r[[1]], r[[i]], stopOnError = FALSE))
                  stop("files are on different grids: ", x[1], " and ", x[i])
      }
      half <- function(z, which){
            n <- terra::nlyr(z) / 2
            z[[if(which == "u") seq_len(n) else n + seq_len(n)]]
      }
      wind_series(c(do.call(c, lapply(r, half, "u")), do.call(c, lapply(r, half, "v"))))
}


#' Get the time of each step in a wind_series
#'
#' Returns the time of each step of a `wind_series`, parsed from its layer names (like
#' `"u 2000-01-01 06:00:00"`, as written by [ncar_download()]) or, if the names hold no times,
#' from the times set with `terra::time()`.
#'
#' @param x A `wind_series`.
#' @return A POSIXct vector (UTC), with one element per time step.
#' @seealso [subset_series()] to select time steps.
#' @examples
#' series <- windscape_example("wind_series")
#' head(wind_times(series))
#' @export
wind_times <- function(x){
      if(!inherits(x, "wind_series")) stop("`x` must be a `wind_series`")
      u <- seq_len(x@n_steps)
      t <- layer_times(names(x)[u])
      if(all(is.na(t))){
            tt <- terra::time(x)[u]
            if(!all(is.na(tt))) t <- as.POSIXct(tt, tz = "UTC")
      }
      if(all(is.na(t))) stop("no times found: layer names should be like \"u 2000-01-01 06:00:00\", ",
                             "or times should be set with terra::time()")
      if(anyNA(t)) warning("times could not be determined for ", sum(is.na(t)), " time step(s)")
      t
}


#' Select time steps from a wind_series
#'
#' Selects time steps from a `wind_series` by month, hour of day, date range, or index, keeping
#' the u and v layers of each step together. Use it to build wind roses for the times that
#' matter for dispersal, such as a flowering season or the hours when propagules are released.
#'
#' Criteria are combined: a time step is kept only if it meets all of them. Times are in UTC.
#' For hours of the day in local solar time, note that local time is about UTC plus longitude /
#' 15 hours (e.g. UTC - 7 hours at 105 degrees west), so `hours` for a large region selects
#' different local times in different places.
#'
#' @param x A `wind_series`, with times in its layer names (see [wind_times()]) if selecting by
#'    `months`, `hours`, `start`, or `end`.
#' @param months Integer vector of months (1-12) to keep.
#' @param hours Integer vector of hours of the day (0-23, UTC) to keep.
#' @param start,end Keep time steps from `start` to `end`, inclusive. Each can be a `Date`, a
#'    POSIXct time, or a character string like `"2000-06-01"` or `"2000-06-01 12:00:00"`
#'    (interpreted as UTC). A date without a time includes that whole day, so
#'    `end = "2000-06-30"` keeps all of June 30.
#' @param steps Time steps to keep, as integer indices or a logical vector with one element per
#'    time step, for selections the other arguments don't cover. Can be combined with them, e.g.
#'    `steps = speed > 10`, given a vector `speed`.
#' @return A `wind_series` containing the selected time steps, in their original order.
#' @seealso [wind_times()]
#' @examples
#' series <- windscape_example("wind_series")
#'
#' # summer only
#' summer <- subset_series(series, months = 6:8)
#' range(wind_times(summer))
#'
#' # afternoons in the first half of the year
#' subset_series(series, hours = 18:23, end = "2000-06-30")
#' @export
subset_series <- function(x, months = NULL, hours = NULL, start = NULL, end = NULL, steps = NULL){
      if(!inherits(x, "wind_series")) stop("`x` must be a `wind_series`")
      keep <- keep_steps(x, months, hours, start, end, steps)
      idx <- which(keep)
      if(length(idx) == 0) stop("no time steps meet the selection criteria")
      n <- x@n_steps
      wind_series(terra::subset(as(x, "SpatRaster"), c(idx, n + idx)))
}

# Logical vector: which time steps of a wind_series meet the selection criteria
keep_steps <- function(x, months = NULL, hours = NULL, start = NULL, end = NULL, steps = NULL){
      n <- x@n_steps
      keep <- rep(TRUE, n)

      if(!is.null(steps)){
            if(is.logical(steps)){
                  if(length(steps) != n) stop("a logical `steps` must have one element per time step (", n, ")")
                  keep <- keep & !is.na(steps) & steps
            }else{
                  if(!is.numeric(steps) || any(steps %% 1 != 0) || any(steps < 1 | steps > n))
                        stop("`steps` must be logical, or integers between 1 and ", n)
                  keep <- keep & seq_len(n) %in% steps
            }
      }

      if(!is.null(months) || !is.null(hours) || !is.null(start) || !is.null(end)){
            t <- wind_times(x)
            if(!is.null(months)){
                  if(any(!months %in% 1:12)) stop("`months` must be integers between 1 and 12")
                  keep <- keep & as.integer(format(t, "%m", tz = "UTC")) %in% months
            }
            if(!is.null(hours)){
                  if(any(!hours %in% 0:23)) stop("`hours` must be integers between 0 and 23")
                  keep <- keep & as.integer(format(t, "%H", tz = "UTC")) %in% hours
            }
            if(!is.null(start)) keep <- keep & t >= as_utc_time(start, "start")
            if(!is.null(end)){
                  e <- as_utc_time(end, "end")
                  keep <- keep & if(attr(e, "whole_day")) t < e + 86400 else t <= e
            }
            keep <- keep & !is.na(t)
      }
      keep
}

# Convert a Date, POSIXct, or character time to POSIXct (UTC), recording whether it was a date
# without a time of day
as_utc_time <- function(z, arg){
      if(length(z) != 1 || is.na(z)) stop("`", arg, "` must be a single date or time")
      whole_day <- FALSE
      if(inherits(z, "POSIXt")){
            out <- as.POSIXct(z, tz = "UTC")
      }else if(inherits(z, "Date")){
            out <- as.POSIXct(format(z), tz = "UTC")
            whole_day <- TRUE
      }else if(is.character(z)){
            whole_day <- grepl("^\\s*\\d{4}-\\d{2}-\\d{2}\\s*$", z)
            out <- as.POSIXct(z, tz = "UTC", tryFormats = c("%Y-%m-%d %H:%M:%OS", "%Y-%m-%d %H:%M",
                                                           "%Y-%m-%d"), optional = TRUE)
      }else{
            stop("`", arg, "` must be a Date, POSIXct, or character string")
      }
      if(is.na(out)) stop("`", arg, "` could not be read as a date or time: ", z)
      attr(out, "whole_day") <- whole_day
      out
}
