#' An S4 object class representing a wind field time series
#'
setClass("wind_series",
         contains = "SpatRaster",
         slots = c(n_steps = "numeric"))

#' Create a wind_series
#'
#' Creates a `wind_series`, a time series of wind fields, from a raster, a list of rasters, or
#' one or more raster files, such as the monthly files saved by [download_wind_data()]. Multiple
#' inputs are combined into one series, in the order given, without reading the data into
#' memory until needed.
#'
#' @param x A multi-layer `SpatRaster` with layers containing u and v wind components; a list of
#'    them; or a character vector of paths to raster files. Each input must contain an equal
#'    number of u and v layers, arranged as given by `order`, and all inputs must share the same
#'    grid. Data must be on a longitude/latitude grid with square cells (equal x and y resolution
#'    in degrees), with u and v components in m/s, oriented to true east and north; see
#'    [wind_rose()].
#' @param order Either `"uuvv"`, the default, indicating that each input has all u components
#'    followed by all v components, or `"uvuv"`, indicating alternating u and v components.
#' @return A `wind_series` object, a `SpatRaster` with all the u layers of all inputs, followed
#'    by all the v layers, in time order.
#' @seealso [subset_series()] and [wind_times()] for working with time steps.
#' @examples
#' series <- windscape_example("wind_series")
#' series
#'
#' \dontrun{
#' # combine monthly files from download_wind_data()
#' series <- wind_series(files)
#' }
#' @export
wind_series <- function(x, order = c("uuvv", "uvuv")){

      order <- match.arg(order)
      if(inherits(x, "SpatRaster")){
            inputs <- list(x)
      }else if(is.character(x)){
            if(length(x) == 0) stop("`x` must contain at least one file path")
            missing <- !file.exists(x)
            if(any(missing)) stop("file(s) not found: ", paste(x[missing], collapse = ", "))
            inputs <- lapply(x, terra::rast)
      }else if(is.list(x)){
            if(length(x) == 0 || !all(vapply(x, inherits, logical(1), "SpatRaster")))
                  stop("a list `x` must contain SpatRasters")
            inputs <- x
      }else{
            inputs <- list(terra::rast(x))
      }

      for(i in seq_along(inputs)){
            if(terra::nlyr(inputs[[i]]) %% 2 != 0)
                  stop("each input must have an even number of layers (u and v)",
                       if(length(inputs) > 1) paste0("; input ", i, " does not"))
            if(i > 1 && !terra::compareGeom(inputs[[1]], inputs[[i]], stopOnError = FALSE))
                  stop("inputs are on different grids (inputs 1 and ", i, ")")
      }
      check_grid(inputs[[1]])

      # split each input into u and v layers, then combine all u's followed by all v's
      half <- function(r, comp){
            r <- as(r, "SpatRaster")
            n <- terra::nlyr(r)
            idx <- if(order == "uuvv"){
                  if(comp == "u") seq_len(n / 2) else n / 2 + seq_len(n / 2)
            }else{
                  if(comp == "u") seq(1, n, 2) else seq(2, n, 2)
            }
            r[[idx]]
      }
      inputs <- unname(inputs) # named inputs would become argument names in c()
      x <- c(do.call(c, lapply(inputs, half, "u")), do.call(c, lapply(inputs, half, "v")))

      y <- as(x, "wind_series")
      y@n_steps <- terra::nlyr(x) / 2
      y
}



#' An S4 object class representing a wind field
#'
setClass("wind_field",
         contains = "SpatRaster",
         slots = c(n_steps = "numeric"))


#' Create a wind_field
#'
#' Creates a `wind_field`, the wind across a grid at a single moment, from a two-layer raster
#' or file, or from one time step of a `wind_series`. To summarize a whole series as a single
#' field, use [mean()][mean,wind_series-method] or [net_flow()].
#'
#' @param x A `SpatRaster` with two layers holding the u and v wind components, the path to a
#'    raster file with those two layers, or a `wind_series`. Data must be on a
#'    longitude/latitude grid, with components in m/s, oriented to true east and north.
#' @param step For a `wind_series` with more than one time step, the time step to use, as an
#'    integer index (see [wind_times()] for the time of each step).
#' @return A `wind_field` object, which is a particular form of `SpatRaster`.
#' @seealso [mean()][mean,wind_series-method] for the time-mean wind of a series.
#' @examples
#' series <- windscape_example("wind_series")
#' field <- wind_field(series, step = 1)
#' @export
wind_field <- function(x, step = NULL){
      if(is.character(x)){
            if(length(x) != 1) stop("`x` must be a single file path")
            if(!file.exists(x)) stop("file not found: ", x)
            x <- terra::rast(x)
      }
      if(!inherits(x, "SpatRaster")) stop("`x` must be a SpatRaster, a file path, or a wind_series")

      # layer subsetting (e.g. series[[5]]) keeps the wind_series class without updating n_steps,
      # so treat a series whose layers don't match its n_steps as a plain raster
      if(inherits(x, "wind_series") && terra::nlyr(x) != 2 * x@n_steps) x <- as(x, "SpatRaster")
      if(inherits(x, "wind_series")){
            n <- x@n_steps
            if(is.null(step)){
                  if(n > 1) stop("`x` has ", n, " time steps; choose one with `step`, ",
                                 "or summarize them with mean()")
                  step <- 1
            }
            if(length(step) != 1 || !is.numeric(step) || step %% 1 != 0 || step < 1 || step > n)
                  stop("`step` must be a single integer between 1 and ", n)
            x <- as(x, "SpatRaster")[[c(step, n + step)]]
      }else if(!is.null(step)){
            stop("`step` applies only to a wind_series")
      }

      if(terra::nlyr(x) != 2) stop("`x` must have two layers (u and v); for wind data with more ",
                                   "time steps, use wind_series()")
      xt <- ext(x)
      if(any(c(xt$xmin < -180, xt$xmax > 360, xt$ymin < -90, xt$ymax > 90))) stop("wind field rasters must be in lon-lat coordinates")
      as(as(x, "SpatRaster"), "wind_field") # via SpatRaster, so subclasses (e.g. wind_series layers) work
}


#' Mean wind of a wind_series
#'
#' Averages the u and v components of a `wind_series` over its time steps, giving the mean wind
#' vector in each cell as a `wind_field`. Long series stored in files are processed in blocks,
#' without loading the whole series into memory.
#'
#' The mean wind vector describes the net drift of the air over the series: where winds blow
#' from many directions, it is short even if winds are strong. It is similar to, but not the same
#' as, the [net_flow()] of a wind rose built from the series. Net flow describes the connectivity
#' model rather than the wind itself: it is shaped by `trans`, and with `trans = 1` it is
#' typically a few percent weaker than the mean wind, because the rose divides each wind between
#' two neighbor directions. Use `mean()` to describe the wind, and `net_flow()` to describe the
#' wind rose used for connectivity modeling.
#'
#' @param x A `wind_series`.
#' @param ... Not used.
#' @param na.rm Logical: ignore missing values when averaging?
#' @return A `wind_field` of the mean u and v components.
#' @seealso [subset_series()] to average over selected time steps, such as one season.
#' @examples
#' series <- windscape_example("wind_series")
#' prevailing <- mean(series)
#' @aliases mean,wind_series-method
#' @exportMethod mean
setMethod("mean", "wind_series", function(x, ..., na.rm = FALSE){
      check_series(x)
      n <- x@n_steps
      r <- as(x, "SpatRaster")
      u <- terra::mean(r[[seq_len(n)]], na.rm = na.rm)
      v <- terra::mean(r[[n + seq_len(n)]], na.rm = na.rm)
      out <- c(u, v)
      names(out) <- c("u", "v")
      wind_field(out)
})



#' Get the time of each step in a wind_series
#'
#' Returns the time of each step of a `wind_series`, parsed from its layer names (like
#' `"u 2000-01-01 06:00:00"`, as written by [download_wind_data()]) or, if the names hold no times,
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
      check_series(x)
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
      check_series(x)
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


# Stop if a wind_series is malformed: layer subsetting (e.g. series[[1:10]]) keeps the class
# without updating n_steps, which would silently corrupt anything that relies on it
check_series <- function(x){
      if(terra::nlyr(x) != 2 * x@n_steps)
            stop("`x` is a malformed wind_series: it has ", terra::nlyr(x), " layers but claims ",
                 x@n_steps, " time steps. This happens when layers are selected with `[[`; use ",
                 "subset_series() to select time steps.", call. = FALSE)
      invisible(x)
}
