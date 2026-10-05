#' An S4 object class representing a wind rose
#'
setClass("wind_rose",
         contains = "SpatRaster",
         slots = c(trans = "function",
                   n_steps = "numeric"))


# Create a wind_rose object from an 8-layer SpatRaster of conductances toward each neighbor,
# clockwise beginning in the southwest. `trans` (a number or function) and `n_steps` are stored
# as metadata. Internal: users create roses with wind_rose().
as_wind_rose <- function(x, trans, n_steps = NA_integer_){
      if(!inherits(x, "SpatRaster")) stop("x must be a SpatRaster")
      if(nlyr(x) != 8) stop("x must have 8 layers")
      names(x) <- c("SW", "W", "NW", "N", "NE", "E", "SE", "S")
      x <- as(x, "wind_rose")
      if(is.numeric(trans)){
            p <- trans
            trans <- function(x) x^p
      }
      x@trans <- trans
      x@n_steps <- n_steps
      x
}


#' Build or load a wind rose
#'
#' A wind rose summarizes a time series of wind fields as the average wind conductance from each
#' grid cell toward each of its eight neighbors. Given a `wind_series`, `wind_rose()` builds one;
#' given a raster or file holding a saved wind rose, it loads it.
#'
#' For each time step, the wind in each cell is divided between the two neighbors whose
#' directions bracket the wind direction, in proportion to how closely the wind points toward
#' each, and its transformed speed (see `trans`) is converted to conductance toward each
#' neighbor; conductance is then averaged over all time steps. Long series are processed in
#' chunks of time steps, so a series spanning many files (see [wind_series()]) can be summarized
#' without loading it all into memory; the result is identical to processing it at once. To
#' build a rose from selected time steps, such as a season or time of day, select them first
#' with [subset_series()].
#'
#' The `trans` parameter defines the transformation function used to convert wind speed
#' into conductance. If a numeric value is supplied, the function speed^trans is used.
#' A value of trans = 0 will ignore speed, assigning weights based on direction only;
#' trans = 1 assumes conductance is proportional to windspeed, trans = 2 assumes it's
#' proportional to aerodynamic drag, and trans = 3 assumes it's proportional to force.
#' Any intermediate value can also be used. Any function that transforms a numeric
#' vector can also be supplied; for example, to model seed dispersal for a species
#' that only releases seeds when winds exceed 10 m/s, we could specify a threshold
#' function `trans = function(x){x[x < 10] <- 0; return(x)}`.
#'
#' Grid geometry: windscape works on longitude/latitude grids with square cells. Distances
#' and bearings to each cell's neighbors are computed on the ellipsoid at that cell's latitude,
#' so conductance accounts for the narrowing of cells toward the poles. Projected grids are not
#' currently supported. If source data are projected, reproject them to longitude/latitude
#' before building a \code{wind_series}, rotating u and v to true east and north if they are
#' defined relative to the projected grid (as in some reanalysis products).
#'
#' @param x A `wind_series`, to build a wind rose from; or a saved wind rose to load, as an
#'    8-layer `SpatRaster` or the path to a raster file (e.g. one written with
#'    `terra::writeRaster()`). A saved rose's layers must hold conductance toward the
#'    southwest, west, northwest, north, northeast, east, southeast, and south neighbors, in
#'    that order.
#' @param trans Either a function, or a positive number indicating the power to raise wind
#'    speeds to; see details. When loading a saved rose, give the `trans` it was built with,
#'    which is recorded with the rose (e.g. for [combine_roses()]).
#' @param n_steps When loading a saved rose, the number of time steps it summarizes, needed to
#'    combine it with other roses using [combine_roses()]. Ignored when building a rose, which
#'    records its number of time steps automatically.
#' @param ... When building a rose, additional arguments passed to `terra::app()`, such as
#'    `cores`. A `filename` (with optional `overwrite`) writes the result to a file.
#' @return A \code{wind_rose} object. This is an 8-layer raster stack, where each layer is wind
#'   conductance from the focal cell to one of its neighbors (clockwise starting in the SW).
#'   If input windspeeds are in m/s and `trans = 1`, values are in (1 / hours)
#' @seealso [combine_roses()] to combine roses built from different time periods.
#' @examples
#' series <- windscape_example("wind_series")
#' rose <- wind_rose(series)
#'
#' # save and reload
#' f <- tempfile(fileext = ".tif")
#' terra::writeRaster(rose, f)
#' rose2 <- wind_rose(f, trans = 1, n_steps = rose@n_steps)
#' @aliases windrose_rasters
#' @export
wind_rose <- function(x, trans = 1, n_steps = NA_integer_, ...){

      # load a saved rose
      if(inherits(x, "wind_rose")) return(x)
      if(!inherits(x, "wind_series")){
            if(is.character(x)){
                  if(length(x) != 1) stop("`x` must be a single file path to load a saved wind rose; ",
                                          "to build a rose from files of wind data, use wind_series() first")
                  if(!file.exists(x)) stop("file not found: ", x)
                  x <- terra::rast(x)
            }
            if(!inherits(x, "SpatRaster"))
                  stop("`x` must be a wind_series, or a saved wind rose (an 8-layer SpatRaster or file path)")
            if(any(grepl("^[uv]( |$)", names(x))))
                  stop("`x` looks like wind data (layer names begin with u or v), not a saved wind rose; ",
                       "use wind_series() first")
            if(terra::nlyr(x) != 8) stop("a saved wind rose must have 8 layers; `x` has ", terra::nlyr(x))
            return(as_wind_rose(x, trans = trans, n_steps = n_steps))
      }

      # build a rose from a wind_series
      check_series(x)
      check_grid(x)
      trn <- trans
      if(is.numeric(trans)) trn <- function(x) x^trans

      dots <- list(...)
      filename <- dots$filename
      overwrite <- isTRUE(dots$overwrite)
      dots$filename <- dots$overwrite <- NULL

      # chunks of time steps, each holding at most about `budget` values
      n <- x@n_steps
      budget <- getOption("windscape.chunk_values", 5e7)
      per_chunk <- max(1, floor(budget / (2 * terra::ncell(x))))
      chunks <- split(seq_len(n), ceiling(seq_len(n) / per_chunk))

      out <- NULL
      for(steps in chunks){
            xs <- if(length(chunks) == 1) x else subset_series(x, steps = steps)
            r <- do.call(rose_chunk, c(list(xs, trn), dots))
            out <- if(is.null(out)) r else combine_roses(out, r)
      }
      if(!is.null(filename)){
            out <- as_wind_rose(terra::writeRaster(out, filename, overwrite = overwrite),
                                trans = trn, n_steps = n)
      }
      out
}

# Build a rose from one chunk of a wind_series
rose_chunk <- function(x, trn, ...){
      rsn <- function(x){
            r <- x[[1]]
            r[] <- mean(res(r[[1]]))
            r
      }
      lat <- function(x){
            l <- x[[1]]
            l[] <- terra::crds(l)[,2]
            l
      }
      n <- x@n_steps
      x <- c(lat(x), rsn(x), x)
      r <- terra::app(x, fun = rose, trans = trn, ...)
      as_wind_rose(r, trn, n)
}



#' Combine wind roses built from different time periods
#'
#' Combines wind roses for the same grid, built from different sets of time steps (e.g. separate
#' months or years), into a single rose representing all of those time steps. Because a wind rose
#' is an average over time steps, the result is the mean of the input roses weighted by their
#' numbers of time steps, and equals the rose that would be built from all the time steps at once.
#'
#' @param ... Two or more `wind_rose` objects, or a single list of them. They must share the same
#'    grid and the same `trans` function, and each must have a known number of time steps
#'    (`n_steps`), as roses built with [wind_rose()] do.
#' @return A `wind_rose` whose `n_steps` is the total across the inputs.
#' @seealso [wind_rose()], which uses this function to build roses from long series in chunks.
#' @export
combine_roses <- function(...){
      x <- list(...)
      if(length(x) == 1 && is.list(x[[1]]) && !inherits(x[[1]], "SpatRaster")) x <- x[[1]]
      if(length(x) < 2) stop("at least two wind roses are needed")
      if(!all(vapply(x, inherits, logical(1), "wind_rose"))) stop("all inputs must be `wind_rose` objects")
      n <- vapply(x, function(r) as.numeric(r@n_steps), numeric(1))
      if(anyNA(n)) stop("all roses must have a known number of time steps (`n_steps`)")
      probe <- c(0, 0.5, 1, 2, 5, 10, 20, 50)
      t1 <- x[[1]]@trans(probe)
      for(r in x[-1]){
            if(!terra::compareGeom(x[[1]], r, stopOnError = FALSE)) stop("roses are on different grids")
            if(!isTRUE(all.equal(t1, r@trans(probe)))) stop("roses use different `trans` functions")
      }
      out <- x[[1]] * n[1]
      for(i in seq_along(x)[-1]) out <- out + x[[i]] * n[i]
      as_wind_rose(out / sum(n), trans = x[[1]]@trans, n_steps = sum(n))
}
