#' An S4 object class representing a wind rose
#'
setClass("wind_rose",
         contains = "SpatRaster",
         slots = c(trans = "function",
                   n_steps = "numeric"))


#' Load a wind_rose from a raster file on disk
#'
#' @param x Path to a raster file containing wind rose data
#' @param trans Transformation to convert speed into conductance. Either a single number, or
#'    a function. See documentation for \link{wind_rose}.
#' @param n_steps An integer giving the number of time steps represented in \code{x}.
#' @return A `wind_rose` object.
#' @export
read_wind_rose <- function(x, trans = 1, n_steps = NA_integer_){
      as_wind_rose(rast(x), trans = trans, n_steps = n_steps)
}

#' Create a wind_rose object from a set of raster layers
#'
#' @param x SpatRaster with 8 layers representing flows in each semi-cardinal direction,
#'    clockwise beginning in southwest.
#' @param trans Transformation to convert speed into conductance. Either a single number, or
#'    a function. See documentation for \link{wind_rose}.
#' @param n_steps An integer giving the number of time steps represented in \code{x}.
#' @return A `wind_rose` object.
#' @export
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


#' Summarize a time series of wind fields into a wind rose
#'
#' This function converts a \code{wind_field_ts} into a \code{wind_rose} object
#' summarizing the distribution of wind speed and direction observations in each grid
#' cell. The result is a set of eight raster layers giving the average wind conductance
#' toward each of a cell's 'queen' neighbors.
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
#' @param x Data set of class `wind_series`, or a character vector of paths to files in
#'    `wind_series` layout, such as those returned by [ncar_download()]. Multiple files are
#'    processed one at a time and combined with [combine_roses()], so a long record can be
#'    summarized without loading it all into memory at once. The result is identical to
#'    building a rose from all files combined with [read_wind_series()].
#' @param trans Either a function, or a positive number indicating the power to raise windspeeds to; see details.
#' @param ... Additional arguments passed to `terra::app`, e.g. 'filename'. When `x` is a
#'    vector of files, these are passed to `terra::app` for each file, and `filename` is not
#'    allowed; use `terra::writeRaster()` on the result instead.
#' @return A \code{wind_rose} object. This is an 8-layer raster stack, where each layer is wind conductance from
#'   the focal cell to one of its neighbors (clockwise starting in the SW).
#'   If input windspeeds are in m/s and `trans = 1`, values are in (1 / hours)
#' @aliases windrose_rasters
#' @export
wind_rose <- function(x, trans = 1, ...){

      if(is.character(x)){
            if(length(x) == 1) return(wind_rose(read_wind_series(x), trans = trans, ...))
            if("filename" %in% names(list(...)))
                  stop("`filename` is not supported when `x` is a vector of files; ",
                       "use terra::writeRaster() on the result.")
            out <- NULL
            for(f in x){
                  r <- wind_rose(read_wind_series(f), trans = trans, ...)
                  out <- if(is.null(out)) r else combine_roses(out, r)
            }
            return(out)
      }
      if(!inherits(x, "wind_series")) stop("`x` must be an object of class `wind_series`, or file paths.")
      check_grid(x)

      trn <- trans
      if(is.numeric(trans)) trn <- function(x) x^trans

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

      x <- c(lat(x), rsn(x), x)
      r <- terra::app(x, fun = rose, trans = trn)
      as_wind_rose(r, trn, x@n_steps)
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
#' @seealso [wind_rose()], which uses this function to build roses from multiple files.
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
