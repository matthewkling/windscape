setClass("wind_rose",
         contains = "SpatRaster",
         slots = c(trans = "function",
                   n_steps = "numeric"))



#' Create a wind_rose object from a set of raster layers
#'
#' @param x SpatRaster with 8 layers representing flows in each semi-cardinal direction,
#'    clockwise beginning in southwest.
#' @param p Windspeed power, see documentation for `wind_rose`.
#' @return A `wind_rose` `SpatRaster` object.
as_wind_rose <- function(x, trans, n_steps = NULL){
      if(!inherits(x, "SpatRaster")) stop("x must be a SpatRaster")
      if(nlyr(x) != 8) stop("x must have 8 layers")
      names(x) <- c("SW", "W", "NW", "N", "NE", "E", "SE", "S")
      x <- as(x, "wind_rose")
      x@trans <- trans
      x@n_steps <- n_steps
      x
}


#' Summarize a set of wind_field rasters into a wind_rose object
#'
#' This function creates a `wind_rose` object summarizing the distribution of
#' wind speed and direction observations in each grid cell. The result is a set
#' of eight raster layers giving the average wind conductance toward each of a
#' cell's 'queen' neighbors.
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
#' @param x Data set of class `wind_field`.
#' @param trans Either a function, or a positive number indicating the power to raise windspeeds to; see details.
#' @param ... Additional arguments passed to `terra::app`, e.g. 'filename'.
#' @return An 8-layer raster stack, where each layer is wind conductance from
#'   the focal cell to one of its neighbors (clockwise starting in the SW).
#'   If input windspeeds are in m/s and `trans = 1`, values are in (1 / hours)
#' @aliases windrose_rasters
wind_rose <- function(x, trans = 1, ...){

      if(!inherits(x, "wind_field")) stop("`x` must be an object of class `wind_field`.")

      xt <- ext(x)
      if(any(c(xt$xmin < -180, xt$xmax > 360, xt$ymin < -90, xt$ymax > 90))) stop("wind field rasters must be in lon-lat coordinates")

      trn <- trans
      if(inherits(trans, "numeric")) trn <- function(x) x^trans

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
