#' An S4 object class representing a wind field time series
#'
setClass("wind_field_ts",
         contains = "SpatRaster",
         slots = c(n_steps = "numeric"))

#' Generate a wind field time series data set from a set of rasters.
#'
#' @param x Multi-layer \code{SpatRaster} with layers containing u and v wind components, or
#'    an object like a file path that can be converted to a \code{SpatRaster}. Note that
#'    raster data must be in lat-long coordinates.
#' @param order Either \code{"uuvv"}, the default indicating `x` has all u components
#'    followed by all v components, or \code{"uuvv"}, indicating the u and v components
#'    of `x` are alternating.
#' @return A `wind_field_ts` object, which is a particular form of \code{SpatRaster}.
wind_field_ts <- function(x, order){

      x <- rast(x)
      if(terra::nlyr(x) %% 2 != 0) stop("if `v` is not specified, `u` must have an even number of layers.")

      # collate layers
      if(order == "uvuv"){
            even <- function(x) x %% 2 == 0
            x <- x[[c(which(!even(1:nlyr(x))),
                      which(even(1:nlyr(x))))]]
      }

      # create wind_field
      y <- as(x, "wind_field_ts")
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
wind_field <- function(x){
      if(terra::nlyr(x) != 2) stop("`x` must have two layers.")
      xt <- ext(x)
      if(any(c(xt$xmin < -180, xt$xmax > 360, xt$ymin < -90, xt$ymax > 90))) stop("wind field rasters must be in lon-lat coordinates")
      as(x, "wind_field")
}
