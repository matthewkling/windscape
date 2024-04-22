setClass("wind_field",
         contains = "SpatRaster",
         slots = c(n_steps = "numeric"))

#' Generate a wind field data set from a set of rasters
#'
#' @param x Multi-layer `SpatRaster` with layers containing u and v wind components;
#'    note that these must be in lat-long coordinates.
#' @param order Either "uuvv", the default indicating `x` has all u components
#'    followed by all v components, or "uvuv", indicating the u and v components
#'    of `x` are alternating.
#' @return A `wind_field` object, which is a particular form of `SpatRaster`.
wind_field <- function(x, order){

      if(terra::nlyr(x) %% 2 != 0) stop("if `v` is not specified, `u` must have an even number of layers.")

      # collate layers
      if(order == "uvuv"){
            even <- function(x) x %% 2 == 0
            x <- x[[c(which(!even(1:nlyr(x))),
                      which(even(1:nlyr(x))))]]
      }

      # create wind_field
      y <- as(x, "wind_field")
      y@n_steps <- nlyr(x)/2
      y
}
