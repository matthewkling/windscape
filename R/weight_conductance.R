#' Weight a wind rose's conductances
#'
#' This function adjusts the conductance values in a wind rose, multiplying them by the values of a secondary raster data set \code{w},
#' which can be used to integrate non-wind factors into a wind connectivity analysis.
#'
#' As a hypothetical example, to incorporate a decreased but nonzero likelihood of dispersal over inhospitable areas, \code{w} could be a raster layer with 0.1 indicating water or mountains and 1.0 elsewhere.
#' This would have the effect of down-weighting conductance over water or mountains by 90%.
#'
#' @param rose A `wind_rose`.
#' @param w A single-layer `SpatRaster` with values to be multiplied by \code{rose}, on the same grid as \code{rose}.
#' @return A `wind_rose`, with conductance values weighted by \code{w}.
#' @export
weight_conductance <- function(rose, w){
      if(!inherits(rose, "wind_rose")) stop("`rose` must be a wind_rose")
      as_wind_rose(as(rose, "SpatRaster") * w, trans = rose@trans, n_steps = rose@n_steps)
}
