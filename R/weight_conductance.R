
#' Weight a windrose conductance raster
#'
#' This function adjusts the conductance values in a windrose raster, multiplying them by the values of a secondary raster data set \code{w},
#' which can be used to integrate non-wind factors into a wind connectivity analysis.
#'
#' As a hypothetical example, to incorporate a decreased but nonzero likelihood of dispersal over inhospitable areas, \code{w} could be a raster layer with 0.1 indicating water or mountains and 1.0 elsewhere.
#' This would have the effect of down-weighting conductance over water or mountains by 90%.
#'
#' @param rose A raster stack created using \code{windrose_rasters()}.
#' @param w A raster layer with values to be multiplied by \code{rose}. This layer must have the same spatial properties as the windrose raster.
#' @return A version of \code{rose}, with conductance values weighted by \code{w}.
#' @export
weight_conductance <- function(rose, w){
      rose * w
}
