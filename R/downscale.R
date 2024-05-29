
#' Downscale a wind rose raster data to higher spatial resolution
#'
#' This function disaggregates a wind_rose raster, interpolating it to higher spatial resolution and then
#' adjusting values to maintain correct connectivity units. Downscaling is useful when connectivity is to
#' be modeled between sites that are separated by few or no grid cells.
#'
#' @param x wind_rose
#' @param fact a single integer giving the factor by which to disaggregate the data; this is passed to \link[terra]{disagg}.
#' @param method disaggregation method; either "near" or "bilinear'; see \link[terra]{disagg} for details.
#' @return A higher resolution version of x
#' @details
#' Larger values of \code{fact} create higher-resolution wind surfaces, which improves connectivity estimates
#' between nearby sites. The trade-off is computational cost.
#'
#' @export
downscale <- function(x, fact, method = "bilinear"){
      if(!inherits(x, "wind_rose")) stop("x must be a wind_rose object")
      if(length(fact) != 1) stop("fact must be a single integer")
      disagg(x, fact, method = method) * fact
}
