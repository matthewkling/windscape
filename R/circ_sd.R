#' Weighted circular standard deviation
#'
#' @param x Vector of bearings (in degrees)
#' @param w Optional vector of weights
#' @param ... Additional arguments to weighted.mean, e.g. `na.rm=TRUE`
#'
#' @return A named vector including the bearing of the mean resultant vector
#'   ("bearing") and the circular standard deviation ("iso")
circ_sd <- function(x, # bearings -- in degrees, not radians
                    w=NULL,
                    ...){
      # a custom weighted version of the "mean resultant vector" from circular stats
      # see: https://doi.org/10.3389/fpsyg.2018.02040 and help(CircStats::circ.disp)
      # my implementation: convert the bearings into XY, then take the mean of the XYs, weight by wind speed
      if(is.null(w)) w <- rep(1, length(x))
      x <- x / 180 * pi
      xy <- cbind(x = sin(x), y = cos(x))
      xy <- apply(xy, 2, weighted.mean, w=w, ...)  # mean resultant vector
      rbar <- sqrt(sum(xy^2)) # mean resultant vector length
      iso <- sqrt(1-rbar) # circular standard deviation
      angle <- bearing(c(0,0), xy)
      return(c(bearing=angle, iso=iso))
}
