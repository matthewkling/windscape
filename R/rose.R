#' Calculate 8-neighbor edge loadings from a time series of u and v windspeeds
#'
#' @param x A vector of wind data containing: latitude, resolution, u
#'   windspeeds, v windspeeds
#' @param p A positive number indicating the power to raise windspeeds to (see details)
#' @return A vector of 8 conductance values to neighboring cells, clockwise
#'   starting with the southwest neighbor. If input windspeeds are in m/s,
#'   values are in (1 / hours ^ p)
#' @details A value of p = 0 will ignore speed, assigning weights based on
#'   direction only. P = 1 assumes conductance is proportional to windspeed, p =
#'   2 assumes it's proportional to aerodynamic drag, and p = 3 assumes it's
#'   proportional to force. Any intermediate value can also be used.
#' @export
rose <- function(x, trans = identity){

      # unpack & restructure: row=timestep, col=u&v components
      lat <- x[1]
      res <- x[2]
      x <- x[3:length(x)]
      uv <- matrix(x, ncol=2, byrow=F)

      # wind speed and direction
      speed <- sqrt(uv[,1]^2 + uv[,2]^2)
      weight <- trans(speed)
      dir <- spin90(windscape::direction(uv[,2], -1*uv[,1]))
      dir[dir<0] <- dir[dir<0] + 360
      dir[dir==0] <- 360

      # bearings to queen neighbors
      nc <- cbind(x = c(0, res, res, res, 0, -res, -res, -res),
                  y = c(res, res, 0, -res, -res, -res, 0, res) + lat)
      nc[,2] <- pmin(nc[,2], 90)
      nc[,2] <- pmax(nc[,2], -90)
      nb <- geosphere::bearingRhumb(c(0, lat), nc)
      nb <- c(nb, 360)

      # allocate weights to neighbors, based on direction
      l <- edge_loadings(dir, weight, nb)

      # adjust conductance weights to account for inter-cell distance and n observations
      nd <- geosphere::distGeo(c(0, lat), nc) # distance, in m
      l <- l * 3600 / nd / nrow(uv) # units are now in 3600 / seconds^p, which is 1 / hours if p = 1

      # reorder, clockwise from SW
      l[c(6:8, 1:5)]
}


#' @useDynLib windscape
#' @importFrom Rcpp sourceCpp
NULL
