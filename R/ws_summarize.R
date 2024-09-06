#' Summary statistics for a windshed
#'
#' A "windshed" is a raster layer representing upwind or downwind connectivity for a given point. This
#' function computes a range of statistics summarizing a windshed, such as the size, centroid location,
#' and isotropy of the wind accessibility surface. It is useful primarily for comparing the windshed
#' properties of multiple different sites.
#'
#' @param x SpatRaster representing wind accessibility, with higher values indicating greater accessibility.
#'    (For example, this could be the output of `random_walk()`, or `least_cost_surface(..., rate = TRUE)`.)
#' @param origin Coordinates of center point (matrix with 2 columns and 1 row).
#' @param radius Optional windshed radius, in km; cells farther from the origin than this will be excluded
#'    from the summary statistics.
#'
#' @return A named vector of summary statistics:
#' \itemize{
#'  \item{"centroid_x": }{Longitude of windshed "centroid". The centroid is the mean of all lon-lat coordinates, weighed by the wind accessibility values in \code{x}.}
#'  \item{"centroid_y": }{Latitude of windshed centroid.}
#'  \item{"centroid_distance": }{Distance from origin to windshed centroid, in km.}
#'  \item{"centroid_bearing": }{Bearing from origin to windshed centroid, in degrees.}
#'  \item{"windshed_distance: }{Mean distance from origin, in km. This is the average distance from the origin to every grid cell, weighted by accessibility.}
#'  \item{"windshed_bearing": }{Mean bearing from origin, in degrees. This is the circular-mean bearing from the origin to every other grid cell, weighted by accessibility.}
#'  \item{"windshed_isotropy": }{Circular standard deviation of the bearing from the origin to every other grid cell, weighted by accessibility. Low values indicate that accessible sites are concentrated in a narrow range of compass directions, while high values indicate a wide distribution of directions. The range of possible values is 0-1.}
#'  \item{"windshed_size": }{Mean accessibility. This is the average of all accessibility values, weighted to correct for grid distortions. Higher values indicate landscapes with greater overall wind accessibility.}
#'  \item{"windshed_landarea": }{A realative measure of the land area covered by grid cells with nonzero wind accessibility.}
#' }
ws_summarize <- function(x, origin, radius = NULL){

      z <- values(x %>% setNames("z"))
      p <- cbind(crds(x), z)

      # cell weights based on latitude
      y <- seq(p[1,2], p[nrow(p),2], length.out=nrow(x))
      xres2 <- res(x)[1] / 2
      yres2 <- res(x)[2] / 2
      weight <- sapply(y, function(y) geosphere::distGeo(c(0-xres2, y), c(0+xres2, y)) / 1000)
      p <- cbind(p, weight = rep(weight, each = ncol(x)))

      # exclude cells outside radius
      if(!is.null(radius)) p <- p[geosphere::distGeo(origin, crds(x)) / 1000 < radius, ]

      # cell weights based on longitude
      # downweight nearby east-west neighbors of origin pixel at higher latitudes
      # to correct bias from latitudinal cell size gradient
      ox <- origin[1]
      oy <- origin[2]
      d0 <- geosphere::distGeo(c(-xres2, 0), c(xres2, 0)) / 2 / 1000 # halfwidth of cell at equator
      dy <- geosphere::distGeo(c(-xres2, oy), c(xres2, oy)) / 2 / 1000 # halfwidth of cell at focal latitude
      ri <- abs(p[,"y"] - oy) < yres2 # cell indices in the center row
      ph <- p[ri, ]
      dp <- abs(ph[,"x"] - ox) * dy / yres2 # km to origin
      bi <- d0 - dy # distance range of cells overlapping hole
      bo <- d0 + dy
      w <- (dp - bi) / (bo - bi)
      w[w > 1] <- 1
      w[w < 0] <- 0
      p[ri,"weight"] <- p[ri,"weight"] * w

      # remove NA values
      use <- is.finite(p[,"z"])
      p <- p[use, ]
      xy <- p[,1:2]
      weights <- p[,"weight"]
      w <- p[,"z"]

      # size of windshed
      land_area <- sum(weights[w > 0])
      ws_size <- weighted.mean(w, weights)

      # distance and bearing to centroid of windshed
      ctd <- apply(xy, 2, function(z) weighted.mean(z, w*weights, na.rm=T))
      ctd_dist <- geosphere::distGeo(origin, ctd) / 1000
      ctd_brng <- geosphere::bearing(origin, ctd)

      # distance and bearing to every cell
      dist <- geosphere::distGeo(origin, xy) / 1000
      brng <- round(geosphere::bearing(origin, xy))

      # mean distance, mean bearing
      ws_dist <- weighted.mean(dist, w*weights, na.rm=T)
      iso <- circ_sd(brng, w*weights, na.rm=T)
      ws_brng <- iso["bearing"]
      ws_iso <- iso["iso"]

      # package results
      names(ctd) <- names(ws_iso) <- names(ws_brng) <- NULL
      c(centroid_x = ctd[1],
        centroid_y = ctd[2],
        centroid_distance = ctd_dist,
        centroid_bearing = ctd_brng,
        windshed_distance = ws_dist,
        windshed_bearing = ws_brng,
        windshed_isotropy = ws_iso,
        windshed_size = ws_size,
        windshed_landarea = land_area)
}
