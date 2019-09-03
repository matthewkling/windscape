#' Summarize a windshed raster
#'
#' @param x raster layer of wind flow (where positive values are more accessible)
#' @param origin coordinates of center point (2-column matrix)
#'
#' @return A named vector of summary stats
ws_summarize <- function(x, origin){

      p <- cbind(coordinates(x), z = values(x))

      # cell weights based on latitude
      y <- seq(p[1,2], p[nrow(p),2], length.out=nrow(x))
      xres2 <- res(x)[1] / 2
      yres2 <- res(x)[2] / 2
      weight <- sapply(y, function(y) distGeo(c(0-xres2, y), c(0+xres2, y)) / 1000)
      p <- cbind(p, weight = rep(weight, each=ncol(x)))

      # cell weights based on longitude
      # downweight nearby east-west neigbors of origin pixel at higher latitudes
      # to correct bias from latitudinal cell size gradient
      ox <- origin[1]
      oy <- origin[2]
      d0 <- distGeo(c(-xres2, 0), c(xres2, 0)) / 2 / 1000 # halfwidth of cell at equator
      dy <- distGeo(c(-xres2, oy), c(xres2, oy)) / 2 / 1000 # halfwidth of cell at focal latitude
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
      ctd_dist <- distGeo(origin, ctd) / 1000
      ctd_brng <- bearing(origin, ctd)

      # distance and bearing to every cell
      dist <- distGeo(origin, xy) / 1000
      brng <- round(bearing(origin, xy))

      # mean distance, mean bearing
      ws_dist <- weighted.mean(dist, w*weights, na.rm=T)
      iso <- circ_sd(brng, w*weights, na.rm=T)
      ws_brng <- iso["bearing"]
      #ws_iso <- 1 - anisotropy(brng, w) # too slow
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
        windshed_size=ws_size,
        windshed_landarea=land_area)
}
