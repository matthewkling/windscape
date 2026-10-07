# Wind rose computation. A rose holds, for each cell, the average over time steps of wind
# conductance toward each of its 8 neighbors: at each step the wind is split between the two
# neighbors whose bearings bracket its direction (see rose_accumulate() in src/rose.cpp),
# weighted by transformed speed, and divided by the distance to the neighbor.


# Neighbor geometry on a longitude/latitude grid, which depends only on latitude. For cell
# centers at latitudes `lat` on square cells of size `res` degrees, returns a list of matrices
# with one row per latitude: `nb`, rhumb-line bearings (degrees) to the N, NE, E, SE, S, SW, W,
# and NW neighbors, then 360; and `nd`, geodesic distances (m) to those neighbors.
neighbor_geometry <- function(lat, res){
      dx <- c(0, res, res, res, 0, -res, -res, -res)
      dy <- c(res, res, 0, -res, -res, -res, 0, res)
      geo <- lapply(lat, function(y){
            nc <- cbind(dx, pmax(pmin(dy + y, 90), -90))
            list(nb = c(geosphere::bearingRhumb(c(0, y), nc), 360),
                 nd = geosphere::distGeo(c(0, y), nc))
      })
      list(nb = do.call(rbind, lapply(geo, `[[`, "nb")),
           nd = do.call(rbind, lapply(geo, `[[`, "nd")))
}

# Add a block of time steps to running neighbor-loading sums `acc` (cells x 8, modified in
# place). `m` is a cells x (2 * k) matrix of k u columns then k v columns; `trans` is a number
# (a power of speed, computed in C++) or an elementwise function of speed.
rose_add <- function(acc, m, k, trans, nb, row){
      if(is.numeric(trans)){
            rose_accumulate(m, k, matrix(0, 0, 0), trans, nb, row, acc)
      }else{
            w <- trans(wind_speeds(m, k))
            if(!is.numeric(w) || length(w) != nrow(m) * k)
                  stop("`trans` must be an elementwise function, returning a numeric vector the ",
                       "same length as its input", call. = FALSE)
            dim(w) <- c(nrow(m), k)
            storage.mode(w) <- "double"
            rose_accumulate(m, k, w, 1, nb, row, acc)
      }
      invisible(acc)
}

# Convert running sums of loadings over `n` time steps into mean conductance: per hour if
# speeds are in m/s and trans = 1. Returns cells x 8 in rose layer order (SW, W, NW, N, NE, E,
# SE, S). `nd` and `row` are as in neighbor_geometry() and rose_add().
rose_finish <- function(acc, n, nd, row){
      out <- acc * 3600 / nd[row + 1L, , drop = FALSE] / n
      out[, c(6:8, 1:5), drop = FALSE]
}


#' Wind rose for a single cell
#'
#' Internal, single-cell version of the computation [wind_rose()] runs over whole grids; used
#' in tests.
#'
#' @param x A vector of wind data containing: latitude, resolution, u
#'   windspeeds, v windspeeds
#' @param trans A number (power of speed) or elementwise function transforming wind speed into
#'   conductance; see [wind_rose()].
#' @return A vector of 8 conductance values to neighboring cells, clockwise
#'   starting with the southwest neighbor. If input windspeeds are in m/s and trans = 1,
#'   values are in 1 / hours
#' @noRd
rose <- function(x, trans = identity){
      geo <- neighbor_geometry(x[1], x[2])
      k <- (length(x) - 2) / 2
      acc <- matrix(0, 1, 8)
      rose_add(acc, matrix(x[-(1:2)], nrow = 1), k, trans, geo$nb, 0L)
      as.vector(rose_finish(acc, k, geo$nd, 0L))
}
