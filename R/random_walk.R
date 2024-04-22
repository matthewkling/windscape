
#' Simulate wind diffusion
#'
#' This function estimates diffusion across a windscape by Markov random walk. On each iteration,
#' the "particle mass" in each grid cell is dispersed within the local 9-cell neighborhood in
#' proportion with wind conductance. Depending on usage, this "mass" could represent probability,
#' number of individuals, etc. As the simulation proceeds, the mass diffuses across the landscape.
#'
#' The input wind rose raster is converted to a simplex of nine probabilities for each grid
#' cell, giving the rates at which particles are retained in a cell or moved to each of its
#' eight neighbors. Probabilities of moving to a neighboring cell are proportional to conductance
#' in the wind rose data set, with the probabilities of remaining in a cell scaled so that the
#' cell with the highest conductance has a zero retention probability; this allows conductance
#' to be normalized locally to a simplex while remaining proportional across cells and
#' maximizing the dispersal occurring at each iteration.
#'
#' By default, the simulation tracks the fleeting diffusion of the initial "pulse" of particle
#' mass, which will drift and spread like a cloud, eventually leaving the modeling domain
#' (unless it's a closed spatial domain). The alternative setting `ratchet = TRUE` instead
#' runs a propagating simulation in which local particle mass never declines, continuing to
#' transmit mass at the cumulative maximum rate; this may be more akin to a biological
#' process in which particles can reproduce locally after establishment.
#'
#' @param x a `wind_rose`.
#' @param init Initial conditions from which to begin diffusion, either a two-column
#'    matrix of coordinates, or a SpatRaster layer with non-negative mass values and the
#'    same spatial properties as `x`. If coordinates, the simulation starts with a
#'    mass of 1 at each coordinate location. If a SpatRaster, diffusion is done
#'    directly on the raster; this could be the result of a prior random_walk, or any
#'    other data representing quantities to be spatially dispersed.
#' @param iter Number of simulation iterations (positive integer).
#' @param record Integer vector specifying which iterations to record. Default is to
#'    record only the final iteration, i.e. `iter`. One raster layer is returned for
#'    each value in `record`.
#' @param ratchet Whether or not to "ratchet up" values during diffusion. Default is
#'    FALSE. See details.
random_walk <- function(x, init, iter = 100, record = iter, ratchet = FALSE){

      p <- transition_prob(x)

      # starting distribution
      if(inherits(init, "matrix")){
            init <- raster::cellFromXY(x, init)
            n <- x[[1]]
            values(n) <- 0
            n[unique(init)] <- 1
      }else{
            n <- init
      }

      w <- walk(n, p, iter, record, ratchet)
      n <- rast(n, nlyrs = length(record), vals = w)
      names(n) <- paste0("iter", record)
      as(n, "SpatRaster")
}

transition_prob <- function(x){
      p <- sum(x)
      p <- minmax(p)[2,] - p
      p <- c(p, x)
      p / sum(p)
}

disperse <- function(n, p){
      a <- p
      for(k in 1:9) a[,,k] <- n * a[,,k]
      da <- dim(a)
      a[,,2] <- rbind(rep(0, da[2]), cbind(a[1:(da[1]-1), 2:(da[2]), 2], rep(0, da[1]-1))) # SW
      a[,,3] <- cbind(a[1:(da[1]), 2:(da[2]), 3], rep(0, da[1])) # W
      a[,,4] <- rbind(cbind(a[2:(da[1]), 2:(da[2]), 4], rep(0, da[1]-1)), rep(0, da[2])) # NW
      a[,,5] <- rbind(a[2:(da[1]), 1:(da[2]), 5], rep(0, da[2])) # N
      a[,,6] <- rbind(cbind(rep(0, da[1]-1), a[2:(da[1]), 1:(da[2]-1), 6]), rep(0, da[2])) # NE
      a[,,7] <- cbind(rep(0, da[1]), a[1:(da[1]), 1:(da[2]-1), 7]) # E
      a[,,8] <- rbind(rep(0, da[2]), cbind(rep(0, da[1]-1), a[1:(da[1]-1), 1:(da[2]-1), 8])) # SE
      a[,,9] <- rbind(rep(0, da[2]), a[1:(da[1]-1), 1:(da[2]), 9]) # S
      apply(a, c(1, 2), sum)
}

walk <- function(n, p, i, rec = i, ratch = F){
      r <- as.array(rast(n, nlyrs = length(rec), vals = 0))
      n <- matrix(n, nrow(n), byrow = T)
      p <- as.array(p)
      pb <- txtProgressBar(min = 0, max = i, initial = 0, style = 3)
      for(j in 1:i){
            if(ratch){
                  n <- pmax(n, disperse(n, p))
            }else{
                  n <- disperse(n, p)
            }
            if(j %in% rec) r[,,match(j, rec)] <- n
            setTxtProgressBar(pb, j+1)
      }
      close(pb)
      r
}
