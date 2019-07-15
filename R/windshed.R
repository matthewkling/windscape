

# function to supply to the transitionFunction arg of transition_stack
windflow <- function(x, direction="downwind"){
      # x is vector of collated windrose values for FROM and TO cells
      # direction is downwind or upwind
      # node row and column indices
      g <- x[17:20]

      # edge wind loadings (clockwise from southwest)
      p <- c(x[1]+x[2],
             x[3]+x[4],
             x[5]+x[6],
             x[7]+x[8],
             x[9]+x[10],
             x[11]+x[12],
             x[13]+x[14],
             x[15]+x[16]) / 2

      if(direction=="upwind") p <- p[c(5:8, 1:4)]

      if(g[1]<g[2] & g[4]<g[3]) return(p[1]) #SW
      if(g[1]==g[2] & g[4]<g[3]) return(p[2]) #W
      if(g[2]<g[1] & g[4]<g[3]) return(p[3]) #NW
      if(g[3]==g[4] & g[2]<g[1]) return(p[4]) #N
      if(g[2]<g[1] & g[3]<g[4]) return(p[5]) #NE
      if(g[1]==g[2] & g[3]<g[4]) return(p[6]) #E
      if(g[1]<g[2] & g[3]<g[4]) return(p[7]) #SE
      if(g[3]==g[4] & g[1]<g[2]) return(p[8]) #S
}


transition_stack <- function(x, transitionFunction, directions, symm, ...){

      brk <- x
      x <- brk[[1]]

      tr <- new("TransitionLayer",
                nrows=as.integer(nrow(x)),
                ncols=as.integer(ncol(x)),
                extent=extent(x),
                crs=projection(x, asText=FALSE),
                transitionMatrix = Matrix(0,ncell(x),ncell(x)),
                transitionCells = 1:ncell(x))
      transitionMatr <- transitionMatrix(tr)
      Cells <- which(!is.na(getValues(x)))
      adj <- adjacent(x, cells=Cells, pairs=TRUE, target=Cells, directions=directions)

      ##### start modificiations #####
      # format raster data layers to feed to transitionFunction
      # col order is x[[1]][from, to] ... x[[n]][from,to]
      dataVals <- lapply(1:nlayers(brk),
                         function(i) cbind(values(brk[[i]])[adj[,1]],
                                           values(brk[[i]])[adj[,2]]))
      dataVals <- do.call("cbind", dataVals)
      ##### end modifications #####

      transition.values <- apply(dataVals, 1, transitionFunction, ...)
      transitionMatr[adj] <- as.vector(transition.values)
      transitionMatrix(tr) <- transitionMatr
      matrixValues(tr) <- "resistance"
      return(tr)
}



ws_summarize <- function(x, # raster layer of wind flow (where positive values are more accessible)
                         origin # coordinates of center point (2-column matrix)
){
      p <- cbind(coordinates(x), values(x))

      # cell weights based on latitude
      y <- seq(p[1,2], p[nrow(p),2], length.out=nrow(x))
      inc <- res(x)[1] / 2
      weight <- sapply(y, function(y) distGeo(c(0-inc, y), c(0+inc, y)))
      p <- cbind(p, weight = rep(weight, each=ncol(x)))

      use <- is.finite(p[,3])
      p <- p[use, ]
      xy <- p[,1:2]
      weights <- p[,4]
      w <- p[,3]

      # size of windshed
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


      names(ctd) <- names(ws_iso) <- names(ws_brng) <- NULL
      c(centroid_x = ctd[1],
        centroid_y = ctd[2],
        centroid_distance = ctd_dist,
        centroid_bearing = ctd_brng,
        windshed_distance = ws_dist,
        windshed_bearing = ws_brng,
        windshed_isotropy = ws_iso,
        windshed_size=ws_size)
}


# weighted circular standard deviation
circ_sd <- function(x, # bearings -- in degrees, not radians
                     w=NULL,
                     ...){
      # a custom weighted version of the "mean resultant vector" from circular stats
      # see: https://doi.org/10.3389/fpsyg.2018.02040 and help(CircStats::circ.disp)
      # my implementation: convert the bearings into XY, then take the mean of the XYs, weight by wind speed
      if(is.null(w)) w <- rep(1, length(x))
      x <- x / 180 * pi
      xy <- cbind(x = sin(x), y = cos(x))
      xy <- apply(xy, 2, weighted.mean, w=w, ...)
      rbar <- sqrt(sum(xy^2)) # mean resultant vector
      iso <- sqrt(1-rbar) # circular standard deviation
      angle <- bearing(c(0,0), xy)
      return(c(bearing=angle, iso=iso))
}


# circular distance measure used internally within isotropy function below
cdist <- function(x, y){
      d <- abs(x - y)
      d[d > 180] <- 360 - d[d > 180]
      sum(d)
}

# circular "earth movers distance" between a set of bearings and a uniform distribution,
# normalized such that 1 represents a single angle and 0 represents a uniform distirbution
anisotropy <- function(x, w=NULL, ...){
      if(is.null(w)) w <- rep(1, length(x))
      if(!is.finite(x[1] + w[1])) return(NA)
      e <- cbind(w/sum(w), x)
      u <- cbind(1/360, 1:360) # uniform distribution
      d <- emdist::emd(e, u, dist=cdist)
      d / 90 # normalize by 90, the max possible value
}


