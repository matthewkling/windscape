
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


# direction from u&v components -- in degrees, clockwise from 12:00
direction <- function(x, y) atan2(x, y) * 180 / pi
spin90 <- function(x){
      x <- x - 90
      x[x<(-180)] <- x[x<(-180)] + 360
      x
}

add_lat <- function(x){
      lat <- x[[1]]
      lat[] <- coordinates(lat)[,2]
      stack(lat, x)
}


neighbor_loadings <- function(b, # wind bearings
                              nb # bearings to neighbors
){
      ni <- max(which(b > nb))
      prop <- (b - nb[ni]) / (nb[ni+1] - nb[ni])
      l <- rep(0, 9)
      l[c(ni, ni+1)] <- c(1-prop, prop)
      l <- c(l[1] + l[9], l[2:8])
      as.vector(l)
}


# velocity-weighted frequency of wind in each quadrant
# same as above, but accounting for non-square grid
windrose_geo <- function(x,
                         res, # raster resolution, in degrees
                         p=2, # 0=time, 1=velocity, 2=drag, 3=force
                         summary_fun = sum # function to summarize across time steps
){

      wfx <- function(x) sqrt(sum(x^2)) ^ p

      # unpack & restructure: row=timestep, col=u&v components
      lat <- x[1]
      x <- x[2:length(x)]
      m <- matrix(x, ncol=2, byrow=F)

      # wind force and direction
      weight <- apply(m, 1, wfx)
      dir <- apply(m, 1, function(x) spin90(direction(x[2], -1*x[1])))
      dir[dir<0] <- dir[dir<0] + 360

      # bearings to queen neighbors
      nc <- cbind(x = c(0, res, res, res, 0, -res, -res, -res),
                  y = c(res, res, 0, -res, -res, -res, 0, res) + lat)
      nc[,2] <- pmin(nc[,2], 90)
      nc[,2] <- pmax(nc[,2], -90)
      nb <- bearingRhumb(c(0, lat), nc)
      nb <- c(nb, 360)

      # allocate weights to neighbors, based on direction
      l <- t(apply(matrix(dir, ncol=1), 1, neighbor_loadings, nb=nb))
      l <- sweep(l, 1, weight, `*`)
      l <- apply(l, 2, summary_fun)

      # adjust conductance weights to account for inter-cell distance
      nd <- distGeo(c(0, lat), nc)
      l <- l / nd

      # reorder, clockwise from SW
      l[c(6:8, 1:6)]
}



# generate windoses for a raster dataset
# must be in lat-long proection
windrose_rasters <- function(w, # raster stack of u and v wind components
                             outfile,
                             order = "uvuv", # either "uvuv" or "uuvv", indicating whether u and v components are alternating
                             ncores = 1,
                             ... # see windrose function
){

      # collate data
      if(order == "uvuv"){
            even <- function(x) x %% 2 == 0
            w <- w[[c(which(!even(1:nlayers(w))),
                      which(even(1:nlayers(w))))]]
      }

      resn <- mean(res(w[[1]]))
      rosefun <- function(x) windrose_geo(x, res=resn, ...)

      if(ncores == 1){
            wr <- raster::calc(w, fun=rosefun, forceapply=TRUE,
                               filename=outfile)
            names(wr) <- c("SW", "W", "NW", "N", "NE", "E", "SE", "S")
            return(wr)
      } else {

            # split dataset into batches
            nlayer <- nlayers(w)/2
            core <- rep(1:ncores, each=floor(nlayer/ncores))
            core <- c(core, rep(ncores, nlayer %% ncores))
            w <- lapply(1:ncores, function(x) list(i=x,
                                                   data=w[[c(which(core == x),
                                                             which(core == x) + nlayer)]]))

            # process batches in parallel
            require(doParallel)
            cl <- makeCluster(ncores)
            registerDoParallel(cl)
            wr <- foreach(x = w,
                          .packages=c("raster", "windscape")) %dopar% {
                                x$data <- add_lat(x$data)
                                raster::calc(x$data, fun=rosefun, forceapply=TRUE,
                                             filename=paste0(substr(outfile, 1, nchar(outfile)-4),
                                                             "_", x$i,
                                                             substr(outfile, nchar(outfile)-3, nchar(outfile))))
                          }
            stopCluster(cl)

      }
}


add_coords <- function(windrose){
      rows <- cols <- windrose[[1]]
      rows[] <- rep(1:nrow(rows), each=ncol(rows))
      cols[] <- rep(1:ncol(rows), nrow(rows))
      windrose <- stack(windrose, rows, cols)
      names(windrose) <- c("SW", "W", "NW", "N", "NE", "E", "SE", "S", "row", "col")
      return(windrose)
}

windrose_names <- function() c("SW", "W", "NW", "N", "NE", "E", "SE", "S")
windrose_bearings <- function() c(225, 270, 315, 0, 45, 90, 135, 180)





wind_stats <- function(wr){
      require(ineq)
      total <- sum(wr)
      directionality <- calc(wr, Gini) # Gini: 1 = directionality, 0 = equality
      mean_u
      mean_v
      mean_direction
      mean_speed
}
