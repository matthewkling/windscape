
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

# convert an angle [0-45] to a loading [0-1] in the 45-degree direction
loading <- function(theta){ # theta in [0,45]
      sin((4*(theta-22.5))*pi/180)/2+.5
}

#velocity-weighted frequency of wind in each quadrant
windrose <- function(x, weighting="velocity"){
      #x0 <- s[277,349]

      wfx <- switch(weighting,
                    time = function(x) 1,
                    velocity = function(x) sqrt(sum(x^2)),
                    force = function(x) sqrt(sum(x^2))^3)

      # restructure: row=timestep, col=u&v components
      m <- matrix(x, ncol=2, byrow=F)

      # calcualte force and direction
      frc <- apply(m, 1, wfx) # fun should be force or velocity
      dir <- apply(m, 1, function(x) spin90(direction(x[2], -1*x[1])))

      # octant bounds for each vector
      d1 <- plyr::round_any(dir, 45, floor)
      d1[d1 == -180] <- 180
      d2 <- plyr::round_any(dir, 45, ceiling)
      d2[d2 == -180] <- 180

      # loadings on each octant bound, based on angle
      theta <- dir %% 45
      l1 <- loading(theta)
      l2 <- 1-l1

      # add some dummy data to ensure all quadrants are represented
      frc <- c(frc, rep(0,8))
      l1 <- c(l1, rep(0,8))
      l2 <- c(l2, rep(0,8))
      d1 <- c(d1, c(-135, -90, -45, 0, 45, 90, 135, 180))
      d2 <- c(d2, c(-135, -90, -45, 0, 45, 90, 135, 180))

      # repackage
      d <- c(d1, d2)
      l <- c(l1, l2)
      frc <- c(frc, frc)

      # summarize
      d <- split(frc*l, d)
      d <- unlist(lapply(d, sum))
      return(as.numeric(d))
}




add_coords <- function(windrose){
      rows <- cols <- windrose[[1]]
      rows[] <- rep(1:nrow(rows), each=ncol(rows))
      cols[] <- rep(1:ncol(rows), nrow(rows))
      windrose <- stack(windrose, rows, cols)
      names(windrose) <- c("SW", "W", "NW", "N", "NE", "E", "SE", "S", "row", "col")
      return(windrose)
}




# generate windoses for a raster dataset
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

      rosefun <- function(x) windrose(x, ...)

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
                          .packages=c("raster", "windshed")) %dopar% {
                                raster::calc(x$data, fun=rosefun, forceapply=TRUE,
                                             filename=paste0(substr(outfile, 1, nchar(outfile)-4),
                                                             "_", x$i,
                                                             substr(outfile, nchar(outfile)-3, nchar(outfile))))
                          }
            stopCluster(cl)
            #wr <- Reduce("+", wr)
      }
}




wind_stats <- function(wr){
      require(ineq)
      total <- sum(wr)
      directionality <- calc(wr, Gini) # Gini: 1 = directionality, 0 = equality
      mean_u
      mean_v
      mean_direction
      mean_speed
}
