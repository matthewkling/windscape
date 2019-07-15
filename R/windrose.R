
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


# direction from u&v components -- in degrees, clockwise from north
direction <- function(x, y) atan2(x, y) * 180 / pi
spin90 <- function(x){
      x <- x - 90
      x[x<(-180)] <- x[x<(-180)] + 360
      x
}



add_coords <- function(windrose){
      rows <- cols <- windrose[[1]]
      rows[] <- rep(1:nrow(rows), each=ncol(rows))
      cols[] <- rep(1:ncol(rows), nrow(rows))
      windrose <- stack(windrose, rows, cols)
      names(windrose) <- c("SW", "W", "NW", "N", "NE", "E", "SE", "S", "row", "col")
      return(windrose)
}

add_res <- function(x){
      r <- x[[1]]
      r[] <- mean(res(r[[1]]))
      stack(r, x)
}


add_lat <- function(x){
      lat <- x[[1]]
      lat[] <- coordinates(lat)[,2]
      stack(lat, x)
}


# velocity-weighted frequency of wind in each quadrant
windrose <- function(x,
                          p=1 # 0=time, 1=velocity, 2=drag, 3=force
){

      # unpack & restructure: row=timestep, col=u&v components
      lat <- x[1]
      res <- x[2]
      x <- x[3:length(x)]
      uv <- matrix(x, ncol=2, byrow=F)

      # wind strength and direction
      weight <- sqrt(uv[,1]^2 + uv[,2]^2)^p
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
      # l <- t(sapply(dir, neighbor_loadings, nb = nb))
      # l <- sweep(l, 1, weight, `*`)
      # l <- apply(l, 2, summary_fun)
      l <- windscapeRCPP::edge_loadings(dir, weight, nb)

      # adjust conductance weights to account for inter-cell distance
      nd <- geosphere::distGeo(c(0, lat), nc)
      l <- l/nd
      # units are now in 1/seconds^p

      # reorder, clockwise from SW
      l[c(6:8, 1:5)]
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


      rosefun <- function(x) windrose(x, ...)

      if(ncores == 1){
            wr <- add_res(w)
            wr <- add_lat(wr)
            wr <- raster::calc(wr, fun=rosefun, forceapply=TRUE, filename=outfile)
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
                                x$data <- add_res(x$data)
                                x$data <- add_lat(x$data)
                                raster::calc(x$data, fun=rosefun, forceapply=TRUE,
                                             filename=paste0(substr(outfile, 1, nchar(outfile)-4),
                                                             "_", x$i,
                                                             substr(outfile, nchar(outfile)-3, nchar(outfile))))
                          }
            stopCluster(cl)

      }
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
