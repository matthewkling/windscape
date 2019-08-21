
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' Calculate the direction of a vector based on its x and y components.
#'
#' @param x horizontal component (numeric vector)
#' @param y vertical component (numeric vector)
#' @return the angle of the vector, in degrees clockwise from 12 o'clock
direction <- function(x, y) atan2(x, y) * 180 / pi


#' Rotate a set of bearings counterclockwise by 90 degrees
#'
#' @param x vector of angles, in degrees
#' @return rotated angles
spin90 <- function(x){
      x <- x - 90
      x[x<(-180)] <- x[x<(-180)] + 360
      x
}



#' Augment a raster object, adding a layer containing the cell resolution.
#'
#' @param x raster layer or stack
#' @return the input layer, with an additional layer of resolution values
add_res <- function(x){
      r <- x[[1]]
      r[] <- mean(res(r[[1]]))
      stack(r, x)
}


#' Augment a raster object, adding a layer containing the latitude of each cell.
#'
#' @param x raster layer or stack
#' @return the input layer, with an additional layer of latitude values
add_lat <- function(x){
      lat <- x[[1]]
      lat[] <- coordinates(lat)[,2]
      stack(lat, x)
}



#' Calculate 8-neighbor edge loadings from a time series of u and v windspeeds
#'
#' @param x A vector of wind data containing: latitude, resolution, u
#'   windspeeds, v windspeeds
#' @param p A positive number indicating the power to raise windspeeds to (see details)
#' @return A vector of 8 conductance values to neighboring cells, clockwise
#'   starting with the southwest neighbor. If input windspeeds are in m/s,
#'   values are in (1 / seconds ^ p)
#' @details A value of p = 0 will ignore speed, assigning weights based on
#'   direction only. P = 1 assumes conductance is proportional to windspeed, p =
#'   2 assumes it's proportional to aerodynamic drag, and p = 3 assumes it's
#'   proportional to force. Any intermediate value can also be used.
windrose <- function(x, p=1){

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
      l <- edge_loadings(dir, weight, nb)

      # adjust conductance weights to account for inter-cell distance
      nd <- geosphere::distGeo(c(0, lat), nc)
      l <- l/nd
      # units are now in 1/seconds^p

      # reorder, clockwise from SW
      l[c(6:8, 1:5)]
}




# generate windoses for a raster dataset
# must be in lat-long proection

#' Generate a
#'
#' @param w Raster stack of u and v wind components; note that these must be in
#'   lat-long coordinates.
#' @param outfile File path to save output raster to.
#' @param order Either "uvuv" (the default) or "uuvv", indicating whether the u
#'   and v components of the w object are alternating.
#' @param ncores Integer indicating the number of computing cores to use for
#'   parallel processing. If multiple cores are used, an output file will be
#'   saved for each core, appending the core number to the user-supplied
#'   filename.
#' @param ... Additional arguments passed to `windrose`
#'
#' @return An 8-layer raster stack, where each layer is wind conductance from
#'   the focal cell to one of its neighbors (clockwise starting in the SW).
windrose_rasters <- function(w, outfile, order = "uvuv",  ncores = 1, ...){

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


#' Get neighbor names for a windrose object
#'
#' @return Names of neighbors; note that actual directions vary by latitude
windrose_names <- function() c("SW", "W", "NW", "N", "NE", "E", "SE", "S")

#' Get neighbor directions for a windrose object
#'
#' @return Bearings to neighbors; note that actual directions vary by latitude
#'   so these values are not accurate
windrose_bearings <- function() c(225, 270, 315, 0, 45, 90, 135, 180)




