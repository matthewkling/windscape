#' Generate a windrose raster stack from a stack of windfield rasters
#'
#' @param w Raster stack of u and v wind components; note that these must be in
#'   lat-long coordinates.
#' @param order Either "uvuv" (the default) or "uuvv", indicating whether the u
#'   and v components of the w object are alternating.
#' @param ncores Integer indicating the number of computing cores to use for
#'   parallel processing. If multiple cores are used, an output file will be
#'   saved for each core, appending the core number to the user-supplied
#'   filename.
#' @param p A positive number indicating the power to raise windspeeds to (see \code{\link{windrose} for details}).
#' @param ... Additional arguments passed to `raster::calc`, e.g. 'filename'.
#'
#' @return An 8-layer raster stack, where each layer is wind conductance from
#'   the focal cell to one of its neighbors (clockwise starting in the SW).
windrose_rasters <- function(w, order = "uvuv",  ncores = 1, p, ...){

      # collate data
      if(order == "uvuv"){
            even <- function(x) x %% 2 == 0
            w <- w[[c(which(!even(1:nlayers(w))),
                      which(even(1:nlayers(w))))]]
      }


      rosefun <- function(x) windrose(x, p)

      if(ncores == 1){
            wr <- add_res(w)
            wr <- add_lat(wr)
            wr <- raster::calc(wr, fun=rosefun, forceapply=TRUE, ...)
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
