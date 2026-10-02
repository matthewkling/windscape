#' Calculate the direction of a vector based on its x and y components.
#'
#' @param x horizontal component (numeric vector)
#' @param y vertical component (numeric vector)
#' @return the angle of the vector, in degrees clockwise from 12 o'clock
#' @noRd
direction <- function(x, y) atan2(x, y) * 180 / pi


#' Rotate a set of bearings counterclockwise by 90 degrees
#'
#' @param x vector of angles, in degrees
#' @return rotated angles
#' @noRd
spin90 <- function(x){
      x <- x - 90
      x[x<(-180)] <- x[x<(-180)] + 360
      x
}

#' Augment a raster object, adding layers containing cell row-column indices.
#'
#' @param x SpatRaster
#' @return x, with additional row and column index layers added
#' @noRd
add_coords <- function(windrose){
      rows <- cols <- windrose[[1]]
      rows[] <- rep(1:nrow(rows), each=ncol(rows))
      cols[] <- rep(1:ncol(rows), nrow(rows))
      windrose <- c(windrose, rows, cols)
      names(windrose) <- c("SW", "W", "NW", "N", "NE", "E", "SE", "S", "row", "col")
      return(windrose)
}


#' Extend a global raster with copies of itself
#'
#' @param x SpatRaster
#' @param width Longitudinal degrees of extension
#' @return a wider version of x, with the eastern side repeated on the west and vice-versa
#' @export
tesselate <- function(x, width = 20){
      x <- x %>% crop(ext(180 - width, 180, -90, 90)) %>% terra::shift(-360) %>% terra::merge(x)
      x <- x %>% crop(ext(-180, -180 + width, -90, 90)) %>% terra::shift(360) %>% terra::merge(x)
      x
}


#' Aspect ratio of grid cell at a given latitude
#'
#' @param y Numeric vector of latitudes
#' @param res Cell resolution, in degrees latitude; any sufficiently small value should yield reasonably precise results.
#' @noRd
aspect <- function(y, res = .00001){
      y1 <- y * pi / 180
      y2 <- (y + res) * pi / 180
      x <- res * pi / 180
      dy <- asin(sqrt((1 - cos(y2 - y1)) / 2))
      dx <- asin(sqrt((cos(y1) * cos(y2) * (1 - cos(x))) / 2))
      dy / dx

      # # confirmation that formulation is correct (the below quantities should be equal):
      # y = 70 # arbitrary
      # aspect(y)
      # geosphere::distHaversine(c(0, y), c(0, y + .0001)) / geosphere::distHaversine(c(0, y), c(0.0001, y))
}


# windscape's neighbor geometry (see rose()) assumes a longitude/latitude grid with square cells.
# A raster with no CRS is accepted if its extent is plausible for longitude/latitude.
check_grid <- function(x){
      if(!isTRUE(terra::is.lonlat(x, perhaps = TRUE, warn = FALSE)))
            stop("windscape currently supports only longitude/latitude grids, but `x` has a ",
                 "projected CRS (or an unknown CRS with a non-geographic extent). Reproject the ",
                 "data to longitude/latitude first; if u and v are defined relative to the ",
                 "projected grid, they must also be rotated to true east and north.", call. = FALSE)
      # rose() uses a single cell size, mean(res), for both axes. A 1% mismatch gives at most ~0.5%
      # error in neighbor distances; this tolerance admits nearly-square grids such as CFSR's
      # Gaussian grid (0.316 x 0.317 degrees).
      rs <- terra::res(x)
      if(abs(rs[1] - rs[2]) > 0.01 * max(rs))
            stop("grid cells must be square in degrees (within 1%), but `x` has resolution ",
                 signif(rs[1], 4), " x ", signif(rs[2], 4), ". Resample to equal x and y resolution.",
                 call. = FALSE)
      invisible(TRUE)
}
