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

#' Augment a raster object, adding layers containing cell row-column indices.
#'
#' @param x SpatRaster
#' @return x, with additional row and column index layers added
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
      x <- x %>% crop(extent(180 - width, 180, -90, 90)) %>% terra::shift(-360) %>% terra::merge(x)
      x <- x %>% crop(extent(-180, -180 + width, -90, 90)) %>% terra::shift(360) %>% terra::merge(x)
      x
}
