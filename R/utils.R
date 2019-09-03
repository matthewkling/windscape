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


#' Get neighbor names for a windrose object
#'
#' @return Names of neighbors; note that actual directions vary by latitude
windrose_names <- function() c("SW", "W", "NW", "N", "NE", "E", "SE", "S")


#' Get neighbor directions for a windrose object
#'
#' @return Bearings to neighbors; note that actual directions vary by latitude
#'   so these values are not accurate
windrose_bearings <- function() c(225, 270, 315, 0, 45, 90, 135, 180)
