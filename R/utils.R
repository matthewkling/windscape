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


# Does a longitude/latitude grid span all 360 degrees of longitude (within half a cell), so
# that its east and west edges meet? Always FALSE for projected grids.
is_global <- function(x){
      if(!isTRUE(terra::is.lonlat(x, perhaps = TRUE, warn = FALSE))) return(FALSE)
      width <- terra::xmax(x) - terra::xmin(x)
      abs(width - 360) <= terra::res(x)[1] / 2
}

# Resolve a `wrap` argument to TRUE or FALSE for grid `x`. NULL (the default everywhere) wraps
# exactly when the grid is global. TRUE on a lon/lat grid that isn't global gets a warning,
# since joining its edges isn't physically meaningful; planar grids can wrap freely (e.g. for
# periodic test domains). Wrapping needs at least three columns: with fewer, a cell's east and
# west neighbors coincide.
resolve_wrap <- function(x, wrap = NULL){
      if(is.null(wrap)) return(is_global(x) && terra::ncol(x) >= 3)
      if(!(is.logical(wrap) && length(wrap) == 1 && !is.na(wrap)))
            stop("`wrap` must be TRUE, FALSE, or NULL (to wrap global grids)", call. = FALSE)
      if(!wrap) return(FALSE)
      if(terra::ncol(x) < 3) stop("`wrap = TRUE` requires a grid at least three cells wide", call. = FALSE)
      if(isTRUE(terra::is.lonlat(x, perhaps = TRUE, warn = FALSE)) && !is_global(x)){
            warning("`wrap = TRUE` joins the east and west edges of the grid, but it spans ",
                    signif(terra::xmax(x) - terra::xmin(x), 4), " degrees of longitude rather than 360",
                    call. = FALSE)
      }
      TRUE
}

# Trail numbers for path data whose paths may cross a wrapped grid's east-west seam: a new trail
# starts with each new path (a change in `id`) and wherever x jumps by more than half the grid's
# `width` between consecutive points, so that no line is drawn across the map. Rows must be
# ordered along each path.
seam_trails <- function(id, x, width){
      n <- length(x)
      if(n == 0) return(integer(0))
      new <- c(TRUE, id[-1] != id[-n] | abs(diff(x)) > width / 2)
      cumsum(new)
}


# Cell number of the neighbor at row offset `dr` and column offset `dc` from every cell of an
# nr x nc grid, in terra cell order; NA where the neighbor is off the grid. With `wrap`, columns
# wrap around, so the first and last columns are neighbors.
neighbor_cells <- function(nr, nc, dr, dc, wrap = FALSE){
      row <- rep(seq_len(nr), each = nc)
      col <- rep(seq_len(nc), times = nr)
      r2 <- row + dr
      c2 <- col + dc
      if(wrap) c2 <- (c2 - 1) %% nc + 1
      j <- (r2 - 1) * nc + c2
      j[r2 < 1 | r2 > nr | c2 < 1 | c2 > nc] <- NA
      j
}
