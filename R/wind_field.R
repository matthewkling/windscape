#' An S4 object class representing a wind field time series
#'
setClass("wind_series",
         contains = "SpatRaster",
         slots = c(n_steps = "numeric"))

#' Generate a wind field time series data set from a set of rasters.
#'
#' @param x Multi-layer \code{SpatRaster} with layers containing u and v wind components, or
#'    an object like a file path that can be converted to a \code{SpatRaster}. Data must be
#'    on a longitude/latitude grid with square cells (equal x and y resolution in degrees),
#'    with u and v components oriented to true east and north; see \link{wind_rose}.
#' @param order Either \code{"uuvv"}, the default, indicating `x` has all u components
#'    followed by all v components, or \code{"uvuv"}, indicating the u and v components
#'    of `x` are alternating.
#' @return A `wind_series` object, which is a particular form of \code{SpatRaster}.
#' @export
wind_series <- function(x, order = c("uuvv", "uvuv")){

      order <- match.arg(order)
      if(!inherits(x, "SpatRaster")) x <- rast(x)
      check_grid(x)
      if(terra::nlyr(x) %% 2 != 0) stop("if `v` is not specified, `u` must have an even number of layers.")

      # collate layers
      if(order == "uvuv"){
            even <- function(x) x %% 2 == 0
            x <- x[[c(which(!even(1:nlyr(x))),
                      which(even(1:nlyr(x))))]]
      }

      # create wind_field
      y <- as(as(x, "SpatRaster"), "wind_series")
      y@n_steps <- nlyr(x)/2
      y
}



#' An S4 object class representing a wind field
#'
setClass("wind_field",
         contains = "SpatRaster",
         slots = c(n_steps = "numeric"))


#' Generate a wind field data set from a pair of rasters.
#'
#' @param x A `SpatRaster` with two layers representing u and v wind components;
#'    note that these must be in lat-long coordinates.
#' @return A `wind_field` object, which is a particular form of `SpatRaster`.
#' @export
wind_field <- function(x){
      if(terra::nlyr(x) != 2) stop("`x` must have two layers.")
      xt <- ext(x)
      if(any(c(xt$xmin < -180, xt$xmax > 360, xt$ymin < -90, xt$ymax > 90))) stop("wind field rasters must be in lon-lat coordinates")
      as(as(x, "SpatRaster"), "wind_field") # via SpatRaster, so subclasses (e.g. wind_series layers) work
}


#' Load a wind_series from one or more raster files on disk
#'
#' Reads files in `wind_series` layout (all u layers followed by all v layers), such as the
#' monthly files saved by [ncar_download()], and combines them into a single `wind_series`.
#' Data are not loaded into memory until needed.
#'
#' @param x Character vector of file paths. Each file must hold an equal number of u and v
#'   layers, with u layers first, and all files must share the same grid. Time steps are combined
#'   in the order the files are given.
#' @return A `wind_series` object whose time steps are those of all the files combined: all
#'   files' u layers, followed by all files' v layers.
#' @seealso [ncar_download()]
#' @export
read_wind_series <- function(x){
      if(!is.character(x) || length(x) == 0) stop("`x` must be a character vector of file paths")
      missing <- !file.exists(x)
      if(any(missing)) stop("file(s) not found: ", paste(x[missing], collapse = ", "))
      r <- lapply(x, terra::rast)
      for(i in seq_along(r)){
            if(terra::nlyr(r[[i]]) %% 2 != 0) stop("file has an odd number of layers: ", x[i])
            if(i > 1 && !terra::compareGeom(r[[1]], r[[i]], stopOnError = FALSE))
                  stop("files are on different grids: ", x[1], " and ", x[i])
      }
      half <- function(z, which){
            n <- terra::nlyr(z) / 2
            z[[if(which == "u") seq_len(n) else n + seq_len(n)]]
      }
      wind_series(c(do.call(c, lapply(r, half, "u")), do.call(c, lapply(r, half, "v"))))
}
