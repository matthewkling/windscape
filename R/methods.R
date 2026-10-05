#' @include wind_field.R wind_rose.R
NULL

#' Selecting and combining layers of windscape objects
#'
#' `wind_series`, `wind_field`, and `wind_rose` objects are `SpatRaster`s whose layers have a
#' fixed structure (u and v layers, or conductance toward eight neighbors). Selecting layers
#' with `[[` or `terra::subset()`, or combining objects with `c()`, generally breaks that
#' structure, so these operations return a plain `SpatRaster`. To select time steps while keeping
#' a series intact, use [subset_series()]; to take one time step as a field, use
#' [wind_field()]; and to combine series, use [wind_series()] with a list of series or files.
#'
#' @param x A `wind_series`, `wind_field`, or `wind_rose`.
#' @param i,subset Layers to select, as for a `SpatRaster`.
#' @param j Not used.
#' @param ... Other arguments passed to the `SpatRaster` method, or for `c()`, other objects
#'    to combine.
#' @return A `SpatRaster`.
#' @examples
#' series <- windscape_example("wind_series")
#' class(series[[1:2]]) # a plain SpatRaster
#' class(subset_series(series, steps = 1:2)) # a wind_series
#' @name windscape-layers
#' @exportMethod [[ subset c
#' @aliases [[,wind_series-method [[,wind_field-method [[,wind_rose-method subset,wind_series-method subset,wind_field-method subset,wind_rose-method c,wind_series-method c,wind_field-method c,wind_rose-method
NULL

plain <- function(x) if(inherits(x, "SpatRaster")) as(x, "SpatRaster") else x

for(cls in c("wind_series", "wind_field", "wind_rose")){
      setMethod("[[", cls, function(x, i, j, ...) as(x, "SpatRaster")[[i]])
      setMethod("subset", cls, function(x, subset, ...) terra::subset(as(x, "SpatRaster"), subset, ...))
      setMethod("c", cls, function(x, ...) do.call(c, c(list(as(x, "SpatRaster")), lapply(list(...), plain))))
}
rm(cls)


#' Changing the grid of a wind rose
#'
#' A wind rose's conductances are rates of flow between neighboring cells, so they depend on cell
#' size. Changing the grid with `terra::aggregate()`, `terra::disagg()`, `terra::resample()`, or
#' `terra::project()` would average or interpolate conductances without accounting for the new
#' cell size, giving wrong results, so these operations are not allowed on a `wind_rose`. To
#' refine a rose's grid, use [downscale()]. To use a coarser or different grid, change the grid of
#' the `wind_series` (e.g. with `terra::aggregate()`) and build a new rose from it.
#'
#' @param x A `wind_rose`.
#' @param y,... Not used.
#' @return These methods signal an error.
#' @name wind_rose-grid
#' @aliases aggregate,wind_rose-method disagg,wind_rose-method resample,wind_rose-method project,wind_rose-method
#' @exportMethod aggregate disagg resample project
NULL

rose_grid_error <- function(fun){
      function(x, ...) stop(fun, "() would give wrong conductances for a wind_rose, because ",
                            "conductance depends on cell size. To refine the grid, use downscale(); for a ",
                            "coarser or different grid, change the grid of the wind_series and build a new ",
                            "rose from it.", call. = FALSE)
}
setMethod("aggregate", "wind_rose", rose_grid_error("aggregate"))
setMethod("disagg", "wind_rose", rose_grid_error("disagg"))
setMethod("resample", signature(x = "wind_rose", y = "ANY"), rose_grid_error("resample"))
setMethod("project", "wind_rose", rose_grid_error("project"))
