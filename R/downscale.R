#' Downscale a wind rose to higher spatial resolution
#'
#' Disaggregates a wind rose to a finer grid with [terra::disagg()], interpolating between cells,
#' and multiplies conductance by `fact`. Conductance is a rate per hop between neighboring cells,
#' and hops on the finer grid are `fact` times shorter, so this rescaling keeps flow (conductance
#' times the distance to each neighbor) unchanged: wind moves material at the same speed on the
#' downscaled grid as on the original. Downscaling adds no wind information. It smooths the wind
#' field between the original cell centers, but does not represent real fine-scale variation
#' within cells, such as that caused by terrain.
#'
#' @param x A `wind_rose`.
#' @param fact A single integer: the factor by which to increase resolution. Each cell becomes
#'    `fact` x `fact` cells.
#' @param method Interpolation method, passed to [terra::disagg()]: `"bilinear"` (the default)
#'    interpolates between the original cell centers; `"near"` gives each new cell the value of
#'    the cell it came from.
#' @return A `wind_rose` with `fact` times as many rows and columns as `x`, and the same `trans`
#'    and `n_steps`. Memory use and run time for connectivity analyses grow with the number of
#'    cells, roughly `fact^2`.
#' @details
#' How much downscaling changes results, and whether it improves them, depends on the
#' connectivity model.
#'
#' **Least-cost analyses.** [pairwise_least_cost()] places sites at their actual locations
#' within cells by default (`snap = FALSE`), so it handles nearby sites without downscaling, and
#' downscaling changes its results little: typically by a few percent, reflecting the smoother
#' interpolated wind field. Downscaling does reduce error for nearby sites when sites are
#' snapped to cell centers, as in `pairwise_least_cost(snap = TRUE)`, [least_cost_surface()], and
#' [least_cost_paths()].
#'
#' **Random walk analyses.** The random walk functions treat each site as its cell, so downscaling
#' separates nearby sites into distinct cells (see [check_cell_distance()]). But the spread of a
#' random walk is numerical diffusion that scales with cell size (see [random_walk()]), so
#' downscaling also narrows the dispersal kernel: in one test, downscaling by a factor of 4 reduced
#' the spread of a stream-mode walk by 30 to 50 percent. Downscaling therefore changes the model
#' as well as the grid. Results at different resolutions are not directly comparable, so choose one
#' resolution for an analysis and treat it as part of the model specification.
#' @examples
#' rose <- windscape_example("wind_rose")
#' fine <- downscale(rose, 2)
#' dim(rose)
#' dim(fine)
#' @export
downscale <- function(x, fact, method = "bilinear"){
      if(!inherits(x, "wind_rose")) stop("x must be a wind_rose object")
      if(length(fact) != 1) stop("fact must be a single integer")
      disagg(x, fact, method = method) * fact
}
