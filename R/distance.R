#' Pairwise distances between points
#'
#' A wrapper around \link[geosphere]{distm}, returning distances in km.
#'
#' @param ll two-column matrix of point coordinates
#' @return Pairwise distances between the points, in km
#' @export
point_distance <- function(ll){
      geosphere::distm(ll) / 1000
}


#' Pairwise distances between cell centroids
#'
#' Calculate the pairwise distances between the centroids of the raster grid cells that a set of points fall into.
#' This is useful as an alternative to the distances between the points themselves, if distances are going to be
#' compared to other metrics that are based on cell centers (e.g. wind connectivity).
#'
#' @param x SpatRaster
#' @param ll two-column matrix of site coordinates
#' @return Pairwise distances between the cell centroids, in km
#' @export
cell_distance <- function(x, ll){
      cc <- terra::xyFromCell(x, terra::cellFromXY(x, ll)) # cell centers
      point_distance(cc)
}


#' Check how grid cells distort distances among sites
#'
#' Compares pairwise distances among sites with distances between the centers of the grid cells
#' the sites fall in, and prints a report. Connectivity models that work from cell to cell treat
#' each site as the center of its cell: the random walk functions ([random_walk()],
#' [pairwise_random_walk()]), and the least-cost functions when sites are snapped to cell centers
#' ([least_cost_surface()], [least_cost_paths()], and `pairwise_least_cost(snap = TRUE)`). For
#' these, sites separated by only a few cells have distorted distances and directions, and sites in
#' the same cell can't be distinguished at all. [pairwise_least_cost()] with its default
#' `snap = FALSE` places sites at their actual locations, so this check does not apply to it.
#'
#' Where many site pairs are affected, options are to use wind data on a finer grid, or to
#' [downscale()] the wind rose (see its documentation for how downscaling changes random walk
#' results).
#'
#' @param x A SpatRaster on the grid used for connectivity modeling, e.g. a `wind_rose`.
#' @param ll A two-column matrix of site coordinates.
#' @param return Logical: return the matrix of ratios? Default `FALSE`, which only prints the
#'    report.
#' @return Prints the number of site pairs, the number (and percentage) in the same grid cell,
#'    and the distribution of discrepancies between cell and point distances, as percentages of
#'    point distance. If `return = TRUE`, also returns a matrix of the ratios of cell distances to
#'    point distances (`NaN` for a site with itself, and 0 for distinct sites in the same cell).
#' @export
check_cell_distance <- function(x, ll, return = FALSE){
      cell <- cell_distance(x, ll)
      r <- cell / point_distance(ll)
      d <- r[upper.tri(r)]
      d <- exp(abs(log(d))) - 1
      n <- sum(upper.tri(cell))
      f <- function(b) paste0("\n\t>= ", b*100, "%: ", sum(d >= b), " (", signif(mean(d >= b), 3), "%)")
      message("Total point pairs: ", n)
      cells <- terra::cellFromXY(x, ll)
      same <- outer(cells, cells, "==")[upper.tri(cell)] # compare cell IDs, not distances (which may not be exactly 0)
      message("Point pairs in the same grid cell: ",
              sum(same),
              " (", signif(mean(same) * 100, 3), "%)")
      m <- apply(matrix(c(0, .01, .01, .025, .025, .05, .05, .1, .1, .25, .25, Inf), ncol = 2, byrow = T),
                 1, function(x){
                       b <- d >= x[1] & d <= x[2]
                       paste0("\n\t", x[1]*100, "--", x[2]*100, "%: ", sum(b), " (", signif(mean(b)*100, 3), "%)")
                 })
      message("Distribution of cell-point distance discrepancies:",
              paste(m))
      if(return) return(r)
}
