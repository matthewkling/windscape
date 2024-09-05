#' Pairwise distances between points
#'
#' A wrapper around \link[geosphere]{distm}, returning distances in km.
#'
#' @param x two-column matrix of point coordinates
#' @return Pairwise distances between the points, in km
#' @export
point_distance <- function(ll){
      require(geosphere)
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
      xx <- x[[1]][[1]] %>% crop(extent(ll)*2) %>% setValues(1:ncell(.))
      cc <- crds(xx)[as.vector(terra::extract(xx, ll))[[1]],] # cell centers
      point_distance(cc)
}


#' Compare point distances to cell distances
#'
#' Check for discrepancies between pairwise distances between a set of points, versus pairwise distances between the
#' centers of the raster grid cells that the points fall within. There will always be some differences due to the
#' continuous versus discrete nature of the data types, but if discrepancies are large as a percentage of the distances,
#' then they will contribute substantial noise that could invalidate wind connectivity estimates. This function reports
#' how many site pairs are in the same grid cell (an infinite difference that makes wind calculations impossible) as
#' well as information about the distribution of discrepancies as percentages of distance. The more site pairs there
#' are with nontrivial discrepancies, the less reliable a wind connectivity model will be for these sites. Problems
#' identified in this check can be solved using the \link{downscale} or \link{vrcd} functions.
#'
#' @param x SpatRaster (e.g. a `wind_rose`)
#' @param ll two-column matrix of site coordinates
#' @return A matrix of the ratios of point distances to cell distances
#' @export
check_cell_distance <- function(x, ll, return = FALSE){
      cell <- cell_distance(x, ll)
      r <- cell / point_distance(ll)
      d <- r[upper.tri(r)]
      d <- exp(abs(log(d))) - 1
      n <- sum(upper.tri(cell))
      f <- function(b) paste0("\n\t>= ", b*100, "%: ", sum(d >= b), " (", signif(mean(d >= b), 3), "%)")
      message("Total point pairs: ", n)
      message("Point pairs in the same grid cell: ",
              sum(cell[upper.tri(cell)] == 0),
              " (", signif(mean(cell[upper.tri(cell)] == 0) * 100, 3), "%)")
      m <- apply(matrix(c(0, .01, .01, .025, .025, .05, .05, .1, .1, .25, .25, Inf), ncol = 2, byrow = T),
            1, function(x){
                  b <- between(d, x[1], x[2])
                  paste0("\n\t", x[1]*100, "--", x[2]*100, "%: ", sum(b), " (", signif(mean(b)*100, 3), "%)")
            })
      message("Distribution of cell-point distance discrepancies:",
              paste(m))
      if(return) return(r)
}
