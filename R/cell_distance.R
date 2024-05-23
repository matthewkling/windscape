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
      cc <- crds(xx)[as.vector(extract(xx, ll))[[1]],] # cell centers
      geosphere::distm(cc) / 1000
}
