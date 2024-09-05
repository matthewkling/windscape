#' Pairwise wind distances between points
#'
#' Calculate pairwise wind cost-distances (e.g. travel times) or flow rates (the inverse of cost-distances)
#' among a set of sites, using the least cost path algorithm. This function is a wrapper around
#' \link[gdistance]{costDistance}.
#'
#' @param graph A \link{wind_graph}.
#' @param sites A two-column matrix of point coordinates.
#' @param adjust Whether to scale results to correct for discrepancies between point-to-point
#' distances and cell-to-cell distances. Default is TRUE.
#' @param rate Whether to return values as "rates" instead of the default "cost distances". Rates
#' are the inverse of cost distances, representing flow rather than travel time.
#' @return A square matrix of wind travel times.
#' @export
least_cost_distance <- function(graph, sites, adjust = TRUE, rate = FALSE){
      d <- gdistance::costDistance(graph, sites)
      if(adjust){
            r <- rast(nrows = graph@nrows, ncols = graph@ncols, extent = graph@extent, crs = graph@crs)
            d <- d * point_distance(sites) / cell_distance(r, sites)
      }
      if(rate){
            d <- 1 / d
      }
      d
}


#' Accumulated wind cost surface
#'
#' Calculate the accumulated wind cost-distance (e.g. travel times) or flow rate (the inverse of cost-distance)
#' from one or more sites to every grid cell across the domain, using the least cost path algorithm. This
#' function is a wrapper around \link[gdistance]{accCost}.
#'
#' @param graph A \link{wind_graph}.
#' @param sites A two-column matrix of point coordinates.
#' @param rate Whether to return values as "rates" instead of the default "cost distances". Rates
#' are the inverse of cost distances, representing flow rather than travel time.
#' @return A SpatRaster of wind connectivity values.
#' @export
least_cost_surface <- function(graph, sites, rate = FALSE){
      if(inherits(sites, "SpatVector")) sites <- crds(sites)
      d <- gdistance::accCost(graph, sites)
      if(rate){
            d <- 1 / d
      }
      rast(d)
}
