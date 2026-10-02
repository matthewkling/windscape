#' Pairwise wind distances between points
#'
#' Calculate pairwise wind cost-distances (e.g. travel times) or flow rates (the inverse of cost-distances)
#' among a set of sites, using the least cost path algorithm. This function is a wrapper around
#' \link[gdistance]{costDistance}.
#'
#' Paths are restricted to the eight neighbor directions, so cost distances are overestimated
#' for routes between neighbor bearings. For uniform wind on square cells the maximum error is
#' about 8 percent (routes 22.5 degrees off a grid axis). On a longitude/latitude grid, cells
#' narrow east-west toward the poles and neighbor bearings become uneven, so the maximum error
#' grows with latitude, to roughly 18 percent at 60 degrees. The \code{adjust} option corrects
#' for point-versus-cell distance discrepancies, not for this directional bias.
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
            pd <- point_distance(sites)
            ratio <- pd / cell_distance(r, sites)
            ratio[pd == 0] <- 1 # a site to itself (or an identical site): no adjustment
            d <- d * ratio
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
#' Paths are restricted to the eight neighbor directions, so cost distances are overestimated
#' for routes between neighbor bearings; see \link{least_cost_distance} for magnitudes.
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
