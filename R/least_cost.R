#' Pairwise least-cost travel times among sites
#'
#' Calculate pairwise wind cost-distances (e.g. travel times) or flow rates (the inverse of
#' cost-distances) among a set of sites, using the least cost path algorithm. This function is a
#' wrapper around \link[gdistance]{costDistance}. For the random walk counterpart, see
#' [pairwise_random_walk()].
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
#' @return A square matrix with one row and column per site: the least-cost travel time (or,
#'    with `rate = TRUE`, its inverse) from the row's site to the column's site. Small travel
#'    times mean strong connectivity. Because wind graphs are directed, the matrix is generally
#'    asymmetric. The diagonal is zero.
#' @export
pairwise_least_cost <- function(graph, sites, adjust = TRUE, rate = FALSE){
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
#' for routes between neighbor bearings; see \link{pairwise_least_cost} for magnitudes.
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


#' Least-cost paths through a wind graph
#'
#' Traces the least-cost (fastest) paths between sites through a wind graph, as sequences of
#' grid cell centers. The result has the same structure as [wind_trails()] output, so it can be
#' drawn with [geom_wind_trail()].
#'
#' @param graph A `wind_graph`, created with [wind_graph()].
#' @param from,to Origin and destination sites: two-column matrices (or data frames) of
#'    longitude and latitude, or `SpatVector`s of points.
#' @param pairs Which origin-destination pairs to trace: `"all"` (the default) traces a path
#'    from every origin to every destination; `"nearest"` traces one path from each origin, to
#'    the destination it can reach at the lowest cost; `"matched"` pairs `from` and `to` row by
#'    row, which requires them to have the same number of sites.
#'
#' @return A data frame with one row per path vertex, ordered from origin to destination along
#'    each path:
#' * `trail`: path id.
#' * `from`, `to`: row numbers of the path's origin in `from` and destination in `to`.
#' * `step`: vertex number along the path, starting at 0 at the origin.
#' * `hours`: cumulative travel cost from the origin, in hours (if `trans = 1` in
#'    [wind_rose()] and wind speeds are in m/s). Its final value on each path equals the
#'    pair's [pairwise_least_cost()] (with `adjust = FALSE`).
#' * `x`, `y`: longitude and latitude of the grid cell center.
#'
#' Pairs with no path (where the destination can't be reached, e.g. because wind never blows
#' toward it) are omitted with a warning. Pairs whose origin and destination fall in the same
#' grid cell are omitted, since their path has no length.
#'
#' @details
#' Paths are computed with [gdistance::shortestPath()]. They follow the graph's direction:
#' with a downwind graph (the default in [wind_graph()]), a path from `from` to `to` is the
#' fastest downwind route. With an upwind graph, it is the fastest downwind route from `to` to
#' `from`, traced in reverse. Because wind graphs are directed, the path from A to B generally
#' differs from the path from B to A.
#'
#' Paths move between neighboring grid cells, in eight directions, so they appear as segments
#' at 45-degree angles, including occasional staircase patterns where the optimal route lies
#' between two neighbor directions. This reflects the model's structure rather than the plot.
#' Paths from one origin often share segments, showing the main transport corridors.
#'
#' @examples
#' rose <- windscape_example("wind_rose")
#' graph <- wind_graph(rose)
#' site <- cbind(-105, 40)
#' destinations <- cbind(c(-95, -115, -100, -110), c(45, 35, 32, 48))
#' paths <- least_cost_paths(graph, site, destinations)
#' head(paths)
#'
#' library(ggplot2)
#' ggplot(paths, aes(x, y)) +
#'   geom_wind_trail(aes(color = hours)) +
#'   coord_quickmap()
#' @export
least_cost_paths <- function(graph, from, to, pairs = c("all", "nearest", "matched")){
      if(!inherits(graph, "TransitionLayer")) stop("`graph` must be a wind_graph")
      pairs <- match.arg(pairs)
      as_sites <- function(s, name){
            if(inherits(s, "SpatVector")) s <- terra::crds(s)
            s <- as.matrix(s)
            if(ncol(s) != 2 || !is.numeric(s)) stop("`", name, "` must be a two-column matrix of coordinates")
            unname(s)
      }
      from <- as_sites(from, "from")
      to <- as_sites(to, "to")

      # origin-destination pairs to trace
      cost <- gdistance::costDistance(graph, from, to)
      cost <- matrix(cost, nrow(from), nrow(to))
      if(pairs == "all"){
            od <- expand.grid(to = seq_len(nrow(to)), from = seq_len(nrow(from)))[, c("from", "to")]
      }else if(pairs == "nearest"){
            best <- apply(cost, 1, function(z) if(all(!is.finite(z))) NA else which.min(z))
            od <- data.frame(from = seq_len(nrow(from)), to = best)
            od <- od[!is.na(od$to), ]
      }else{
            if(nrow(from) != nrow(to)) stop("with `pairs = \"matched\"`, `from` and `to` must have the same number of sites")
            od <- data.frame(from = seq_len(nrow(from)), to = seq_len(nrow(to)))
      }
      template <- raster::raster(graph)
      same_cell <- raster::cellFromXY(template, from[od$from, , drop = FALSE]) ==
            raster::cellFromXY(template, to[od$to, , drop = FALSE])
      od <- od[!same_cell, , drop = FALSE]
      unreachable <- !is.finite(cost[cbind(od$from, od$to)])
      if(any(unreachable)) warning(sum(unreachable), " origin-destination pair(s) have no path and were omitted")
      od <- od[!unreachable, , drop = FALSE]
      if(nrow(od) == 0){
            return(data.frame(trail = integer(0), from = integer(0), to = integer(0),
                              step = integer(0), hours = numeric(0), x = numeric(0), y = numeric(0)))
      }

      tm <- gdistance::transitionMatrix(graph)
      out <- list()
      for(i in unique(od$from)){
            dest <- od$to[od$from == i]
            sl <- gdistance::shortestPath(graph, from[i, , drop = FALSE], to[dest, , drop = FALSE],
                                          output = "SpatialLines")
            for(k in seq_along(dest)){
                  xy <- sl@lines[[k]]@Lines[[1]]@coords
                  cells <- raster::cellFromXY(template, xy)
                  step_cost <- 1 / tm[cbind(cells[-length(cells)], cells[-1])]
                  out[[length(out) + 1]] <- data.frame(from = i, to = dest[k],
                                                       step = seq_len(nrow(xy)) - 1,
                                                       hours = c(0, cumsum(step_cost)),
                                                       x = xy[, 1], y = xy[, 2])
            }
      }
      out <- do.call(rbind, out)
      out$trail <- as.integer(factor(paste(out$from, out$to), levels = unique(paste(out$from, out$to))))
      out <- out[, c("trail", "from", "to", "step", "hours", "x", "y")]
      rownames(out) <- NULL
      out
}
