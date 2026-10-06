#' An S4 object class representing a wind graph
#'
setClass("wind_graph",
         contains = "TransitionLayer",
         slots = c(p = "numeric",
                   direction = "character"))

#' Build a wind connectivity graph
#'
#' Constructs the directed network that the least-cost functions ([least_cost()],
#' [least_cost_paths()], and [pairwise_least_cost()]) find paths through. Those functions build
#' the graph from a wind rose automatically, so you don't normally need to call this function.
#' Building the graph yourself can save time when making many least-cost calls on a large grid:
#' pass the graph in place of the rose.
#'
#' The graph links each grid cell to its eight neighbors. Each link's conductance is the mean,
#' over its two cells, of the rose's conductance in the link's direction. In a downwind graph,
#' links point the way the wind carries material; an upwind graph is the same network with every
#' link reversed, for measuring travel toward a site.
#'
#' @param x A `wind_rose`.
#' @param direction Either `"downwind"` (the default) or `"upwind"`: whether links follow the
#'    wind, or run against it.
#' @param wrap Logical: join the east and west edges of the grid, linking cells across them? The
#'    default, `NULL`, does so if `x` is a global grid spanning all 360 degrees of longitude,
#'    where -180 and 180 are the same meridian, and not otherwise. `TRUE` or `FALSE` overrides
#'    this; `TRUE` on a longitude/latitude grid that isn't global gives a warning.
#' @return A `wind_graph`: a gdistance \link[gdistance]{Transition-class} object, with its
#'    direction recorded.
#' @examples
#' rose <- windscape_example("wind_rose")
#' graph <- wind_graph(rose)
#' sites <- cbind(c(-110, -100), c(40, 40))
#' pairwise_least_cost(graph, sites)
#' @export
#' @aliases build_wind_graph
wind_graph <- function(x, direction = "downwind", wrap = NULL){
      if(!inherits(x, "wind_rose")) stop("`x` must be a wind_rose")
      direction <- match.arg(direction, c("downwind", "upwind"))
      wrap <- resolve_wrap(x, wrap)
      g <- as(wind_transition(x, direction = direction, wrap = wrap), "wind_graph")
      g@p <- g@p
      g@direction <- direction
      g
}


#' geoCorrection is not applicable to wind graphs
#'
#' Wind graph conductances already account for the distances between cell centers (see
#' \link{wind_rose}), so \link[gdistance]{geoCorrection} must not be applied to them.
#'
#' @param x A \code{wind_graph}.
#' @param type,... Ignored.
#' @return Always an error.
#' @export
setMethod("geoCorrection", "wind_graph", function(x, type, ...){
      stop("`geoCorrection()` should not be used on `wind_graph` objects, because their ",
           "conductances already account for distances between cells.", call. = FALSE)
})


# Build the transition layer of a wind graph from a wind rose, vectorized over edges. Edges
# connect each non-NA cell to its eight neighbors (plus left-right "suture" edges between the
# first and last columns, if `wrap`). Each edge's conductance is the mean, over the edge's two
# cells, of the rose's conductance in the edge's direction (or, for an upwind graph, in the
# opposite direction).
wind_transition <- function(x, direction = "downwind", wrap = FALSE){
      v <- terra::values(x) # cells x 8 layers: SW, W, NW, N, NE, E, SE, S
      nr <- terra::nrow(x)
      nc <- terra::ncol(x)
      n <- nr * nc
      e <- terra::ext(x)
      template <- raster::raster(nrows = nr, ncols = nc, crs = terra::crs(x),
                                 ext = raster::extent(e$xmin, e$xmax, e$ymin, e$ymax)) # geometry only

      tr <- new("TransitionLayer",
                nrows = as.integer(nr), ncols = as.integer(nc),
                extent = raster::extent(template), crs = raster::projection(template, asText = FALSE),
                transitionMatrix = Matrix::Matrix(0, n, n), transitionCells = 1:n)

      # edges from each non-NA cell to its eight non-NA neighbors, built directly rather than with
      # raster::adjacent(), which joins the east and west edges of global grids on its own
      cells <- which(!is.na(v[, 1]))
      dr <- c(1, 0, -1, -1, -1, 0, 1, 1)
      dc <- c(-1, -1, -1, 0, 1, 1, 1, 0)
      adj <- do.call(rbind, lapply(1:8, function(k){
            j <- neighbor_cells(nr, nc, dr[k], dc[k], wrap)[cells]
            ok <- !is.na(j) & !is.na(v[j, 1])
            cbind(cells[ok], j[ok])
      }))

      from <- adj[, 1]
      to <- adj[, 2]
      rf <- (from - 1) %/% nc + 1; cf <- (from - 1) %% nc + 1
      rt <- (to - 1) %/% nc + 1;   ct <- (to - 1) %% nc + 1

      # direction of each edge; suture edges (columns not adjacent) have reversed east-west sense
      suture <- abs(cf - ct) != 1
      west <- ifelse(suture, ct > cf, ct < cf)
      east <- ifelse(suture, ct < cf, ct > cf)
      south <- rf < rt
      north <- rt < rf
      k <- ifelse(south & west, 1, ifelse(!south & !north & west, 2, ifelse(north & west, 3,
           ifelse(north & !west & !east, 4, ifelse(north & east, 5, ifelse(!south & !north & east, 6,
           ifelse(south & east, 7, 8)))))))
      if(direction == "upwind") k <- c(5:8, 1:4)[k]

      # mean of the two cells' conductance in the edge's direction
      cond <- (v[cbind(from, k)] + v[cbind(to, k)]) / 2

      # zero-conductance edges (wind never blows that way) are left out, not stored as explicit zeros
      gdistance::transitionMatrix(tr) <- Matrix::drop0(Matrix::sparseMatrix(i = from, j = to, x = cond,
                                                                             dims = c(n, n)))
      gdistance::matrixValues(tr) <- "resistance"
      tr
}
