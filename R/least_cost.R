#' Pairwise least-cost travel times among sites
#'
#' Calculate pairwise wind cost-distances (e.g. travel times) or flow rates (the inverse of
#' cost-distances) among a set of sites, using the least cost path algorithm. For the random
#' walk counterpart, see [pairwise_random_walk()].
#'
#' Least-cost travel on a wind graph moves between neighboring cell centers, so a grid alone
#' places each site at the center of its cell. For sites a few cells apart or less, that
#' distorts both the distance and the direction between them, and sites in the same cell would
#' be zero hours apart. By default (`snap = FALSE`), sites are instead added to the graph at
#' their actual locations. Each site is linked to the centers of its own and the eight
#' surrounding cells, and directly to any other site in those cells, by edges whose travel time
#' is computed exactly for the wind where the edge starts, treated as uniform along the edge.
#' That wind is the cell's own at a cell center, and is interpolated between cell centers at a
#' site. Paths then run from site to site through these edges and the grid. Results change
#' continuously as sites move, rather than jumping at cell boundaries. With `snap = TRUE`,
#' sites are snapped to cell centers, as in [gdistance::costDistance()].
#'
#' The exact travel time through uniform wind is the continuum limit of least-cost travel on
#' an eight-neighbor grid: the time to cover a displacement using the cheapest combination of
#' the cell's eight flow vectors (conductance times the displacement to each neighbor; see
#' [net_flow()]), which uses at most two of them. Equivalently, the region reachable in one hour
#' is the convex hull of the flow vectors. Site edges therefore share the grid's metric: for
#' sites at cell centers, results are close to those from `snap = TRUE`, and slightly lower
#' where a site edge combines two directions more efficiently than the grid can.
#'
#' Paths are restricted to the eight neighbor directions, so cost distances are overestimated
#' for routes between neighbor bearings. For uniform wind on square cells the maximum error is
#' about 8 percent (routes 22.5 degrees off a grid axis). On a longitude/latitude grid, cells
#' narrow east-west toward the poles and neighbor bearings become uneven, so the maximum error
#' grows with latitude, to roughly 18 percent at 60 degrees. Site edges have the same
#' directional bias, so that it is consistent across distances.
#'
#' @param graph A \link{wind_graph}.
#' @param sites A two-column matrix of point coordinates.
#' @param snap Logical: snap sites to the centers of their grid cells? Default `FALSE`, which
#'    places sites at their actual locations within cells; see details.
#' @param rate Whether to return values as "rates" instead of the default "cost distances". Rates
#' are the inverse of cost distances, representing flow rather than travel time.
#' @return A square matrix with one row and column per site: the least-cost travel time (or,
#'    with `rate = TRUE`, its inverse) from the row's site to the column's site, in hours if
#'    `trans = 1` in [wind_rose()] and wind speeds are in m/s. Small travel times mean strong
#'    connectivity. Because wind graphs are directed, the matrix is generally asymmetric. The
#'    diagonal is zero. Sites outside the graph's extent get `NA`, with a warning.
#' @export
pairwise_least_cost <- function(graph, sites, snap = FALSE, rate = FALSE){
      if(inherits(sites, "SpatVector")) sites <- crds(sites)
      sites <- unname(as.matrix(sites))
      if(snap){
            d <- gdistance::costDistance(graph, sites)
      }else{
            d <- lc_sites(graph, sites)
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
#' @return A single-layer SpatRaster. With `rate = FALSE`, the layer is named `hours` and gives
#'   the travel cost from the nearest site to each cell, in hours if `trans = 1` in [wind_rose()]
#'   and wind speeds are in m/s (otherwise in relative units). With `rate = TRUE`, it is named
#'   `rate` and gives the inverse.
#' @export
least_cost_surface <- function(graph, sites, rate = FALSE){
      if(inherits(sites, "SpatVector")) sites <- crds(sites)
      d <- gdistance::accCost(graph, sites)
      if(rate){
            d <- 1 / d
      }
      d <- rast(d)
      names(d) <- if(rate) "rate" else "hours"
      d
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
#'    pair's [pairwise_least_cost()] (with `snap = TRUE`).
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
# ---- site handling for least-cost distances (snap = FALSE) ----

# Outgoing edge conductance from every cell toward each of its eight neighbors: an ncell x 8
# matrix, columns in ROSE_DIRS order, read from a transition matrix `tm` for an nr x nc grid
# (zero where there is no edge).
lc_out_conductance <- function(tm, nr, nc){
      cells <- seq_len(nr * nc)
      row <- (cells - 1) %/% nc + 1
      col <- (cells - 1) %% nc + 1
      dc <- c(-1, -1, -1, 0, 1, 1, 1, 0)
      dr <- c(1, 0, -1, -1, -1, 0, 1, 1) # rows increase southward
      out <- matrix(0, nr * nc, 8)
      for(k in 1:8){
            r2 <- row + dr[k]
            c2 <- col + dc[k]
            ok <- r2 >= 1 & r2 <= nr
            if(nc > 2) c2 <- (c2 - 1) %% nc + 1 else ok <- ok & c2 >= 1 & c2 <= nc # wrapped edges, if any
            out[ok, k] <- tm[cbind(cells[ok], (r2[ok] - 1) * nc + c2[ok])]
      }
      out
}

# East and north displacement (km on longitude/latitude grids, map units otherwise) from points
# `a` to points `b` (two-column matrices).
lc_displacement <- function(a, b, lonlat){
      if(!lonlat) return(b - a)
      d <- geosphere::distGeo(a, b) / 1000
      brg <- geosphere::bearingRhumb(a, b) * pi / 180
      brg[d == 0] <- 0 # bearing is undefined for coincident points
      cbind(d * sin(brg), d * cos(brg))
}

# Least-cost travel time across displacements D (m x 2) through locally uniform wind with flow
# vectors Vx, Vy (m x 8; conductance times neighbor displacement). This is the continuum limit
# of least-cost travel on an eight-neighbor grid: the minimum total time sum(t) over t >= 0 with
# sum_k t_k V_k = D, a linear program whose optimum uses at most two directions.
lc_leg_time <- function(D, Vx, Vy){
      m <- nrow(D)
      best <- rep(Inf, m)
      dn <- sqrt(rowSums(D^2))
      best[dn == 0] <- 0
      vn <- sqrt(Vx^2 + Vy^2)
      for(i in 1:8){ # single directions
            cr <- Vx[, i] * D[, 2] - Vy[, i] * D[, 1]
            dt <- Vx[, i] * D[, 1] + Vy[, i] * D[, 2]
            ok <- vn[, i] > 0 & dn > 0 & dt > 0 & abs(cr) <= 1e-9 * vn[, i] * dn
            best[ok] <- pmin(best[ok], dn[ok] / vn[ok, i])
      }
      for(i in 1:7) for(j in (i + 1):8){ # pairs of directions
            det <- Vx[, i] * Vy[, j] - Vx[, j] * Vy[, i]
            ok <- abs(det) > 1e-12 * vn[, i] * vn[, j] & vn[, i] > 0 & vn[, j] > 0
            ti <- (D[, 1] * Vy[, j] - Vx[, j] * D[, 2]) / det
            tj <- (Vx[, i] * D[, 2] - Vy[, i] * D[, 1]) / det
            ok <- ok & ti >= 0 & tj >= 0
            best[ok] <- pmin(best[ok], ti[ok] + tj[ok])
      }
      best
}

# Pairwise least-cost travel times with sites as nodes of their own: each site is joined to the
# centers of its own and neighboring cells, and directly to sites in its own and neighboring
# cells, by edges whose cost is the continuum travel time (lc_leg_time()) through the wind where
# the edge starts: at a cell center, that cell's flow vectors, so a leg along a neighbor direction
# costs at most the grid edge in that direction; at a site, flow vectors interpolated between cell
# centers. Returns an n x n matrix. Upwind graphs are the transpose of downwind graphs, so they are
# solved as downwind and transposed back.
lc_sites <- function(graph, sites){
      upwind <- identical(graph@direction, "upwind")
      tmat <- gdistance::transitionMatrix(graph)
      if(upwind) tmat <- Matrix::t(tmat)
      nr <- graph@nrows
      nc <- graph@ncols
      nn <- nr * nc
      ns <- nrow(sites)
      r <- rast(nrows = nr, ncols = nc, extent = graph@extent, crs = graph@crs)
      if(is.na(graph@crs)) terra::crs(r) <- "EPSG:4326" # wind graphs don't keep their CRS; windscape grids are lon/lat
      lonlat <- isTRUE(terra::is.lonlat(r, perhaps = TRUE, warn = FALSE))
      cond <- lc_out_conductance(tmat, nr, nc)
      D <- rw_cell_displacements(r)
      Vx <- cond * D$dx # flow vectors of each cell, km/h
      Vy <- cond * D$dy
      # flow vectors at sites: bilinear interpolation between cell centers, so they vary
      # continuously as a site moves (falling back to the site's own cell at the grid margin)
      site_flows <- function(xy, cell){
            vr <- terra::rast(r, nlyrs = 16)
            terra::values(vr) <- cbind(Vx, Vy)
            v <- terra::extract(vr, xy, method = "bilinear")
            v <- as.matrix(v[, setdiff(names(v), "ID"), drop = FALSE])
            bad <- !stats::complete.cases(v)
            v[bad, ] <- cbind(Vx, Vy)[cell[bad], , drop = FALSE]
            list(x = v[, 1:8, drop = FALSE], y = v[, 9:16, drop = FALSE])
      }

      site_cell <- terra::cellFromXY(r, sites)
      ok <- !is.na(site_cell)
      if(!all(ok)) warning("some sites are outside the graph's extent; their distances are NA")
      rc <- terra::rowColFromCell(r, site_cell)
      sv <- site_flows(sites[ok, , drop = FALSE], site_cell[ok])
      SVx <- matrix(NA_real_, ns, 8)
      SVy <- matrix(NA_real_, ns, 8)
      SVx[ok, ] <- sv$x
      SVy[ok, ] <- sv$y

      # site <-> cell-center legs, for the 3 x 3 block of cells around each site
      off <- expand.grid(dr = -1:1, dc = -1:1)
      legs <- do.call(rbind, lapply(which(ok), function(s){
            rr <- rc[s, 1] + off$dr
            cc <- rc[s, 2] + off$dc
            keep <- rr >= 1 & rr <= nr & cc >= 1 & cc <= nc
            data.frame(site = s, cell = (rr[keep] - 1) * nc + cc[keep])
      }))
      ctr <- terra::xyFromCell(r, legs$cell)
      sxy <- sites[legs$site, , drop = FALSE]
      t_out <- lc_leg_time(lc_displacement(sxy, ctr, lonlat), SVx[legs$site, , drop = FALSE],
                           SVy[legs$site, , drop = FALSE])
      t_in <- lc_leg_time(lc_displacement(ctr, sxy, lonlat), Vx[legs$cell, , drop = FALSE],
                          Vy[legs$cell, , drop = FALSE])

      # direct site-to-site legs, for sites in the same or neighboring cells
      near <- function(v) outer(v, v, function(a, b) abs(a - b) <= 1)
      pr <- which(outer(ok, ok, "&") & near(rc[, 1]) & near(rc[, 2]), arr.ind = TRUE)
      pr <- pr[pr[, 1] != pr[, 2], , drop = FALSE]
      t_dir <- lc_leg_time(lc_displacement(sites[pr[, 1], , drop = FALSE], sites[pr[, 2], , drop = FALSE], lonlat),
                           SVx[pr[, 1], , drop = FALSE], SVy[pr[, 1], , drop = FALSE])

      # augmented graph: cells 1..nn, then sites
      tm <- Matrix::summary(methods::as(tmat, "generalMatrix"))
      tm <- tm[tm$x > 0, ]
      from <- c(tm$i, nn + legs$site, legs$cell, nn + pr[, 1])
      to <- c(tm$j, legs$cell, nn + legs$site, nn + pr[, 2])
      w <- c(1 / tm$x, t_out, t_in, t_dir)
      keep <- is.finite(w)
      g <- igraph::make_graph(rbind(from[keep], to[keep]), n = nn + ns, directed = TRUE)
      d <- igraph::distances(g, v = nn + which(ok), to = nn + which(ok), mode = "out",
                             weights = w[keep], algorithm = "dijkstra")
      out <- matrix(NA_real_, ns, ns)
      out[ok, ok] <- d
      if(upwind) out <- t(out)
      out
}
