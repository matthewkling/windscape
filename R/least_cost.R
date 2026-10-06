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
#' Unlike [least_cost()] and [least_cost_paths()], this function has no `direction` argument:
#' the matrix holds travel in both directions between every pair of sites. Row `i` gives travel
#' downwind from site `i`, and column `j` gives travel upwind of site `j`.
#'
#' @param rose A [wind_rose()]. The least-cost functions build a [wind_graph()] from it
#'    internally. Alternatively, a `wind_graph` built in advance, to save rebuilding it when
#'    making many calls on a large grid; its direction must match the analysis (always
#'    `"downwind"` here).
#' @param sites A two-column matrix (or data frame) of point coordinates, or a `SpatVector` of
#'    points.
#' @param snap Logical: snap sites to the centers of their grid cells? Default `FALSE`, which
#'    places sites at their actual locations within cells; see details.
#' @param rate Whether to return values as "rates" instead of the default "cost distances". Rates
#' are the inverse of cost distances, representing flow rather than travel time.
#' @param ... Further arguments passed to [wind_graph()], such as `wrap`. Not allowed when
#'    `rose` is already a `wind_graph`.
#' @return A square matrix with one row and column per site: the least-cost travel time (or,
#'    with `rate = TRUE`, its inverse) from the row's site to the column's site, in hours if
#'    `trans = 1` in [wind_rose()] and wind speeds are in m/s. Small travel times mean strong
#'    connectivity. Because wind is directional, the matrix is generally asymmetric. The
#'    diagonal is zero. Sites outside the grid get `NA`, with a warning.
#' @examples
#' rose <- windscape_example("wind_rose")
#' sites <- cbind(lon = c(-110, -105, -100, -95), lat = c(40, 42, 38, 44))
#' pairwise_least_cost(rose, sites) # travel time from row site to column site, in hours
#' @export
pairwise_least_cost <- function(rose, sites, snap = FALSE, rate = FALSE, ...){
      if(inherits(rose, "wind_graph") && identical(rose@direction, "upwind")){
            stop("pairwise_least_cost() needs a downwind wind_graph. Its result already covers ",
                 "both directions: element [i, j] is travel from site i to site j, so travel ",
                 "toward a site is read down its column.", call. = FALSE)
      }
      graph <- lc_graph(rose, "downwind", ...)
      sites <- as_site_matrix(sites, "sites")
      if(snap){
            ok <- lc_inside(graph, sites)
            if(!all(ok)) warning("some sites are outside the grid; their distances are NA")
            d <- matrix(NA_real_, nrow(sites), nrow(sites))
            if(sum(ok) == 1) d[ok, ok] <- 0
            if(sum(ok) > 1) d[ok, ok] <- as.matrix(gdistance::costDistance(graph, sites[ok, , drop = FALSE]))
      }else{
            d <- lc_sites(graph, sites)
      }
      if(rate){
            d <- 1 / d
      }
      d
}


#' Least-cost travel time surface
#'
#' Map the least-cost wind travel time (or flow rate, its inverse) between one or more sites and
#' every grid cell across the domain: downwind, from the sites to each cell, or upwind, from each
#' cell to the sites. This function is a wrapper around \link[gdistance]{accCost}.
#'
#' Paths are restricted to the eight neighbor directions, so cost distances are overestimated
#' for routes between neighbor bearings; see \link{pairwise_least_cost} for magnitudes. Sites
#' are snapped to the centers of their grid cells.
#'
#' @inheritParams pairwise_least_cost
#' @param direction Either `"downwind"` (the default), for travel from the sites to each cell, or
#'    `"upwind"`, for travel from each cell to the sites.
#' @param rose A [wind_rose()], or a `wind_graph` built in advance whose direction matches
#'    `direction`; see [pairwise_least_cost()].
#' @return A single-layer SpatRaster. With `rate = FALSE`, the layer is named `hours` and gives
#'   the travel time between each cell and the nearest site (in the given direction), in hours if
#'   `trans = 1` in [wind_rose()] and wind speeds are in m/s (otherwise in relative units). With
#'   `rate = TRUE`, it is named `rate` and gives the inverse. Sites outside the grid are ignored,
#'   with a warning.
#' @examples
#' rose <- windscape_example("wind_rose")
#' site <- cbind(-105, 40)
#' from_site <- least_cost(rose, site)
#' to_site <- least_cost(rose, site, direction = "upwind")
#' terra::plot(c(from_site, to_site), main = c("from the site", "to the site"))
#' @export
least_cost <- function(rose, sites, direction = "downwind", rate = FALSE, ...){
      direction <- match.arg(direction, c("downwind", "upwind"))
      graph <- lc_graph(rose, direction, ...)
      sites <- as_site_matrix(sites, "sites")
      ok <- lc_inside(graph, sites)
      if(!any(ok)) stop("all sites are outside the grid")
      if(!all(ok)) warning(sum(!ok), " site(s) outside the grid were ignored")
      d <- gdistance::accCost(graph, sites[ok, , drop = FALSE])
      if(rate){
            d <- 1 / d
      }
      d <- rast(d)
      names(d) <- if(rate) "rate" else "hours"
      d
}


#' Least-cost paths
#'
#' Traces the least-cost (fastest) wind paths between sites, as sequences of grid cell centers.
#' The result has the same structure as [wind_trails()] output, so it can be drawn with
#' [geom_wind_path()]. Paths run between `sites` and the points in `to`: downwind, from the
#' sites to the points, or upwind, from the points to the sites. Without `to`, the points are a
#' regular grid across the domain, so the paths show the network of fastest routes downwind
#' (or upwind) of the sites.
#'
#' @inheritParams least_cost
#' @param sites,to Sites, and the points to trace paths to (downwind) or from (upwind):
#'    two-column matrices (or data frames) of longitude and latitude, or `SpatVector`s of points.
#'    If `to` is `NULL` (the default), the points are about `n` points on a regular grid across
#'    the domain, and points that can't be reached are dropped silently.
#' @param direction Either `"downwind"` (the default), for paths from `sites` to `to`, or
#'    `"upwind"`, for paths from `to` to `sites`. With explicit `to` points, the two give the
#'    same paths for swapped arguments; `direction` matters most with the grid of points and
#'    with `pairs = "nearest"`.
#' @param rose A [wind_rose()], or a `wind_graph` built in advance whose direction matches
#'    `direction`; see [pairwise_least_cost()].
#' @param pairs Which pairs to trace: `"all"` (the default) traces a path between every site
#'    and every point in `to`; `"nearest"` traces one path for each site, to or from the point
#'    in `to` with the lowest travel time; `"matched"` pairs `sites` and `to` row by row, which
#'    requires them to have the same number of rows; it isn't available without `to`.
#' @param n Approximate number of grid points, when `to` is `NULL`.
#'
#' @return A data frame with one row per path vertex, ordered along each path in the direction
#'    of travel (from the site to the point for downwind paths, and from the point to the site
#'    for upwind paths):
#' * `trail`: path id.
#' * `site`, `to`: row numbers of the path's site in `sites` and point in `to`.
#' * `step`: vertex number along the path, starting at 0 at its upwind end.
#' * `hours`: cumulative travel time along the path, in hours (if `trans = 1` in
#'    [wind_rose()] and wind speeds are in m/s). Its final value on each path equals the
#'    pair's travel time from [pairwise_least_cost()] (with `snap = TRUE`).
#' * `x`, `y`: longitude and latitude of the grid cell center.
#'
#' Pairs with no path (e.g. because wind never blows from one toward the other), and sites and
#' points outside the grid, are omitted with a warning. Pairs whose site and point fall in the same grid cell are omitted, since
#' their path has no length.
#'
#' @details
#' Paths are computed with [gdistance::shortestPath()]. Because wind is directional, the
#' fastest path from A to B generally differs from the fastest path from B to A.
#'
#' Paths move between neighboring grid cells, in eight directions, so they appear as segments
#' at 45-degree angles, including occasional staircase patterns where the optimal route lies
#' between two neighbor directions. This reflects the model's structure rather than the plot.
#' Paths from one site form a tree: once two paths meet, they share the rest of their route back
#' to the site, so shared segments show the main transport corridors. Because diagonal steps
#' between cell centers can cross without passing through a common cell, two branches of the
#' tree occasionally appear to cross.
#'
#' @examples
#' rose <- windscape_example("wind_rose")
#' site <- cbind(-105, 40)
#' destinations <- cbind(c(-95, -115, -100, -110), c(45, 35, 32, 48))
#' paths <- least_cost_paths(rose, site, destinations)
#' head(paths)
#'
#' library(ggplot2)
#' ggplot(paths, aes(x, y)) +
#'   geom_wind_path(aes(color = hours)) +
#'   coord_quickmap()
#'
#' # the network of fastest routes from the site, and to it
#' down <- least_cost_paths(rose, site, n = 200)
#' up <- least_cost_paths(rose, site, n = 200, direction = "upwind")
#' ggplot(rbind(cbind(down, direction = "downwind"), cbind(up, direction = "upwind")),
#'        aes(x, y)) +
#'   geom_wind_path(aes(color = hours), arrow = NULL) +
#'   facet_wrap(~direction) +
#'   coord_quickmap()
#' @export
least_cost_paths <- function(rose, sites, to = NULL, direction = "downwind",
                             pairs = c("all", "nearest", "matched"), n = 50, ...){
      direction <- match.arg(direction, c("downwind", "upwind"))
      pairs <- match.arg(pairs)
      graph <- lc_graph(rose, direction, ...)
      sites <- as_site_matrix(sites, "sites")
      grid <- is.null(to)
      if(grid){
            if(pairs == "matched") stop("`pairs = \"matched\"` requires `to`")
            if(length(n) != 1 || !is.numeric(n) || n < 1) stop("`n` must be a positive number")
            template <- if(inherits(rose, "SpatRaster")) rose[[1]] else terra::rast(raster::raster(graph))
            to <- grid_points(template, n)
      }
      to <- as_site_matrix(to, "to")

      # site-point pairs to trace. On an upwind graph, paths run from the sites against the
      # wind, so cost[i, j] is the travel time from point j to site i. Sites and points outside
      # the grid get no paths.
      in_s <- lc_inside(graph, sites)
      in_t <- lc_inside(graph, to)
      if(!all(in_s)) warning(sum(!in_s), " site(s) outside the grid were omitted")
      if(!all(in_t)) warning(sum(!in_t), " point(s) in `to` outside the grid were omitted")
      cost <- matrix(Inf, nrow(sites), nrow(to))
      if(any(in_s) && any(in_t)){
            cost[in_s, in_t] <- matrix(gdistance::costDistance(graph, sites[in_s, , drop = FALSE],
                                                               to[in_t, , drop = FALSE]),
                                       sum(in_s), sum(in_t))
      }
      if(pairs == "all"){
            od <- expand.grid(to = seq_len(nrow(to)), site = seq_len(nrow(sites)))[, c("site", "to")]
      }else if(pairs == "nearest"){
            best <- apply(cost, 1, function(z) if(all(!is.finite(z))) NA else which.min(z))
            od <- data.frame(site = seq_len(nrow(sites)), to = best)
            od <- od[!is.na(od$to), ]
      }else{
            if(nrow(sites) != nrow(to)) stop("with `pairs = \"matched\"`, `sites` and `to` must have the same number of rows")
            od <- data.frame(site = seq_len(nrow(sites)), to = seq_len(nrow(to)))
      }
      od <- od[in_s[od$site] & in_t[od$to], , drop = FALSE]
      template <- raster::raster(graph)
      same_cell <- raster::cellFromXY(template, sites[od$site, , drop = FALSE]) ==
            raster::cellFromXY(template, to[od$to, , drop = FALSE])
      od <- od[!same_cell, , drop = FALSE]
      unreachable <- !is.finite(cost[cbind(od$site, od$to)])
      if(any(unreachable) && !grid) warning(sum(unreachable), " site-point pair(s) have no path and were omitted")
      od <- od[!unreachable, , drop = FALSE]
      if(nrow(od) == 0){
            return(data.frame(trail = integer(0), site = integer(0), to = integer(0),
                              step = integer(0), hours = numeric(0), x = numeric(0), y = numeric(0)))
      }

      upwind <- direction == "upwind"
      tm <- gdistance::transitionMatrix(graph)
      out <- list()
      for(i in unique(od$site)){
            dest <- od$to[od$site == i]
            sl <- gdistance::shortestPath(graph, sites[i, , drop = FALSE], to[dest, , drop = FALSE],
                                          output = "SpatialLines")
            for(k in seq_along(dest)){
                  xy <- sl@lines[[k]]@Lines[[1]]@coords
                  cells <- raster::cellFromXY(template, xy)
                  step_cost <- 1 / tm[cbind(cells[-length(cells)], cells[-1])]
                  if(upwind){ # traced against the wind; report in the direction of travel
                        xy <- xy[nrow(xy):1, , drop = FALSE]
                        step_cost <- rev(step_cost)
                  }
                  out[[length(out) + 1]] <- data.frame(site = i, to = dest[k],
                                                       step = seq_len(nrow(xy)) - 1,
                                                       hours = c(0, cumsum(step_cost)),
                                                       x = xy[, 1], y = xy[, 2])
            }
      }
      out <- do.call(rbind, out)
      out$trail <- as.integer(factor(paste(out$site, out$to), levels = unique(paste(out$site, out$to))))
      out <- out[, c("trail", "site", "to", "step", "hours", "x", "y")]
      rownames(out) <- NULL
      out
}


# Which points (rows of a two-column matrix) fall inside a wind graph's grid.
lc_inside <- function(graph, xy){
      !is.na(raster::cellFromXY(raster::raster(graph), xy))
}

# Wind graph for the least-cost functions: built from a wind rose, or a prebuilt wind_graph
# checked against the analysis direction.
lc_graph <- function(rose, direction, ...){
      if(inherits(rose, "wind_graph")){
            if(length(list(...)) > 0){
                  stop("arguments in `...` are passed to wind_graph(), so they can't be used ",
                       "when `rose` is already a wind_graph", call. = FALSE)
            }
            if(!identical(rose@direction, direction)){
                  stop("`rose` is ", if(rose@direction == "upwind") "an " else "a ", rose@direction,
                       " wind_graph, but the analysis direction ",
                       "is \"", direction, "\". Pass the wind_rose instead, or a wind_graph built with ",
                       "direction = \"", direction, "\".", call. = FALSE)
            }
            return(rose)
      }
      if(!inherits(rose, "wind_rose")){
            stop("`rose` must be a wind_rose (or a wind_graph built from one)", call. = FALSE)
      }
      wind_graph(rose, direction = direction, ...)
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
      # bearing is undefined for coincident points; treat points within 1 mm as coincident, since
      # coordinates rebuilt from a grid's extent can differ from the originals by round-off
      zero <- !is.na(d) & d < 1e-6
      d[zero] <- 0
      brg[zero] <- 0
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
# centers. Returns an n x n matrix. `graph` must be a downwind graph.
lc_sites <- function(graph, sites){
      tmat <- gdistance::transitionMatrix(graph)
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
      if(!all(ok)) warning("some sites are outside the grid; their distances are NA")
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
      out
}
