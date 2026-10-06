#' Paths of material through a random walk
#'
#' Traces the routes by which material moves in a random walk model: for a downwind walk, the
#' paths from the source(s) to where material is deposited (or leaves the domain); for an upwind
#' walk, the paths from where material originates to the receptor(s). The result is path data,
#' like the output of [least_cost_paths()], for drawing with [geom_wind_path()]. Where
#' [least_cost_paths()] gives the single fastest route between places, these paths show the full
#' spread of routes, weighted by how much material uses each.
#'
#' **Despite the name, these paths involve no randomness.** They are not simulated random walks:
#' individual particles in a random walk zigzag unpredictably, and no path here follows one.
#' Instead, each path is a streamline of the walk's net flux of material (see `flux` in
#' [random_walk()]), computed deterministically from the stream-mode solution. A streamline is
#' the mean route of the material that travels along it, so the paths show where material goes
#' on average. Running the function twice gives identical results.
#'
#' **How paths are placed.** For a downwind walk, `n` end points are placed in proportion to
#' where material ends up: deposition in each cell, plus material leaving the domain across its
#' edges (with `half_life = Inf`, there is no deposition, so all paths end at the edges). End
#' points are placed deterministically, by dividing the total into `n` equal shares and putting
#' one end point at the middle of each share; several end points falling in one cell are spread
#' across it on a regular grid. Each path is then traced backward along the flux until it reaches
#' a source. Because each end point represents an equal share of the material, each path carries
#' an equal share too: the density of paths shows where material goes, and each path's end shows
#' where its share lands. For an upwind walk, start points are placed in proportion to the
#' `origin` of material deposited at the receptor(s), and paths are traced forward along the flux
#' to a receptor.
#'
#' **Paths to chosen points.** With `to`, paths are traced from the given points instead, like
#' [least_cost_paths()] to chosen destinations. Each path is still the mean route of the material
#' moving through its end point, but the paths no longer carry equal shares of material, so their
#' density has no meaning. Points the flux doesn't reach (e.g. upwind of the source, with a short
#' half-life) can't be traced, and are dropped with a warning.
#'
#' Since most material is typically deposited near its source, most sampled paths are short; use a
#' larger `n` to resolve long-distance routes. With `half_life = Inf`, paths depend on where the
#' domain edges are, as do all random walk results without deposition.
#'
#' @param rose A `wind_rose`.
#' @param init Source(s) for a downwind walk, or receptor(s) for an upwind walk, as for
#'    [random_walk()]: a two-column matrix of point coordinates, or a raster of weights.
#' @param to Optional points to trace paths to, as a two-column matrix (or data frame) of
#'    longitude and latitude, or a `SpatVector` of points. For a downwind walk, each path runs
#'    from a source to one of these points: the mean route by which material deposited there
#'    arrives. For an upwind walk, each path runs from one of these points to a receptor. If
#'    `NULL` (the default), `n` points are placed in proportion to where material ends up (see
#'    details).
#' @param n Number of paths, when `to` is `NULL`.
#' @param ... Further arguments passed to [random_walk()], such as `direction`, `half_life`,
#'    `latitude_correction`, `timescale`, or `wrap`. The walk is always run in stream mode with flux, so
#'    `mode`, `flux`, `density`, `iter`, and `record` can't be supplied.
#' @return A data frame with one row per point along each path, ordered from upstream to
#'    downstream: `trail` (line ID, for drawing), `path` (path ID; with `to`, the row of `to` the
#'    path was traced from), `step` (position along the path, from 0 at its upstream end), and
#'    `x` and `y` (coordinates). Each path is one trail, unless it crosses the east-west seam of a
#'    wrapped global grid (see `wrap` in [random_walk()]), where a new trail starts so that lines
#'    don't cross the map. For downwind walks, each path starts at a source and
#'    ends where its material is deposited or leaves the domain; for upwind walks, it starts
#'    where its material originates and ends at a receptor.
#' @seealso [random_walk()]; [least_cost_paths()] for fastest routes; [geom_wind_path()] to draw
#'    the paths.
#' @examples
#' \donttest{
#' library(ggplot2)
#' rose <- windscape_example("wind_rose")
#' site <- cbind(-105, 40)
#' p <- random_walk_paths(rose, site, n = 60, half_life = 48)
#'
#' # paths, with a dot where each path's share of material is deposited
#' ends <- p[!duplicated(p$path, fromLast = TRUE), ]
#' ggplot(p, aes(x, y)) +
#'   geom_wind_path(arrow = NULL, alpha = 0.6) +
#'   geom_point(data = ends, size = 0.8) +
#'   coord_quickmap()
#' }
#' @export
random_walk_paths <- function(rose, init, to = NULL, n = 50, ...){

      if(!inherits(rose, "wind_rose")) stop("`rose` must be a wind_rose")
      if(length(n) != 1 || !is.numeric(n) || n < 1 || n %% 1 != 0) stop("`n` must be a positive integer")
      dots <- list(...)
      fixed <- intersect(names(dots), c("mode", "flux", "density", "iter", "record"))
      if(length(fixed)) stop("`", fixed[1], "` can't be supplied: random_walk_paths() always runs ",
                             "a stream-mode walk with flux")
      direction <- match.arg(if(is.null(dots$direction)) "downwind" else dots$direction,
                             c("downwind", "upwind"))
      dots$direction <- direction
      wrap <- suppressWarnings(resolve_wrap(rose, dots$wrap)) # random_walk() gives any warning
      if(!is.null(to)){
            to <- as_site_matrix(to, "to")
            if(anyNA(terra::cellFromXY(rose, to))) stop("some `to` points fall outside the extent of `rose`")
      }

      w <- do.call(random_walk, c(list(rose = rose, init = init, mode = "stream", flux = TRUE,
                                       density = FALSE), dots))

      # where paths end (downwind) or start (upwind): cells weighted by material, as mass per cell
      if(direction == "downwind"){
            weight <- terra::values(w$deposition)[, 1] + rw_edge_loss(rose, w, dots)
      }else{
            weight <- terra::values(w$origin)[, 1]
      }
      weight[!is.finite(weight) | weight < 0] <- 0
      if(sum(weight) <= 0) stop("no material moves in this random walk, so there are no paths")
      pts <- if(is.null(to)) equal_share_points(w$flux, weight, n) else to
      n_paths <- nrow(pts)

      # the source or receptor cells where paths begin or end
      n0 <- terra::values(rw_init(rose, init))[, 1]
      targets <- which(!is.na(n0) & n0 > 0)
      tip <- target_coords(rose, init, targets)
      # a path arrives when it enters a target cell or one of its eight neighbors (near a source,
      # the flux spreads out from the source cell, so paths converge on it without always
      # entering it); each such cell is assigned to its nearest target
      adj <- rw_adjacent(rose, targets, wrap)
      dd <- (wrap_dx(terra::xFromCell(rose, adj[, 2]) - terra::xFromCell(rose, adj[, 1]), rose, wrap))^2 +
            (terra::yFromCell(rose, adj[, 2]) - terra::yFromCell(rose, adj[, 1]))^2
      adj <- adj[order(dd), , drop = FALSE]
      adj <- adj[!duplicated(adj[, 2]), , drop = FALSE]
      arrive <- stats::setNames(match(adj[, 1], targets), adj[, 2])

      # trace along the flux: backward to a source (downwind walks) or forward to a receptor.
      # Trace a quarter of the domain diagonal at a time, extending only paths that haven't yet
      # arrived, up to three diagonals in all. With `wrap`, paths continue across the east-west
      # seam; they are split into separate trails there at the end.
      e <- terra::ext(rose)
      span <- geosphere::distGeo(c(e$xmin, e$ymin), c(e$xmax, e$ymax)) / 1000
      cell_km <- mean(terra::res(rose)) * 111.32 * cos(mean(c(e$ymin, e$ymax)) * pi / 180)
      leg <- span / 4
      steps <- ceiling(leg / (cell_km / 4))
      trace_dir <- if(direction == "downwind") "upwind" else "downwind"
      xy <- lapply(seq_len(nrow(pts)), function(i) pts[i, , drop = FALSE])
      out <- vector("list", nrow(pts))
      active <- seq_len(nrow(pts))
      for(pass in 1:12){
            if(length(active) == 0) break
            start <- do.call(rbind, lapply(xy[active], function(p) p[nrow(p), , drop = FALSE]))
            tr <- wind_trails(w$flux, start, distance = leg, steps = steps, direction = trace_dir,
                              wrap = if(wrap) "horizontal" else "neither")
            tr <- split(tr, tr$particle) # one particle per path, even if it wraps
            still <- integer(0)
            for(k in seq_along(active)){
                  i <- active[k]
                  d <- tr[[as.character(k)]]
                  d <- d[order(abs(d$step)), c("x", "y")]
                  path <- rbind(xy[[i]], as.matrix(d[-1, , drop = FALSE]))
                  cells <- terra::cellFromXY(rose, path)
                  hit <- which(as.character(cells) %in% names(arrive))[1]
                  if(is.na(hit)){ xy[[i]] <- path; still <- c(still, i); next }
                  path <- rbind(path[seq_len(hit), , drop = FALSE],
                                tip[arrive[[as.character(cells[hit])]], , drop = FALSE])
                  if(direction == "downwind") path <- path[nrow(path):1, , drop = FALSE] # source first
                  out[[i]] <- data.frame(step = seq_len(nrow(path)) - 1, x = path[, 1], y = path[, 2])
            }
            active <- still
      }
      missed <- sum(vapply(out, is.null, logical(1)))
      if(missed > 0) warning(missed, " of ", n_paths, " paths could not be traced to a ",
                             if(direction == "downwind") "source" else "receptor",
                             " and were dropped", call. = FALSE)
      ids <- which(!vapply(out, is.null, logical(1)))
      if(length(ids) == 0) return(data.frame(trail = integer(0), path = integer(0), step = integer(0),
                                             x = numeric(0), y = numeric(0)))
      # with `to`, path IDs are rows of `to`; otherwise they number the traced paths
      if(is.null(to)) path_id <- seq_along(ids) else path_id <- ids
      out <- do.call(rbind, Map(function(d, i) cbind(path = i, d), out[ids], path_id))
      out$trail <- seam_trails(out$path, out$x, terra::xmax(rose) - terra::xmin(rose))
      out <- out[, c("trail", "path", "step", "x", "y")]
      rownames(out) <- NULL
      out
}


# Mass leaving the domain across its edges from each cell, in the units of stream-mode
# deposition with density = FALSE (residence = n (1 - lambda) t, and each step a fraction
# 1 - rowSums(P) of the remaining airborne mass leaves the domain)
rw_edge_loss <- function(rose, w, dots){
      lc <- if(is.null(dots$latitude_correction)) TRUE else dots$latitude_correction
      ts <- if(is.null(dots$timescale)) 1 else dots$timescale
      su <- suppressWarnings(rw_setup(rose, timescale = ts, latitude_correction = lc, wrap = dots$wrap))
      t <- su$t
      stay <- Matrix::rowSums(su$P)
      res <- terra::values(w$residence)[, 1]
      loss <- res * (1 - stay) / t
      loss[!is.finite(loss)] <- 0
      loss
}

# Pairs of (target cell, cell in its queen neighborhood, including itself), as a two-column
# matrix, with the first and last columns adjacent if `wrap`
rw_adjacent <- function(r, cells, wrap = FALSE){
      nr <- terra::nrow(r)
      nc <- terra::ncol(r)
      dr <- c(0, 1, 0, -1, -1, -1, 0, 1, 1)
      dc <- c(0, -1, -1, -1, 0, 1, 1, 1, 0)
      nb <- vapply(1:9, function(k) neighbor_cells(nr, nc, dr[k], dc[k], wrap)[cells],
                   numeric(length(cells)))
      adj <- cbind(rep(cells, 9), as.vector(nb))
      adj <- unique(adj[!is.na(adj[, 2]), , drop = FALSE])
      adj[order(adj[, 1], adj[, 2]), , drop = FALSE]
}

# East-west coordinate differences, taking the short way around if the grid wraps
wrap_dx <- function(dx, r, wrap = FALSE){
      if(!wrap) return(dx)
      w <- terra::xmax(r) - terra::xmin(r)
      (dx + w / 2) %% w - w / 2
}

# n points placed deterministically in proportion to cell weights: one at the middle of each of
# n equal shares of the total, with several points in a cell spread on a regular sub-grid
equal_share_points <- function(r, weight, n){
      cum <- cumsum(weight) / sum(weight)
      cells <- vapply((seq_len(n) - 0.5) / n, function(q) which(cum >= q)[1], integer(1))
      xy <- terra::xyFromCell(r, cells)
      res <- terra::res(r)
      for(cl in unique(cells)){
            k <- which(cells == cl)
            m <- length(k)
            if(m == 1) next
            g <- ceiling(sqrt(m))
            off <- (expand.grid(i = seq_len(g), j = seq_len(g))[seq_len(m), ] - 0.5) / g - 0.5
            xy[k, 1] <- xy[k, 1] + off$i * res[1]
            xy[k, 2] <- xy[k, 2] + off$j * res[2]
      }
      xy
}

# Coordinates where paths end at each source or receptor cell: the supplied point in that cell,
# for point input, or the cell center
target_coords <- function(rose, init, targets){
      xy <- terra::xyFromCell(rose, targets)
      if(is.matrix(init) || is.data.frame(init)){
            p <- as.matrix(init)[, 1:2, drop = FALSE]
            pc <- terra::cellFromXY(rose, p)
            for(i in seq_along(targets)){
                  k <- which(pc == targets[i])[1]
                  if(!is.na(k)) xy[i, ] <- p[k, ]
            }
      }
      xy
}

# About n points on a regular grid over a raster's extent, spaced roughly equally in km and
# snapped to the centers of non-NA cells
grid_points <- function(r, n){
      e <- terra::ext(r)
      lat <- mean(c(e$ymin, e$ymax))
      w <- (e$xmax - e$xmin) * cos(lat * pi / 180)
      h <- e$ymax - e$ymin
      nx <- max(1, round(sqrt(n * w / h)))
      ny <- max(1, round(n / nx))
      xs <- e$xmin + (seq_len(nx) - 0.5) * (e$xmax - e$xmin) / nx
      ys <- e$ymin + (seq_len(ny) - 0.5) * (e$ymax - e$ymin) / ny
      cells <- unique(terra::cellFromXY(r, as.matrix(expand.grid(xs, ys))))
      if(terra::hasValues(r)) cells <- cells[!is.na(terra::values(r[[1]])[cells, 1])]
      terra::xyFromCell(r, cells)
}

# Validate site coordinates: a two-column numeric matrix, data frame, or SpatVector of points
as_site_matrix <- function(s, name){
      if(inherits(s, "SpatVector")) s <- terra::crds(s)
      s <- as.matrix(s)
      if(ncol(s) != 2 || !is.numeric(s)) stop("`", name, "` must be a two-column matrix of coordinates")
      unname(s)
}
