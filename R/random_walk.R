#' Simulate wind dispersal by random walk
#'
#' Estimates the dispersal of particles across a windscape with a Markov random walk on the
#' grid. Each iteration, the "particle mass" in each grid cell is dispersed within its 9-cell
#' neighborhood in proportion to wind conductance. Depending on the application, this mass could
#' represent probability, numbers of seeds or spores, etc. Particles can be released once
#' (`mode = "pulse"`) or continuously (`mode = "stream"`), and can be deposited along the way
#' (`half_life`).
#'
#' @param rose A `wind_rose`.
#' @param init Where particles are released: either a two-column matrix of coordinates, which
#'    places a mass of 1 in each cell containing a coordinate, or a single-layer SpatRaster of
#'    non-negative mass values with the same geometry as `rose`. In stream mode, `init` can be
#'    read either as a release rate (mass per hour) or as a one-time release; see Value. With
#'    `direction = "upwind"`, `init` instead gives the receptor locations (and weights) whose
#'    sources are traced.
#' @param mode Release mode. `"pulse"` (the default) releases `init` once and tracks the
#'    drifting, spreading cloud through time. `"stream"` releases particles continuously and
#'    returns the steady state, where release is balanced by deposition and loss across domain
#'    edges; equivalently, this is the all-time total for a single release (see Value). It is
#'    computed directly rather than by simulating forward.
#' @param direction `"downwind"` (the default) simulates where particles released at `init`
#'    go: their outbound windshed. `"upwind"` simulates where particles reaching `init` come
#'    from: its inbound windshed, or catchment. See Value and Details.
#' @param half_life Half-life of airborne mass, in hours of transport time (assuming
#'    `trans = 1` and wind speeds in m/s; see [wind_rose()]). Particles are deposited at the
#'    continuous rate `k = log(2) / half_life`, in the cell where they are airborne. Default
#'    `Inf`: no deposition.
#' @param timescale A value between 0 and 1 that scales the time step length (see
#'    [rw_max_step()]). The default of 1 uses the longest possible step, advancing the simulation
#'    farthest per iteration. Smaller values give more accurate pulse-mode dynamics or set the
#'    step to a desired duration. Stream-mode results do not depend on `timescale`.
#' @param latitude_correction Logical: correct the compression of east-west spread on
#'    longitude/latitude grids, where cells narrow toward the poles? Default `TRUE`. This adds
#'    conductance to east and west edges without changing drift, which shortens the time step
#'    (pulse mode needs about 1.3 times as many iterations at 45 degrees, 1.9 at 60). Has no effect
#'    on projected rasters. See Details.
#' @param density Logical: express results per km^2 rather than per grid cell? Default
#'    `TRUE`, matching [pairwise_random_walk()]. On a longitude/latitude grid, cell area shrinks
#'    toward the poles, so per-cell values are biased toward lower latitudes; per-km^2 values
#'    remove that bias and are comparable across grid resolutions. Downwind, each cell's values
#'    are divided by its own area. Upwind, values are per particle released from each origin,
#'    so they are instead divided by the area of the receptor (each receptor's own, if several);
#'    with a single receptor this rescales the whole map by a constant. `flux` is divided by each
#'    cell's area in either direction. Use `density = FALSE` for per-cell masses, which sum to
#'    totals (e.g. in the pulse-mode mass balance); equivalently, multiply per-km^2 results by
#'    [terra::cellSize()].
#' @param iter Pulse mode only: number of iterations (time steps) to simulate.
#' @param record Pulse mode only: integer vector of iterations to return, between 0 (the initial
#'    state) and `iter`. Default is the final iteration only.
#' @param method Stream mode only: how to compute the steady state. `"solve"` is exact, using a
#'    sparse LU solve. `"iterate"` repeats the stream update rule until the estimated error is
#'    below `tol`; it uses less memory for very large grids. `"auto"` (the default) uses `"solve"`
#'    for grids of up to 5e5 cells and `"iterate"` otherwise.
#' @param tol Stream mode with `method = "iterate"` only: convergence tolerance, as estimated
#'    L1 error relative to total airborne mass. For finite `half_life` the error estimate is a
#'    rigorous bound; for `half_life = Inf` it is extrapolated from the observed convergence rate.
#' @param max_iter Stream mode with `method = "iterate"` only: maximum number of iterations,
#'    after which a warning is given.
#' @param flux Stream mode only: logical: also return the net flux of material, the direction
#'    and rate at which it moves through each cell? Default `FALSE`. For upwind walks, requires
#'    a finite `half_life`. See Value.
#' @param source Upwind stream mode only: where material is released, used for `origin` and
#'    `flux`: either a two-column matrix of coordinates (a unit release in each cell containing a
#'    coordinate) or a single-layer SpatRaster of non-negative release per grid cell, like
#'    `init`. Default `NULL` releases uniformly per km^2 (each cell in proportion to its area).
#'    To use a map of release per km^2, multiply it by [terra::cellSize()].
#'
#' @return A named list of two `wind_walk` rasters: `airborne` and `deposition` in pulse mode,
#'    or `residence` and `deposition` in stream mode. Masses are in the units of `init`, per
#'    km^2 by default or per grid cell with `density = FALSE` (see `density`), and hours assume `trans = 1` and wind speeds in m/s. Cells that are NA in `rose` are NA in the
#'    output.
#'
#' **Pulse mode.** Each element has one layer per iteration in `record` (named `iter0`,
#' `iter1`, ...):
#' * `airborne`: the mass airborne over each cell at that iteration.
#' * `deposition`: the cumulative mass deposited in each cell up to that iteration. Zero
#'   everywhere if `half_life = Inf`.
#'
#' At every iteration, airborne mass plus deposition plus mass lost across domain edges equals
#' the released mass (summing per-cell values, i.e. with `density = FALSE`). Use [iter_length()] to convert iterations to hours.
#'
#' **Stream mode.** Each element is a single layer, and neither depends on `timescale`. The two
#' differ in units by a factor of time:
#' * `residence`: airborne mass integrated over time, in mass-hours (units of `init` times
#'   hours).
#' * `deposition`: deposited mass, in units of `init`. It equals `k * residence`, where
#'   `k = log(2) / half_life` is the deposition rate per hour, so it is zero everywhere if
#'   `half_life = Inf`.
#'
#' Stream results have two equivalent readings. They are the steady state when `init` is
#' released continuously as a rate (mass per hour), in which case `residence` is the
#' steady-state airborne mass and `deposition` the steady-state deposition rate (mass per hour).
#' They are also the all-time totals for a single release of `init`, in which case `residence`
#' is the total mass-hours spent airborne over each cell and `deposition` the total mass
#' eventually deposited; for a unit release from one cell, `deposition` is the probability
#' distribution of where a particle lands (a probability density per km^2 by default, or a
#' probability per cell with `density = FALSE`).
#'
#' **Net flux.** With `flux = TRUE`, the list also contains `flux`, a [wind_field()] whose `u`
#' and `v` layers are the eastward and northward components of the net transport of material
#' through each cell: the flow from the cell to each neighbor minus the flow back, times the
#' distance to that neighbor, halved (each flow is shared between two cells). Units are mass
#' times km per hour (if `init` is a release rate, `trans = 1`, and wind speeds are in m/s);
#' direction and relative magnitude are usually what matter. Draw it with
#' [geom_wind_arrow()], or [stat_wind_trail()] with `fixed_length = TRUE`. Unlike the gradient
#' of `residence` or `deposition`, which describes the shape of the windshed, net flux shows the
#' transport that produces it: across a plume's flanks, for example, material moves mostly
#' downwind, not sideways down the gradient. Each cell's net outflow equals its release minus
#' its deposition.
#'
#' For upwind walks, `flux` is the transport of the material that ends up at the receptor:
#' material released per `source` (by default, uniformly, one unit per km^2 per hour) and
#' eventually deposited at the receptor (weighted by `init`). It shows the routes by which the
#' receptor's material arrives. Each cell's net outflow equals its release of such material (its
#' release times its probability of eventual deposition at the receptor) minus the amount
#' deposited in the cell, which is nonzero only at the receptor.
#'
#' **Relationship between modes.** For the same `rose`, `init`, `half_life`, and
#' `latitude_correction`, stream results are the totals that a pulse walk accumulates over all
#' time, at any `timescale`. Stream `deposition` equals pulse `deposition` once the pulse has run
#' until negligible mass remains airborne, and stream `residence` equals pulse `airborne` mass
#' integrated over time (exactly, `(1 - lambda) * t` times its sum over iterations 0, 1, 2, ...;
#' see Details). This holds because the walk is linear (particles move independently, so their
#' contributions add) and time-invariant (the same `rose` applies at every time step), which is
#' also why continuous and single releases give the same numbers.
#'
#' **Upwind walks.** With `direction = "upwind"`, each cell's values describe particles released
#' from that cell and their contribution to the receptor at `init` (weighted by its values, if a
#' raster):
#' * Stream `deposition`: the probability that a particle released from the cell is deposited at
#'   the receptor. This maps where the receptor's deposited material comes from.
#' * Stream `residence`: the time, in hours, that a unit release from the cell spends airborne
#'   over the receptor.
#' * Pulse `airborne` at iteration `k`: the probability that a particle released from the cell is
#'   airborne over the receptor `k` time steps later; and pulse `deposition`: the probability
#'   it has been deposited there by then.
#' * Stream `origin` (with a finite `half_life`): the probability that a particle deposited at
#'   the receptor was released from the cell. Where `deposition` asks "if a particle were
#'   released here, would it reach the receptor?", `origin` asks "of the particles reaching the
#'   receptor, what share came from here?" By Bayes' rule, `origin` is `deposition` times
#'   release (see `source`), normalized to sum to 1 over cells, or with `density = TRUE`, to
#'   integrate to 1 per km^2. With the default uniform release per km^2, `origin` with
#'   `density = TRUE` is proportional to `deposition`. Shares are among sources inside the
#'   domain: material from beyond its edges is not represented. For a receptor spanning several
#'   cells, `origin` covers particles landing anywhere in it.
#'
#' @details
#' ## Transition probabilities
#'
#' The wind rose is converted to a simplex of nine probabilities for each grid cell, giving the
#' rates at which particles remain in the cell or move to each of its eight neighbors.
#' Probabilities of moving to a neighbor are proportional to conductance, and retention
#' probabilities are scaled so that the cell with the highest total conductance has zero
#' retention. This keeps conductance proportional across cells while maximizing the dispersal
#' occurring at each iteration.
#'
#' This is uniformization of the continuous-time Markov chain defined by the conductances,
#' which are rates (in units of 1 / hours, if `trans = 1` and wind speeds are in m/s). With time
#' step `t` and rate matrix `Q`, the transition matrix is `P = I + tQ`. Pulse iterations are
#' therefore a first-order approximation of the continuous-time process at times `t`, `2t`,
#' ...; reducing `timescale` improves their accuracy.
#'
#' ## Deposition and the two release modes
#'
#' Each time step, airborne mass disperses and a fraction `lambda = kt / (1 + kt)` of it is
#' deposited, where `k = log(2) / half_life`. This is the uniformization of transport and
#' deposition together, so stream-mode results are exact solutions of the continuous-time
#' process and independent of `timescale`. In pulse mode, mass halves after approximately
#' `half_life` hours (exactly in the limit of small `timescale`).
#'
#' The pulse update rule is `n <- (1 - lambda) * disperse(n)`, with `lambda * n` deposited
#' in place each step before dispersal. The stream update rule is
#' `n <- n0 + (1 - lambda) * disperse(n)`, with fixed point
#' `n* = (I - (1 - lambda) P')^-1 n0`; residence is `(1 - lambda) * t * n*`. The two modes are
#' linked exactly: stream-mode residence equals `(1 - lambda) * t` times the sum of the
#' pulse-mode surfaces over all time steps, starting from step 0.
#'
#' ## Upwind walks
#'
#' Write the downwind stream solution as `n* = G n0`, where entry `G[r, s]` is the airborne mass
#' at cell `r` per unit released at cell `s`. A downwind walk from a source `s` gives column `s`
#' of `G`: everywhere that source's particles go. An upwind walk to a receptor `r` gives row
#' `r`: every source's contribution to that receptor. It is computed with the adjoint process,
#' using the transpose of the downwind transition matrix (`P` in place of `P'`), not by
#' reversing the wind, so for any source `s` and receptor `r`, the upwind value at `s` equals
#' the downwind value at `r`. The same holds for each pulse iteration.
#'
#' ## Domain edges
#'
#' Domain edges are absorbing: mass that disperses off the grid, or into NA cells, is lost and
#' never returns. Values near edges are therefore biased low, because they receive no inflow
#' from beyond the edge. In stream mode with finite `half_life`, the fraction of released mass
#' lost across edges is reported as a message; the domain should be buffered by several decay
#' lengths beyond the area of interest, until this fraction is small or the bias in the area of
#' interest is acceptable. With `half_life = Inf`, edges are the only place stream-mode mass can
#' leave, so the result measures connectivity within the chosen domain and depends on its
#' extent; every valid cell must have a path to an edge or NA cell, or an error is raised.
#'
#' ## Latitude correction
#'
#' Conductances account for latitude (see [wind_rose()]), so the speed at which mass drifts is
#' correct at all latitudes for winds aligned with a neighbor direction. However, the spread of
#' mass in a nearest-neighbor walk is numerical diffusion whose magnitude scales with hop
#' length. Because longitude/latitude cells narrow east-west toward the poles, east-west hops are
#' shorter than north-south hops, and without correction east-west spread is compressed: for
#' isotropic wind, to roughly 0.86, 0.70, and 0.49 times north-south spread at 30, 45, and 60
#' degrees latitude.
#'
#' The correction scales each cell's east-west spread rate up to what the same rose would
#' produce on square cells of the same north-south size, by adding equal conductance to its east
#' and west edges. Because these two neighbors are exactly opposite, drift is unchanged. In tests
#' against the same winds on square cells, the correction brings east-west spread to within about
#' 8 percent of the square-cell value across a range of wind regimes at 45 and 60 degrees. It
#' does not correct the orientation (east-west/north-south covariance) of spread for oblique
#' winds, which remains underestimated at high latitude, nor the small drift bias for winds
#' blowing between neighbor directions (up to about 8 percent at the equator, more at high
#' latitude). The added conductance shortens the time step by an amount set by the
#' highest-latitude cells in the domain. Spread remains dependent on grid resolution with or
#' without the correction.
#'
#' @export
random_walk <- function(rose, init, mode = c("pulse", "stream"), direction = c("downwind", "upwind"),
                        half_life = Inf, timescale = 1, latitude_correction = TRUE,
                        density = TRUE, iter = 100, record = iter,
                        method = c("auto", "solve", "iterate"), tol = 1e-8, max_iter = 1e5,
                        flux = FALSE, source = NULL){

      mode <- match.arg(mode)
      direction <- match.arg(direction)
      method <- match.arg(method)

      disperse <- function(n, p){
            a <- p
            for(k in 1:9) a[,,k] <- n * a[,,k]
            da <- dim(a)
            x0 <- rep(0, da[1])
            y0 <- rep(0, da[2])
            x1 <- rep(0, da[1]-1)
            a[,,2] <- rbind(y0, cbind(a[1:(da[1]-1), 2:(da[2]), 2], x1)) # SW
            a[,,3] <- cbind(a[1:(da[1]), 2:(da[2]), 3], x0) # W
            a[,,4] <- rbind(cbind(a[2:(da[1]), 2:(da[2]), 4], x1), y0) # NW
            a[,,5] <- rbind(a[2:(da[1]), 1:(da[2]), 5], y0) # N
            a[,,6] <- rbind(cbind(x1, a[2:(da[1]), 1:(da[2]-1), 6]), y0) # NE
            a[,,7] <- cbind(x0, a[1:(da[1]), 1:(da[2]-1), 7]) # E
            a[,,8] <- rbind(y0, cbind(x1, a[1:(da[1]-1), 1:(da[2]-1), 8])) # SE
            a[,,9] <- rbind(y0, a[1:(da[1]-1), 1:(da[2]), 9]) # S
            apply(a, c(1, 2), sum)
      }

      diffuse <- function(n, p, i, rec = i, lambda = 0){
            rec <- sort(rec)
            air <- terra::as.array(rast(n, nlyrs = length(rec), vals = 0))
            dep <- air
            n <- matrix(n, nrow(n), byrow = T)
            d <- n * 0
            if(0 %in% rec) air[,,1] <- n
            p <- terra::as.array(p)
            pb <- txtProgressBar(min = 0, max = i, initial = 0, style = 3)
            for(j in 1:i){
                  d <- d + lambda * n # airborne mass deposits in place, before dispersing
                  n <- (1 - lambda) * disperse(n, p)
                  if(j %in% rec){
                        k <- match(j, rec)
                        air[,,k] <- n
                        dep[,,k] <- d
                  }
                  setTxtProgressBar(pb, j+1)
            }
            close(pb)
            list(air = air, dep = dep)
      }

      if(!(timescale > 0 && timescale <= 1)) stop("'timescale' must be greater than 0 and less than or equal to 1.")
      if(!is.null(source) && !(mode == "stream" && direction == "upwind"))
            stop("`source` is only used by upwind stream walks")
      raw_init <- NULL
      if(density){
            area <- rw_cell_area(rose)
            if(direction == "upwind") raw_init <- init
            # upwind values are per particle released from each origin, so normalize by receptor
            # area (each receptor's own), via the receptor weights; exact by linearity
            if(direction == "upwind") init <- rw_init(rose, init) / area
      }
      if(latitude_correction) rose <- rw_latitude_correction(rose)
      t <- rw_max_step(rose) * timescale
      lambda <- rw_decay(half_life, t)
      if(isTRUE(flux) && mode != "stream") stop("`flux = TRUE` is only available in stream mode")
      if(isTRUE(flux) && direction == "upwind" && lambda == 0)
            stop("`flux = TRUE` for upwind walks requires a finite `half_life`: it traces material ",
                 "deposited at the receptor, and with no deposition there is none")
      if(mode == "stream"){
            if(!is.null(source)){
                  source <- rw_init(rose, source)
                  if(any(terra::values(source) < 0, na.rm = TRUE)) stop("`source` values must be non-negative")
            }
            out <- rw_stream(rose, init, t, lambda, method, tol, max_iter, direction,
                             flux = isTRUE(flux), source = source, raw_init = raw_init)
            if(density) out <- rw_per_area(out, area, direction)
            return(out)
      }

      if(any(record < 0 | record > iter | record != round(record)))
            stop("`record` values must be integers between 0 and `iter`")
      record <- sort(unique(record))
      message("\titeration timestep: ~", signif(t, 3),
              "\n\tsimulation duration: ~", signif(t * iter, 3),
              if(lambda > 0) paste0("\n\tdecay per step (lambda): ~", signif(lambda, 3)),
              "\n(these values are in hours, IF trans == 1 and wind_field units are m/s)")
      p <- rw_prob(rose, t)

      n <- rw_init(rose, init)

      out <- if(direction == "downwind") diffuse(n, p, iter, record, lambda) else
            diffuse_upwind(n, p, iter, record, lambda)
      walk <- function(a){
            x <- rast(n, nlyrs = length(record), vals = a)
            names(x) <- paste0("iter", record)
            as_wind_walk(x, mode = mode, n_iter = iter, iter_length = t, decay = lambda,
                         direction = direction)
      }
      out <- structure(list(airborne = walk(out$air), deposition = walk(out$dep)),
                       class = c("random_walk", "list"))
      if(density) out <- rw_per_area(out, area, direction)
      out
}


setClass("wind_walk",
         contains = "SpatRaster",
         slots = c(mode = "character",
                   direction = "character",
                   n_iter = "numeric",
                   iter_length = "numeric",
                   decay = "numeric"),
         prototype = list(decay = NA_real_, direction = "downwind"))


as_wind_walk <- function(x, mode, n_iter, iter_length, decay = NA_real_, direction = "downwind"){
      if(!inherits(x, "SpatRaster")) stop("x must be a SpatRaster")
      x <- as(as(x, "SpatRaster"), "wind_walk")
      x@mode <- mode
      x@n_iter <- n_iter
      x@iter_length <- iter_length
      x@decay <- decay
      x@direction <- direction
      x
}

#' Iteration step length of a random walk
#'
#' @param x The list returned by [random_walk()], or one of its elements.
#' @return Length of one iteration (time step), in hours if `trans = 1` and wind speeds are in m/s.
#'
#' @export
iter_length <- function(x){
      if(is.list(x)) x <- x[[1]]
      x@iter_length
}


#' Maximum iteration duration for a random walk
#'
#' Calculate the maximum possible iteration length for a random walk simulation, which is the
#' residence time of the grid cell with the greatest total conductance. Walks can proceed no
#' faster than this.
#'
#' @param rose A \code{wind_rose}.
#'
#' @export
rw_max_step <- function(rose){
      1 / minmax(sum(rose))[2, ]
}


#' Self-contribution of each source cell in a stream-mode random walk
#'
#' Computes `G_cc`, the stream-mode residence time in cell `c` per unit of release from cell
#' `c` (see [random_walk()]). This is each source's contribution to its own value, which can be
#' subtracted for a leave-one-out surface, e.g. to avoid a site predicting itself.
#'
#' The approximation `G_cc ~ (1 - lambda) * t / (1 - (1 - lambda) * p_cc)` counts only paths in
#' which mass never leaves the cell (`p_cc` is the retention probability; see [random_walk()]).
#' It omits mass that leaves and returns, which is substantial in a nearest-neighbor walk,
#' especially in fast cells where `p_cc` is near zero.
#'
#' @param rose A `wind_rose`.
#' @param half_life Half-life of airborne mass, in hours. Must match the value used for the
#'    stream-mode walk; see [random_walk()].
#' @param timescale Time step scaling factor; see [random_walk()]. Affects only the approximate
#'    values; exact values are independent of `timescale`.
#' @param latitude_correction Logical. Must match the value used for the stream-mode walk; see
#'    [random_walk()].
#' @param cells Cell numbers, or a two-column matrix of coordinates, at which to compute
#'    `G_cc`. Optional for the approximation (which defaults to every cell); required if
#'    `exact = TRUE`.
#' @param exact Logical: compute exact values rather than the approximation? Exact values are a
#'    sparse solve per cell (sharing one factorization), practical for a set of occupied cells
#'    but not for every cell of a large grid. The approximation is a lower bound. Default `FALSE`.
#' @param chunk With `exact = TRUE`: number of cells to solve for at once. Larger values are
#'    faster but use more memory.
#' @return `G_cc` in hours per unit of release: a SpatRaster of approximate values for every
#'    cell if `cells` is NULL, otherwise a numeric vector with one value per cell. For a
#'    leave-one-out surface, `residence_loo = residence - init * G_cc` and
#'    `deposition_loo = log(2) / half_life * residence_loo`, using the `residence` layer from
#'    [random_walk()] and `init` values at the source cells.
#' @export
rw_self_retention <- function(rose, half_life = Inf, timescale = 1, latitude_correction = TRUE,
                              cells = NULL, exact = FALSE, chunk = 200){
      if(!(timescale > 0 && timescale <= 1)) stop("'timescale' must be greater than 0 and less than or equal to 1.")
      if(latitude_correction) rose <- rw_latitude_correction(rose)
      t <- rw_max_step(rose) * timescale
      lambda <- rw_decay(half_life, t)
      p <- rw_prob(rose, t)
      if(inherits(cells, "matrix")) cells <- terra::cellFromXY(rose, cells)

      scale <- (1 - lambda) * t # converts per-step airborne mass to residence time
      if(!exact){
            g <- scale / (1 - (1 - lambda) * p[[1]])
            names(g) <- "self_retention"
            if(is.null(cells)) return(g)
            return(terra::values(g)[cells, 1])
      }

      if(is.null(cells)) stop("`cells` must be supplied when `exact = TRUE`")
      if(anyNA(cells)) stop("some `cells` are NA or outside the extent of `rose`")
      P <- rw_matrix(p)
      if(lambda == 0) rw_check_drainage(P)
      N <- nrow(P)
      A <- Matrix::Diagonal(N) - (1 - lambda) * Matrix::t(P)
      lu <- Matrix::lu(A)
      g <- numeric(length(cells))
      for(s in split(seq_along(cells), ceiling(seq_along(cells) / chunk))){
            B <- Matrix::sparseMatrix(i = cells[s], j = seq_along(s), x = 1, dims = c(N, length(s)))
            X <- Matrix::solve(lu, B)
            g[s] <- X[cbind(cells[s], seq_along(s))]
      }
      g <- g * scale
      g[!attr(P, "valid")[cells]] <- NA
      g
}



# Internal helpers ---------------------

rw_prob <- function(x, t){
      p <- sum(x)
      p <- 1 / t - p
      p <- c(p, x)
      p / sum(p)
}

# Displacement vectors (km) from a cell center at latitude `lat` to its 8 neighbors, in rose
# layer order (SW, W, NW, N, NE, E, SE, S), using the same geometry as rose().
rw_neighbor_displacements <- function(lat, res){
      nc <- cbind(c(-1, -1, -1, 0, 1, 1, 1, 0) * res,
                  pmax(pmin(lat + c(-1, 0, 1, 1, 1, 0, -1, -1) * res, 90), -90))
      d <- geosphere::distGeo(c(0, lat), nc) / 1000
      b <- geosphere::bearingRhumb(c(0, lat), nc) * pi / 180
      cbind(dx = d * sin(b), dy = d * cos(b))
}


# Latitude correction for random walks on longitude/latitude grids. East-west hops are shorter
# than north-south hops, so a walk spreads mass less east-west than north-south. For each cell,
# scale the east-west spread rate, sum_k c_k dx_k^2, by the ratio of north-south to east-west
# cell size (what the same rose would give on square cells of the north-south size), by adding
# equal conductance to the W and E edges. W and E neighbors are exactly opposite, so drift,
# sum_k c_k (dx_k, dy_k), is unchanged. No-op for non-lonlat rasters.
rw_latitude_correction <- function(rose){
      if(!isTRUE(terra::is.lonlat(rose, perhaps = TRUE, warn = FALSE))) return(rose)
      res <- mean(terra::res(rose)) # as in wind_rose()
      nr <- terra::nrow(rose)
      nc <- terra::ncol(rose)
      lats <- terra::yFromRow(rose, seq_len(nr))
      pad <- vapply(lats, function(lat){
            D <- rw_neighbor_displacements(lat, res)
            dx <- D[6, "dx"] # east neighbor
            dy <- D[4, "dy"] # north neighbor
            # per unit of east-west variance rate: conductance to add to each of W and E
            max(dy / dx - 1, 0) / (2 * dx^2)
      }, numeric(1))
      ewvar <- vapply(lats, function(lat) rw_neighbor_displacements(lat, res)[, "dx"]^2, numeric(8))
      v <- terra::values(rose)
      row <- rep(seq_len(nr), each = nc)
      vx <- rowSums(v * t(ewvar)[row, ]) # east-west spread rate of each cell
      delta <- vx * pad[row]
      v[, c(2, 6)] <- v[, c(2, 6)] + delta
      out <- terra::rast(rose, nlyrs = 8)
      terra::values(out) <- v
      names(out) <- names(rose)
      as_wind_rose(out, trans = rose@trans, n_steps = rose@n_steps)
}


# Cell areas in km^2, as a single-layer SpatRaster. For rasters without a usable CRS (e.g.
# planar test grids), falls back to the product of the resolution, in squared map units.
rw_cell_area <- function(rose){
      r <- as(rose, "SpatRaster")[[1]]
      tryCatch(terra::cellSize(r, unit = "km"),
               error = function(e) terra::init(r, prod(terra::res(r))))
}


# Convert random walk results to per-km^2 values: downwind, every element is divided by each
# cell's area; upwind (residence and deposition already normalized by receptor area via
# `init`), only origin and flux are.
rw_per_area <- function(out, area, direction){
      cls <- class(out)
      for(nm in names(out)){
            if(direction == "upwind" && !nm %in% c("origin", "flux")) next
            x <- out[[nm]]
            v <- terra::values(x) / terra::values(area)[, 1]
            terra::values(x) <- v
            out[[nm]] <- x
      }
      class(out) <- cls
      out
}


# Net flows (mass per hour) from each cell to each of its 8 neighbors, in rose layer order
# (SW, W, NW, N, NE, E, SE, S), given continuous-time steady-state mass (residence) `res`, the
# transition probabilities `p` (from rw_prob()), and step length `t`. The flow from c to k is
# res_c * p_ck / t minus the reverse flow res_k * p_kc / t. Flows across domain edges or into
# NA cells count as outflow only. If `h` is given (the probability, from each cell, of eventual
# deposition at a receptor), flows are restricted to material destined for the receptor:
# res_c * p_ck * h_k / t minus res_k * p_kc * h_c / t, and nothing destined leaves the domain.
rw_edge_flows <- function(res, p, t, h = NULL){
      nr <- terra::nrow(p)
      nc <- terra::ncol(p)
      v <- terra::values(p)
      valid <- stats::complete.cases(v)
      v[!valid, ] <- 0
      row <- rep(seq_len(nr), each = nc)
      col <- rep(seq_len(nc), times = nr)
      dr <- c(1, 0, -1, -1, -1, 0, 1, 1)  # SW, W, NW, N, NE, E, SE, S
      dc <- c(-1, -1, -1, 0, 1, 1, 1, 0)
      opposite <- c(5, 6, 7, 8, 1, 2, 3, 4)
      J <- matrix(0, nr * nc, 8)
      for(k in 1:8){
            r2 <- row + dr[k]
            c2 <- col + dc[k]
            inside <- r2 >= 1 & r2 <= nr & c2 >= 1 & c2 <= nc
            j <- (r2 - 1) * nc + c2
            ok <- inside & valid
            ok[ok] <- valid[j[ok]]
            back <- numeric(nr * nc)
            if(is.null(h)){
                  out <- res * v[, k + 1] / t
                  back[ok] <- res[j[ok]] * v[j[ok], opposite[k] + 1] / t
            }else{
                  hk <- numeric(nr * nc)
                  hk[ok] <- h[j[ok]]
                  out <- res * v[, k + 1] * hk / t
                  back[ok] <- res[j[ok]] * v[j[ok], opposite[k] + 1] * h[ok] / t
            }
            J[, k] <- ifelse(valid, out - back, 0)
      }
      J
}

# Net flux vector at each cell, half the sum over neighbors of (net flow x displacement in km),
# as (u, v) columns: eastward and northward components, in mass x km per hour.
rw_flux_vectors <- function(J, rose){
      nr <- terra::nrow(rose)
      nc <- terra::ncol(rose)
      lats <- terra::yFromRow(rose, seq_len(nr))
      cell <- mean(terra::res(rose))
      lonlat <- isTRUE(terra::is.lonlat(rose, perhaps = TRUE, warn = FALSE))
      D <- lapply(lats, function(l){
            if(lonlat) rw_neighbor_displacements(l, cell) else
                  cbind(dx = c(-1, -1, -1, 0, 1, 1, 1, 0) * cell, dy = c(-1, 0, 1, 1, 1, 0, -1, -1) * cell)
      })
      dx <- t(vapply(D, function(d) d[, 1], numeric(8)))[rep(seq_len(nr), each = nc), ]
      dy <- t(vapply(D, function(d) d[, 2], numeric(8)))[rep(seq_len(nr), each = nc), ]
      cbind(u = 0.5 * rowSums(J * dx), v = 0.5 * rowSums(J * dy))
}


# Upwind (adjoint) pulse: n <- (1 - lambda) P n, with lambda * n deposited each step, recording
# the iterations in `rec`. Returns arrays matching diffuse() in random_walk().
diffuse_upwind <- function(n, p, i, rec = i, lambda = 0){
      rec <- sort(rec)
      nr <- terra::nrow(n)
      nc <- terra::ncol(n)
      P <- rw_matrix(p)
      v <- terra::values(n)[, 1]
      v[is.na(v) | !attr(P, "valid")] <- 0
      d <- v * 0
      air <- dep <- array(0, c(nr, nc, length(rec)))
      as_grid <- function(z) matrix(z, nr, nc, byrow = TRUE)
      if(0 %in% rec) air[, , 1] <- as_grid(v)
      for(j in seq_len(i)){
            d <- d + lambda * v
            v <- (1 - lambda) * as.vector(P %*% v)
            if(j %in% rec){
                  k <- match(j, rec)
                  air[, , k] <- as_grid(v)
                  dep[, , k] <- as_grid(d)
            }
      }
      list(air = air, dep = dep)
}


# Starting distribution: coordinates -> unit mass per cell, or a SpatRaster used as-is.
rw_init <- function(rose, init){
      if(inherits(init, "matrix")){
            cells <- terra::cellFromXY(rose, init)
            if(anyNA(cells)) stop("some `init` coordinates fall outside the extent of `rose`")
            n <- rose[[1]]
            terra::values(n) <- 0
            n[unique(cells)] <- 1
            names(n) <- "init"
            return(n)
      }
      if(!inherits(init, "SpatRaster")) stop("`init` must be a two-column matrix or a SpatRaster")
      if(terra::nlyr(init) != 1) stop("`init` must have a single layer")
      if(!terra::compareGeom(init, rose, stopOnError = FALSE)) stop("`init` must match the geometry of `rose`")
      init
}


# Per-step deposition fraction from a half-life and step length (same time units).
# Uniformizing transport (P = I + tQ) and first-order deposition at rate k together gives
# (1 - lambda) * P with lambda = kt / (1 + kt), so stream steady states are exact.
# half_life = Inf gives lambda = 0 (no deposition).
rw_decay <- function(half_life, t){
      if(!is.numeric(half_life) || length(half_life) != 1 || is.na(half_life) || half_life <= 0)
            stop("`half_life` must be a single positive number (hours of transport time), or Inf")
      k <- log(2) / half_life
      k * t / (1 + k * t)
}


# Sparse row-stochastic (substochastic at edges) transition matrix, P[from, to], in terra
# cell order. Built from the 9-layer simplex returned by rw_prob() (stay, SW, W, NW, N, NE,
# E, SE, S), using the same neighbor offsets as disperse() in random_walk(). Cells with NA
# in any layer are treated as outside the domain: they neither release nor receive mass.
rw_matrix <- function(p){
      nr <- terra::nrow(p)
      nc <- terra::ncol(p)
      N <- nr * nc
      v <- terra::values(p)
      valid <- stats::complete.cases(v)
      v[!valid, ] <- 0
      row <- rep(seq_len(nr), each = nc)
      col <- rep(seq_len(nc), times = nr)
      dr <- c(0,  1,  0, -1, -1, -1, 0, 1, 1)  # stay, SW, W, NW, N, NE, E, SE, S
      dc <- c(0, -1, -1, -1,  0,  1, 1, 1, 0)
      ii <- jj <- xx <- vector("list", 9)
      for(k in 1:9){
            r2 <- row + dr[k]
            c2 <- col + dc[k]
            from <- which(valid & v[, k] > 0 & r2 >= 1 & r2 <= nr & c2 >= 1 & c2 <= nc)
            to <- (r2[from] - 1) * nc + c2[from]
            keep <- valid[to]
            ii[[k]] <- from[keep]
            jj[[k]] <- to[keep]
            xx[[k]] <- v[from[keep], k]
      }
      P <- Matrix::sparseMatrix(i = unlist(ii), j = unlist(jj), x = unlist(xx), dims = c(N, N))
      attr(P, "valid") <- valid
      P
}


# Without decay, a steady state exists only if every valid cell can reach a cell that leaks
# mass (to a domain edge or NA cell). Otherwise I - P' is singular.
rw_check_drainage <- function(P){
      valid <- attr(P, "valid")
      reach <- valid & (1 - Matrix::rowSums(P)) > 1e-12
      repeat{
            new <- valid & (reach | as.vector(P %*% as.numeric(reach)) > 0)
            if(all(new == reach)) break
            reach <- new
      }
      stuck <- sum(valid & !reach)
      if(stuck > 0) stop(stuck, " cell(s) have no path to a domain edge or NA cell, so with ",
                         "`half_life = Inf` mass accumulates there without limit and there is no ",
                         "steady state. Use a finite `half_life`.")
      invisible(TRUE)
}


# Steady-state airborne mass under constant per-step release with per-step decay.
rw_stream <- function(rose, init, t, lambda, method = "auto", tol = 1e-8, max_iter = 1e5,
                      direction = "downwind", flux = FALSE, source = NULL, raw_init = NULL){

      n0r <- rw_init(rose, init)
      n0 <- terra::values(n0r)[, 1]
      n0[is.na(n0)] <- 0
      if(any(n0 < 0)) stop("`init` values must be non-negative")

      P <- rw_matrix(rw_prob(rose, t))
      valid <- attr(P, "valid")
      n0[!valid] <- 0
      N <- length(n0)
      if(method == "auto") method <- if(N <= 5e5) "solve" else "iterate"
      if(lambda == 0) rw_check_drainage(P)

      message("\titeration timestep: ~", signif(t, 3),
              "\n\tdecay per step (lambda): ~", signif(lambda, 3),
              "\n(timestep and half_life are in hours IF trans == 1 and wind_field units are m/s)")

      # downwind: n = n0 + (1 - lambda) P' n. Upwind (the adjoint): n = n0 + (1 - lambda) P n.
      Pt <- if(direction == "downwind") Matrix::t(P) else P
      if(method == "solve"){
            A <- Matrix::Diagonal(N) - (1 - lambda) * Pt
            n <- as.vector(Matrix::solve(A, n0))
            n <- pmax(n, 0) # remove LU round-off
            iters <- NA_real_
      }else{
            n <- n0
            d_prev <- NA
            converged <- FALSE
            for(iters in seq_len(max_iter)){
                  n_new <- n0 + (1 - lambda) * as.vector(Pt %*% n)
                  d <- sum(abs(n_new - n))
                  n <- n_new
                  if(lambda > 0){
                        # rigorous L1 contraction bound
                        err <- (1 - lambda) / lambda * d
                  }else{
                        # extrapolate from observed contraction rate
                        rho <- if(is.na(d_prev) || d_prev == 0) NA else d / d_prev
                        err <- if(is.na(rho) || rho >= 1) Inf else rho / (1 - rho) * d
                  }
                  d_prev <- d
                  if(err <= tol * sum(n)){
                        converged <- TRUE
                        break
                  }
            }
            if(!converged) warning("`max_iter` reached before convergence; estimated error ",
                                   signif(err / sum(n), 3), " (relative L1)")
      }

      # mass balance: release = deposition + edge loss (downwind only; upwind values are
      # per-origin probabilities, which don't sum to a released mass)
      if(sum(n0) > 0 && direction == "upwind" && lambda == 0){
            message("\thalf_life is Inf: results depend on domain extent")
      }else if(sum(n0) > 0 && direction == "downwind"){
            if(lambda > 0){
                  edge_loss <- (1 - lambda) * sum(n * (1 - Matrix::rowSums(P)))
                  message("\tfraction of released mass lost across domain edges: ",
                          signif(edge_loss / sum(n0), 3))
            }else{
                  message("\thalf_life is Inf: there is no deposition, and all released mass ",
                          "eventually leaves across domain edges, so results depend on domain extent")
            }
      }

      n[!valid] <- NA
      walk <- function(v, name){
            x <- terra::rast(rose, nlyrs = 1)
            terra::values(x) <- v
            names(x) <- name
            as_wind_walk(x, mode = "stream", n_iter = iters, iter_length = t, decay = lambda,
                         direction = direction)
      }
      out <- list(residence = walk(n * (1 - lambda) * t, "residence"), # continuous-time residence
                  deposition = walk(n * lambda, "deposition"))          # = k * residence

      if(direction == "upwind" && lambda > 0){
            # probability of deposition at the receptor(s) from each origin, from the receptor
            # weights as given (not normalized by receptor area), so that a multi-cell receptor
            # counts all particles landing anywhere in it
            if(is.null(raw_init)){
                  h <- n * lambda
            }else{
                  h <- terra::values(suppressMessages(
                        rw_stream(rose, raw_init, t, lambda, method, tol, max_iter, "upwind"))$deposition)[, 1]
            }
            h[!valid] <- 0
            # release per cell: `source`, or uniform per km^2 by default
            q <- if(is.null(source)) terra::values(rw_cell_area(rose))[, 1] else
                  terra::values(source)[, 1]
            q[is.na(q) | !valid] <- 0
            o <- h * q
            if(sum(o) > 0){
                  o <- o / sum(o)
            }else{
                  warning("no released material reaches the receptor, so `origin` is undefined")
                  o[] <- NaN
            }
            out$origin <- walk(o, "origin") # origins of particles deposited at the receptor
      }

      if(flux){
            res <- n * (1 - lambda) * t
            res[!valid] <- 0
            if(direction == "downwind"){
                  J <- rw_edge_flows(res, rw_prob(rose, t), t)
            }else{
                  # transport of material released per `source` (uniform per km^2 by default) that
                  # is eventually deposited at the receptor: mass from a downwind solve with that
                  # release, times the probability of deposition at the receptor
                  release <- terra::rast(rose, nlyrs = 1)
                  terra::values(release) <- q
                  m <- suppressMessages(rw_stream(rose, release, t, lambda, method, tol, max_iter,
                                                  "downwind"))
                  m <- terra::values(m$residence)[, 1]
                  m[!valid] <- 0
                  J <- rw_edge_flows(m, rw_prob(rose, t), t, h = h)
            }
            fv <- rw_flux_vectors(J, rose)
            fv[!valid, ] <- NA
            fx <- terra::rast(rose, nlyrs = 2)
            terra::values(fx) <- fv
            names(fx) <- c("u", "v")
            out$flux <- as(as(fx, "SpatRaster"), "wind_field")
      }
      structure(out, class = c("random_walk", "list"))
}
