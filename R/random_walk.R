#' Simulate diffusion by wind advection
#'
#' This function estimates diffusion across a windscape by Markov random walk. On each iteration,
#' the "particle mass" in each grid cell is dispersed within the local 9-cell neighborhood in
#' proportion with wind conductance. Depending on usage, this "mass" could represent probability,
#' number of individuals, etc. As the simulation proceeds, the mass diffuses across the landscape.
#'
#' The input wind rose raster is converted to a simplex of nine probabilities for each grid
#' cell, giving the rates at which particles are retained in a cell or moved to each of its
#' eight neighbors. Probabilities of moving to a neighboring cell are proportional to conductance
#' in the wind rose data set, with the probabilities of remaining in a cell scaled so that the
#' cell with the highest conductance has a zero retention probability; this allows conductance
#' to be normalized locally to a simplex while remaining proportional across cells and
#' maximizing the dispersal occurring at each iteration. This is uniformization of the
#' continuous-time Markov chain defined by the conductances (which are rates, in units of
#' 1 / time): with time step t and rate matrix Q, the transition matrix is P = I + tQ. Pulse
#' iterations are therefore a first-order approximation of the continuous-time process at
#' times t, 2t, ...; reducing `timescale` improves their accuracy.
#'
#' The `mode` argument determines how particles are released. With `mode = "pulse"` (the
#' default), the mass in `init` is released once, at the start of the simulation, and drifts
#' and spreads like a cloud. The function returns the airborne mass at the iterations listed in
#' `record`. With `mode = "stream"`, the mass in `init` is released at every time step, and the
#' function returns the steady state: the airborne mass in each cell once release and loss
#' (by deposition and across domain edges) are in balance. `iter` and `record` are ignored in
#' "stream" mode.
#'
#' In both modes, airborne mass can be lost by deposition, controlled by `half_life`. In
#' continuous time, deposition is a first-order process at rate k = ln(2) / half_life. Each time
#' step, airborne mass disperses and a fraction lambda = kt / (1 + kt) of it is deposited
#' (removed from the air), where t is the time step length (see \link{rw_max_step}). This is
#' the uniformization of transport and deposition together, so stream-mode results are exact
#' solutions of the continuous-time process, independent of `timescale`. In pulse mode, mass
#' halves after approximately `half_life` hours (exactly in the limit of small `timescale`).
#' The pulse update rule is n <- (1 - lambda) * disperse(n). The stream update rule is
#' n <- n0 + (1 - lambda) * disperse(n), and the returned surface is its fixed point,
#' n* = (I - (1 - lambda) P')^-1 n0. The two modes are linked exactly: the stream steady state
#' equals the sum of the pulse surfaces over all time steps, starting from step 0. The default
#' `half_life = Inf` means no deposition.
#'
#' Deposition per time step is `lambda` times airborne mass (see \link{rw_deposition}), which
#' treats each cell's resident airborne mass as depositing in that cell. Because of linearity,
#' stream-mode deposition per unit of per-step release is also the probability distribution of
#' where a single particle released from `init` is deposited. Stream-mode deposition is
#' invariant to `timescale`. Airborne mass scales with 1 / t, because shorter steps mean more
#' release events per unit time; `airborne * (1 - lambda) * t` is the continuous-time
#' steady-state airborne mass per unit release rate, and is invariant to `timescale`.
#'
#' In stream mode with `method = "solve"`, the steady state is computed exactly with a sparse
#' LU solve. With `method = "iterate"`, the update rule is iterated until the estimated L1 error
#' falls below `tol` relative to the total airborne mass. For finite `half_life` the estimate is
#' the rigorous bound ((1 - lambda) / lambda) * |n_k - n_(k-1)|; for `half_life = Inf` it is
#' extrapolated from the observed rate of convergence, which is driven by loss across domain
#' edges.
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
#' Grid geometry: conductances account for latitude (see \link{wind_rose}), so the speed at which
#' mass drifts is correct at all latitudes for winds aligned with a neighbor direction. Two biases
#' remain, both inherent to a nearest-neighbor walk on a longitude/latitude grid. First, the
#' spread of mass is numerical diffusion whose magnitude scales with cell dimensions, so it
#' depends on grid resolution, and because cells narrow east-west toward the poles, spread is
#' compressed east-west at high latitude: for isotropic wind, east-west spread is roughly 0.86,
#' 0.70, and 0.49 times north-south spread at 30, 45, and 60 degrees latitude. Second, drift is
#' slightly too slow for winds blowing between neighbor directions (up to about 8 percent at the
#' equator, more at high latitude, e.g. about 14 percent for a northeast wind at 60 degrees).
#'
#' @param rose A \code{wind_rose}.
#' @param init Initial conditions from which to begin diffusion, either a two-column
#'    matrix of coordinates, or a SpatRaster layer with non-negative mass values and the
#'    same spatial properties as \code{rose}. If coordinates, the simulation starts with a
#'    mass of 1 at each coordinate location. If a SpatRaster, diffusion is done
#'    directly on the raster; this could be the ouput from a prior random_walk, or any
#'    other data representing quantities to be spatially dispersed. In stream mode, this is
#'    the mass released per time step.
#' @param iter Number of simulation iterations (positive integer). Pulse mode only.
#' @param record Integer vector specifying which iterations to record, between 0 (the initial
#'    state) and \code{iter}. Default is to record only the final iteration, i.e. \code{iter}.
#'    One raster layer is returned for each value in \code{record}. Pulse mode only.
#' @param mode Release mode: "pulse" (default) or "stream". See details.
#' @param timescale A value between 0 and 1, giving the factor by which to
#'    scale the default time step length (which is calculated from the data; see \link{rw_max_step}).
#'    At each iteration, particle mass either exits cells at the rates given in the wind rose object,
#'    or remains in the cell. The default timescale value of 1 sets the timestep length to the
#'    maximum possible value, allowing the simulation to advance as far as possible in space given the
#'    number of iterations, which is computationally optimal. Reducing this value may be useful for
#'    smoothing the simulation dynamics, and/or setting the timestep to a desired duration.
#' @param half_life Half-life of airborne mass, in hours of transport time (assuming
#'    \code{trans = 1} and wind speeds in m/s, as for \link{rw_max_step}). Default \code{Inf}
#'    (no deposition).
#' @param method Stream mode only: "auto" (default), "solve", or "iterate". "auto" uses "solve"
#'    for grids with up to 5e5 cells and "iterate" otherwise.
#' @param tol Stream mode with \code{method = "iterate"} only: relative convergence tolerance.
#' @param max_iter Stream mode with \code{method = "iterate"} only: maximum iterations.
#'
#' @return A \code{wind_walk} raster object of airborne mass. In pulse mode, one layer per
#'    value of \code{record}; in stream mode, a single layer named "airborne".
#'
#' @export
random_walk <- function(rose, init, iter = 100, record = iter, mode = c("pulse", "stream"),
                        timescale = 1, half_life = Inf,
                        method = c("auto", "solve", "iterate"), tol = 1e-8, max_iter = 1e5){

      mode <- match.arg(mode)
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
            r <- terra::as.array(rast(n, nlyrs = length(rec), vals = 0))
            n <- matrix(n, nrow(n), byrow = T)
            if(0 %in% rec) r[,,1] <- n
            p <- terra::as.array(p)
            pb <- txtProgressBar(min = 0, max = i, initial = 0, style = 3)
            for(j in 1:i){
                  n <- (1 - lambda) * disperse(n, p)
                  if(j %in% rec) r[,,match(j, rec)] <- n
                  setTxtProgressBar(pb, j+1)
            }
            close(pb)
            r
      }

      if(!(timescale > 0 && timescale <= 1)) stop("'timescale' must be greater than 0 and less than or equal to 1.")
      t <- rw_max_step(rose) * timescale
      lambda <- rw_decay(half_life, t)
      if(mode == "stream") return(rw_stream(rose, init, t, lambda, method, tol, max_iter))

      if(any(record < 0 | record > iter | record != round(record)))
            stop("`record` values must be integers between 0 and `iter`")
      record <- sort(unique(record))
      message("\titeration timestep: ~", signif(t, 3),
              "\n\tsimulation duration: ~", signif(t * iter, 3),
              if(lambda > 0) paste0("\n\tdecay per step (lambda): ~", signif(lambda, 3)),
              "\n(these values are in hours, IF trans == 1 and wind_field units are m/s)")
      p <- rw_prob(rose, t)

      n <- rw_init(rose, init)

      w <- diffuse(n, p, iter, record, lambda)
      n <- rast(n, nlyrs = length(record), vals = w)
      names(n) <- paste0("iter", record)
      as_wind_walk(n, mode = mode, n_iter = iter, iter_length = t, decay = lambda)
}


setClass("wind_walk",
         contains = "SpatRaster",
         slots = c(mode = "character",
                   n_iter = "numeric",
                   iter_length = "numeric",
                   decay = "numeric"),
         prototype = list(decay = NA_real_))


as_wind_walk <- function(x, mode, n_iter, iter_length, decay = NA_real_){
      if(!inherits(x, "SpatRaster")) stop("x must be a SpatRaster")
      x <- as(as(x, "SpatRaster"), "wind_walk")
      x@mode <- mode
      x@n_iter <- n_iter
      x@iter_length <- iter_length
      x@decay <- decay
      x
}

#' Iteration step length of a random walk
#'
#' @param x A \code{wind_walk} generated by \code{random_walk()}.
#'
#' @export
iter_length <- function(x) x@iter_length


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


#' Deposition from a random walk
#'
#' Converts airborne mass from \code{random_walk()} into deposition per time step,
#' \code{lambda * airborne}, where lambda is the per-step decay fraction implied by
#' \code{half_life}. For a stream-mode walk this is steady-state deposition per step; for a
#' pulse-mode walk it is the deposition during each recorded time step.
#'
#' @param x A \code{wind_walk} generated by \code{random_walk()}.
#' @return A SpatRaster of deposition per time step, with one layer per layer of \code{x}.
#' @export
rw_deposition <- function(x){
      if(!inherits(x, "wind_walk")) stop("`x` must be a wind_walk generated by random_walk()")
      if(is.na(x@decay)) stop("`x` has no stored decay rate; regenerate it with random_walk()")
      if(x@decay == 0) warning("`x` was generated with `half_life = Inf`, so deposition is zero")
      d <- as(x, "SpatRaster") * x@decay
      names(d) <- if(x@mode == "stream") "deposition" else sub("^iter", "deposition_iter", names(x))
      d
}


#' Self-contribution of each source cell in a stream-mode random walk
#'
#' Computes G_cc, the steady-state airborne mass in cell c per unit of per-step release
#' from cell c, i.e. the diagonal of G = (I - (1 - lambda) P')^-1. Multiplying by a cell's
#' release gives its contribution to its own value, which can be subtracted for a
#' leave-one-out surface: \code{n_loo = n - n0 * G_cc}.
#'
#' The approximation \code{G_cc ~ 1 / (1 - (1 - lambda) p_cc)} counts only paths in which
#' mass never leaves the cell, and is therefore a lower bound. It omits mass that leaves and
#' returns, which is substantial in a nearest-neighbor walk, especially in fast cells where
#' the retention probability p_cc is near zero. The exact value requires one sparse solve
#' per cell (sharing a single factorization), so is practical for a set of occupied cells
#' but not for every cell of a large grid.
#'
#' @param rose A \code{wind_rose}.
#' @param half_life Half-life of airborne mass, in hours; see \link{random_walk}. Must match
#'    the value used for the stream-mode walk.
#' @param cells Optional cell numbers, or a two-column coordinate matrix. Required if
#'    \code{exact = TRUE}.
#' @param exact Logical: compute exact values rather than the diagonal approximation?
#' @param timescale Time step scaling factor; see \link{random_walk}. Must match the value
#'    used for the stream-mode walk.
#' @param chunk Number of right-hand sides per solve when \code{exact = TRUE}.
#' @return If \code{cells} is NULL, a SpatRaster of approximate G_cc. Otherwise a numeric
#'    vector with one value per cell.
#' @export
rw_self_retention <- function(rose, half_life = Inf, cells = NULL, exact = FALSE, timescale = 1, chunk = 200){
      if(!(timescale > 0 && timescale <= 1)) stop("'timescale' must be greater than 0 and less than or equal to 1.")
      t <- rw_max_step(rose) * timescale
      lambda <- rw_decay(half_life, t)
      p <- rw_prob(rose, t)
      if(inherits(cells, "matrix")) cells <- terra::cellFromXY(rose, cells)

      if(!exact){
            g <- 1 / (1 - (1 - lambda) * p[[1]])
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
rw_stream <- function(rose, init, t, lambda, method = "auto", tol = 1e-8, max_iter = 1e5){

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

      Pt <- Matrix::t(P)
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

      # mass balance: release = deposition + edge loss
      if(sum(n0) > 0){
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
      out <- terra::rast(rose, nlyrs = 1)
      terra::values(out) <- n
      names(out) <- "airborne"
      as_wind_walk(out, mode = "stream", n_iter = iters, iter_length = t, decay = lambda)
}
