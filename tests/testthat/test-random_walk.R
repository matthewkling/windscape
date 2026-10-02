# lambda for a 24-hour half-life: uniformized first-order deposition
lambda_24 <- function(t){ k <- log(2) / 24; k * t / (1 + k * t) }

# random_walk() ------------------------

test_that("transition matrix reproduces one step of pulse-mode dispersal", {
      r <- noisy_rose()
      set.seed(2)
      n <- terra::rast(r, nlyrs = 1, vals = runif(terra::ncell(r)))
      t <- rw_max_step(r)
      P <- rw_matrix(rw_prob(r, t))
      one_step <- quietly(random_walk(r, n, iter = 1, mode = "pulse"))
      expect_equal(as.vector(Matrix::t(P) %*% vals(n)), vals(one_step), tolerance = 1e-12)
})

test_that("returns a single-layer stream-mode wind_walk with decay stored", {
      r <- noisy_rose()
      w <- quietly(random_walk(r, cbind(-90, 35), mode = "stream", half_life = 24))
      expect_s4_class(w, "wind_walk")
      expect_equal(terra::nlyr(w), 1)
      expect_equal(names(w), "airborne")
      expect_equal(w@mode, "stream")
      expect_equal(w@decay, lambda_24(w@iter_length))
      expect_true(is.na(w@n_iter))
      expect_true(all(vals(w) >= 0))
})

test_that("direct solve and iteration agree", {
      r <- noisy_rose()
      pts <- cbind(c(-95, -90, -88), c(33, 38, 35))
      a <- quietly(random_walk(r, pts, mode = "stream", half_life = 24, method = "solve"))
      b <- quietly(random_walk(r, pts, mode = "stream", half_life = 24, method = "iterate", tol = 1e-12))
      expect_equal(vals(a), vals(b), tolerance = 1e-9)
      expect_gt(b@n_iter, 1)
})

test_that("solution is a fixed point of the update rule", {
      r <- noisy_rose()
      n0 <- point_raster(r, c(20, 100), c(1, 2))
      w <- quietly(random_walk(r, n0, mode = "stream", half_life = 12))
      lam <- w@decay
      n_star <- terra::rast(r, nlyrs = 1, vals = vals(w))
      step <- quietly(random_walk(r, n_star, iter = 1, mode = "pulse"))
      expect_equal(vals(n0) + (1 - lam) * vals(step), vals(w), tolerance = 1e-10)
})

test_that("release equals deposition plus edge loss", {
      r <- noisy_rose()
      n0 <- point_raster(r, c(50, 90))
      w <- quietly(random_walk(r, n0, mode = "stream", half_life = 48))
      P <- rw_matrix(rw_prob(r, w@iter_length))
      edge <- (1 - w@decay) * sum(vals(w) * (1 - Matrix::rowSums(P)))
      dep <- terra::global(rw_deposition(w), "sum")[[1]]
      expect_equal(dep + edge, sum(vals(n0)), tolerance = 1e-10)
      expect_gt(edge, 0)
})

test_that("edge-loss fraction is reported", {
      r <- noisy_rose()
      msgs <- capture_messages(utils::capture.output(
            random_walk(r, cbind(-90, 35), mode = "stream", half_life = 24)))
      expect_true(any(grepl("lost across domain edges", msgs)))
})

test_that("output is linear in release", {
      r <- noisy_rose()
      a <- point_raster(r, 30)
      b <- point_raster(r, 140, 3)
      ea <- quietly(random_walk(r, a, mode = "stream", half_life = 24))
      eb <- quietly(random_walk(r, b, mode = "stream", half_life = 24))
      eab <- quietly(random_walk(r, a + b, mode = "stream", half_life = 24))
      expect_equal(vals(eab), vals(ea) + vals(eb), tolerance = 1e-10)
})

test_that("coordinate init matches a unit-mass raster init", {
      r <- noisy_rose()
      pts <- cbind(c(-95, -90), c(33, 38))
      cells <- terra::cellFromXY(r, pts)
      a <- quietly(random_walk(r, pts, mode = "stream", half_life = 24))
      b <- quietly(random_walk(r, point_raster(r, cells), mode = "stream", half_life = 24))
      expect_equal(vals(a), vals(b))
})

test_that("very short half-life leaves mass at the source", {
      r <- noisy_rose()
      n0 <- point_raster(r, 90)
      w <- quietly(random_walk(r, n0, mode = "stream", half_life = 1e-9))
      expect_equal(vals(w), vals(n0), tolerance = 1e-8)
})

test_that("longer half-life carries mass farther from the source", {
      r <- noisy_rose()
      n0 <- point_raster(r, 90)
      short <- quietly(random_walk(r, n0, mode = "stream", half_life = 6))
      long <- quietly(random_walk(r, n0, mode = "stream", half_life = 48))
      share_at_source <- function(w) vals(w)[90] / sum(vals(w))
      expect_lt(share_at_source(long), share_at_source(short))
})

test_that("stream deposition is exactly invariant to timescale", {
      r <- noisy_rose()
      run <- function(ts) quietly(random_walk(r, cbind(-90, 35), mode = "stream",
                                              half_life = 24, timescale = ts))
      w1 <- run(1)
      w2 <- run(0.3)
      expect_equal(vals(rw_deposition(w1)), vals(rw_deposition(w2)), tolerance = 1e-10)
      scaled <- function(w) vals(w) * (1 - w@decay) * w@iter_length
      expect_equal(scaled(w1), scaled(w2), tolerance = 1e-10)
})

test_that("stream deposition equals the continuous-time solution", {
      r <- noisy_rose()
      t <- rw_max_step(r)
      Q <- (rw_matrix(rw_prob(r, t)) - Matrix::Diagonal(terra::ncell(r))) / t
      k <- log(2) / 24
      n0 <- point_raster(r, c(20, 100), c(1, 2))
      ref <- as.vector(k * Matrix::solve(k * Matrix::Diagonal(terra::ncell(r)) - Matrix::t(Q), vals(n0)))
      w <- quietly(random_walk(r, n0, mode = "stream", half_life = 24, timescale = 0.5))
      expect_equal(vals(rw_deposition(w)), ref, tolerance = 1e-10)
})

test_that("uniform wind displaces mass downwind; calm wind does not", {
      src <- 221 # center of a 21 x 21 grid
      calm <- uniform_rose()
      west <- uniform_rose(u = 3)
      south <- uniform_rose(v = -3)
      run <- function(r) quietly(random_walk(r, point_raster(r, src), mode = "stream", half_life = 24))
      expect_equal(unname(centroid(run(calm))), c(11, 11), tolerance = 1e-6)
      cw <- centroid(run(west))
      cs <- centroid(run(south))
      expect_gt(cw[["col"]], 11.5)                 # westerly wind -> mass moves east
      expect_equal(cw[["row"]], 11, tolerance = 1e-6)
      expect_gt(cs[["row"]], 11.5)                 # northerly wind (v < 0) -> mass moves south
})

test_that("NA cells are outside the domain and block transport", {
      r <- uniform_rose(nr = 9, nc = 9, u = 3)
      r[terra::cellFromRowCol(r, 1:9, 5)] <- NA     # full-height barrier in column 5
      w <- quietly(random_walk(r, point_raster(r, terra::cellFromRowCol(r, 5, 2)),
                               mode = "stream", half_life = 24))
      v <- matrix(vals(w), 9, 9, byrow = TRUE)
      expect_true(all(is.na(v[, 5])))
      expect_true(all(v[, 6:9] == 0))
      expect_gt(sum(v[, 1:4]), 0)
})

test_that("iteration warns when max_iter is reached", {
      r <- noisy_rose()
      expect_warning(utils::capture.output(suppressMessages(
            random_walk(r, cbind(-90, 35), mode = "stream", half_life = 1000,
                        method = "iterate", max_iter = 5))), "max_iter")
})

test_that("invalid inputs are rejected", {
      r <- noisy_rose()
      xy <- cbind(-90, 35)
      expect_error(random_walk(r, xy, mode = "stream", half_life = -1), "half_life")
      expect_error(random_walk(r, xy, mode = "stream", half_life = c(1, 2)), "half_life")
      expect_error(random_walk(r, xy, mode = "stream", half_life = NA), "half_life")
      expect_error(random_walk(r, xy, mode = "stream", half_life = NULL), "half_life")
      expect_error(random_walk(r, point_raster(r, 1, -1), mode = "stream", half_life = 24), "non-negative")
      expect_error(random_walk(r, cbind(0, 0), mode = "stream", half_life = 24), "outside")
      bad <- terra::rast(nrows = 3, ncols = 3, vals = 1)
      expect_error(random_walk(r, bad, mode = "stream", half_life = 24), "geometry")
      expect_error(random_walk(r, xy, mode = "stream", half_life = 24, timescale = 2), "timescale")
      expect_error(random_walk(r, xy, iter = 5, record = 6), "record")
      expect_error(random_walk(r, xy, iter = 5, record = 2.5), "record")
      expect_error(random_walk(r, xy, mode = "ratchet"))
      expect_error(random_walk(r, xy, mode = "emit"))
})

test_that("rw_deposition requires a wind_walk with a stored decay rate", {
      r <- noisy_rose()
      expect_error(rw_deposition(r), "wind_walk")
      w <- quietly(random_walk(r, cbind(-90, 35), mode = "stream", half_life = 24))
      w@decay <- NA_real_
      expect_error(rw_deposition(w), "decay")
})

test_that("pulse mode runs and conserves mass away from edges", {
      r <- uniform_rose(u = 1)
      w <- quietly(random_walk(r, point_raster(r, 221), iter = 3))
      expect_equal(terra::global(w, "sum")[[1]], 1, tolerance = 1e-12)
})

# pulse mode with decay ------------------------

test_that("pulse decay removes a fraction lambda of airborne mass per step", {
      r <- uniform_rose(u = 1)
      w <- quietly(random_walk(r, point_raster(r, 221), iter = 3, record = 0:3, half_life = 24))
      lam <- w@decay
      expect_equal(lam, lambda_24(w@iter_length))
      expect_equal(unname(terra::global(w, "sum")[[1]]), (1 - lam)^(0:3), tolerance = 1e-12)
})

test_that("pulse with no decay stores lambda = 0", {
      r <- noisy_rose()
      w <- quietly(random_walk(r, cbind(-90, 35), iter = 2))
      expect_equal(w@decay, 0)
})

test_that("record = 0 returns the initial state", {
      r <- noisy_rose()
      n0 <- point_raster(r, c(20, 90), c(1, 2))
      w <- quietly(random_walk(r, n0, iter = 3, record = c(3, 0), half_life = 24))
      expect_equal(names(w), c("iter0", "iter3"))
      expect_equal(vals(w[[1]]), vals(n0))
})

test_that("pulse deposition is lambda times airborne mass in each recorded layer", {
      r <- noisy_rose()
      w <- quietly(random_walk(r, cbind(-90, 35), iter = 4, record = c(2, 4), half_life = 24))
      d <- rw_deposition(w)
      expect_equal(names(d), c("deposition_iter2", "deposition_iter4"))
      expect_equal(terra::values(d), terra::values(w) * w@decay, ignore_attr = TRUE)
})

test_that("rw_deposition warns when there is no decay", {
      r <- noisy_rose()
      w <- quietly(random_walk(r, cbind(-90, 35), iter = 2))
      expect_warning(d <- rw_deposition(w), "Inf")
      expect_true(all(terra::values(d) == 0))
})


# pulse-stream equivalence ------------------------

test_that("stream steady state equals the sum of pulse snapshots", {
      r <- noisy_rose()
      n0 <- point_raster(r, c(20, 100), c(1, 2))
      s <- quietly(random_walk(r, n0, mode = "stream", half_life = 6))
      K <- ceiling(log(1e-13) / log(1 - s@decay))
      p <- quietly(random_walk(r, n0, iter = K, record = 0:K, half_life = 6))
      expect_equal(rowSums(terra::values(p)), vals(s), tolerance = 1e-10)
})


# stream mode without decay ------------------------

test_that("stream with infinite half-life is a fixed point of n0 + disperse(n)", {
      r <- noisy_rose()
      n0 <- point_raster(r, 90)
      w <- quietly(random_walk(r, n0, mode = "stream"))
      expect_equal(w@decay, 0)
      n_star <- terra::rast(r, nlyrs = 1, vals = vals(w))
      step <- quietly(random_walk(r, n_star, iter = 1))
      expect_equal(vals(n0) + vals(step), vals(w), tolerance = 1e-10)
})

test_that("stream with infinite half-life is the limit of long half-lives", {
      r <- noisy_rose()
      a <- quietly(random_walk(r, cbind(-90, 35), mode = "stream"))
      b <- quietly(random_walk(r, cbind(-90, 35), mode = "stream", half_life = 1e9))
      expect_equal(vals(a), vals(b), tolerance = 1e-6)
})

test_that("stream iteration converges without decay", {
      r <- noisy_rose()
      a <- quietly(random_walk(r, cbind(-90, 35), mode = "stream", method = "solve"))
      b <- quietly(random_walk(r, cbind(-90, 35), mode = "stream", method = "iterate", tol = 1e-10))
      expect_equal(vals(a), vals(b), tolerance = 1e-7)
})

test_that("stream with infinite half-life notes dependence on domain extent", {
      r <- noisy_rose()
      msgs <- capture_messages(random_walk(r, cbind(-90, 35), mode = "stream"))
      expect_true(any(grepl("domain extent", msgs)))
      expect_false(any(grepl("lost across domain edges", msgs)))
})

test_that("cells with no path to an edge are an error only without decay", {
      r <- noisy_rose()
      r[100] <- 0 # zero conductance: mass that arrives never leaves
      expect_error(quietly(random_walk(r, cbind(-90, 35), mode = "stream")), "no path")
      expect_error(rw_self_retention(r, cells = 100, exact = TRUE), "no path")
      w <- quietly(random_walk(r, cbind(-90, 35), mode = "stream", half_life = 24))
      expect_true(all(is.finite(vals(w))))
})


# rw_self_retention() ------------------------

test_that("exact self-retention matches a unit-source stream run", {
      r <- noisy_rose()
      cells <- c(20, 95)
      g <- rw_self_retention(r, half_life = 24, cells = cells, exact = TRUE)
      brute <- vapply(cells, function(cc){
            w <- quietly(random_walk(r, point_raster(r, cc), mode = "stream", half_life = 24))
            vals(w)[cc]
      }, numeric(1))
      expect_equal(g, brute, tolerance = 1e-10)
})

test_that("diagonal approximation is a lower bound on the exact value", {
      r <- noisy_rose()
      cells <- c(20, 95, 150)
      approx <- rw_self_retention(r, half_life = 24, cells = cells)
      exact <- rw_self_retention(r, half_life = 24, cells = cells, exact = TRUE)
      expect_true(all(approx >= 1))
      expect_true(all(approx < exact))
})

test_that("approximation raster equals 1 / (1 - (1 - lambda) * p_stay)", {
      r <- noisy_rose()
      g <- rw_self_retention(r, half_life = 24)
      expect_s4_class(g, "SpatRaster")
      t <- rw_max_step(r)
      lam <- lambda_24(t)
      p_stay <- vals(rw_prob(r, t)[[1]])
      expect_equal(vals(g), 1 / (1 - (1 - lam) * p_stay))
})

test_that("leave-one-out correction removes exactly a source's own contribution", {
      r <- noisy_rose()
      src <- c(20, 95, 150)
      e <- c(1, 2, 0.5)
      all_src <- quietly(random_walk(r, point_raster(r, src, e), mode = "stream", half_life = 24))
      g <- rw_self_retention(r, half_life = 24, cells = src, exact = TRUE)
      loo <- vals(all_src)[src] - e * g
      others <- vapply(seq_along(src), function(i){
            w <- quietly(random_walk(r, point_raster(r, src[-i], e[-i]), mode = "stream", half_life = 24))
            vals(w)[src[i]]
      }, numeric(1))
      expect_equal(loo, others, tolerance = 1e-10)
})

test_that("coordinate cells are accepted and exact mode requires cells", {
      r <- noisy_rose()
      xy <- cbind(-90, 35)
      expect_equal(rw_self_retention(r, 24, cells = xy, exact = TRUE),
                   rw_self_retention(r, 24, cells = terra::cellFromXY(r, xy), exact = TRUE))
      expect_error(rw_self_retention(r, 24, exact = TRUE), "cells")
      expect_error(rw_self_retention(r, -5), "half_life")
})

test_that("exact self-retention works without decay", {
      r <- noisy_rose()
      g <- rw_self_retention(r, cells = 95, exact = TRUE)
      w <- quietly(random_walk(r, point_raster(r, 95), mode = "stream"))
      expect_equal(g, vals(w)[95], tolerance = 1e-10)
      expect_lt(rw_self_retention(r, cells = 95), g)
})
