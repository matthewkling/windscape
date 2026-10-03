# lambda for a 24-hour half-life: uniformized first-order deposition
lambda_24 <- function(t){ k <- log(2) / 24; k * t / (1 + k * t) }

# per-step airborne mass n* behind a stream-mode result (residence = (1 - lambda) * t * n*)
per_step <- function(w) vals(w$residence) / ((1 - w$residence@decay) * w$residence@iter_length)

# random_walk() ------------------------

test_that("transition matrix reproduces one step of pulse-mode dispersal", {
      r <- noisy_rose()
      set.seed(2)
      n <- terra::rast(r, nlyrs = 1, vals = runif(terra::ncell(r)))
      rc <- rw_latitude_correction(r)
      t <- rw_max_step(rc)
      P <- rw_matrix(rw_prob(rc, t))
      one_step <- quietly(random_walk(r, n, iter = 1, mode = "pulse", density = FALSE))$airborne
      expect_equal(as.vector(Matrix::t(P) %*% vals(n)), vals(one_step), tolerance = 1e-12)
})

test_that("stream mode returns residence and deposition layers", {
      r <- noisy_rose()
      w <- quietly(random_walk(r, cbind(-90, 35), mode = "stream", half_life = 24))
      expect_type(w, "list")
      expect_named(w, c("residence", "deposition"))
      for(x in w) expect_s4_class(x, "wind_walk")
      expect_equal(names(w$residence), "residence")
      expect_equal(vals(w$deposition), log(2) / 24 * vals(w$residence))
      expect_equal(w$residence@mode, "stream")
      expect_equal(w$residence@decay, lambda_24(w$residence@iter_length))
      expect_equal(iter_length(w), w$residence@iter_length)
      expect_true(is.na(w$residence@n_iter))
      expect_true(all(vals(w$residence) >= 0))
})

test_that("direct solve and iteration agree", {
      r <- noisy_rose()
      pts <- cbind(c(-95, -90, -88), c(33, 38, 35))
      a <- quietly(random_walk(r, pts, mode = "stream", half_life = 24, method = "solve"))$residence
      b <- quietly(random_walk(r, pts, mode = "stream", half_life = 24, method = "iterate", tol = 1e-12))$residence
      expect_equal(vals(a), vals(b), tolerance = 1e-9)
      expect_gt(b@n_iter, 1)
})

test_that("solution is a fixed point of the update rule", {
      r <- noisy_rose()
      n0 <- point_raster(r, c(20, 100), c(1, 2))
      w <- quietly(random_walk(r, n0, mode = "stream", half_life = 12, density = FALSE))
      lam <- w$residence@decay
      n_star <- terra::rast(r, nlyrs = 1, vals = per_step(w))
      step <- quietly(random_walk(r, n_star, iter = 1, mode = "pulse", density = FALSE))$airborne
      expect_equal(vals(n0) + (1 - lam) * vals(step), per_step(w), tolerance = 1e-10)
})

test_that("release equals deposition plus edge loss", {
      r <- noisy_rose()
      n0 <- point_raster(r, c(50, 90))
      w <- quietly(random_walk(r, n0, mode = "stream", half_life = 48, density = FALSE))
      P <- rw_matrix(rw_prob(rw_latitude_correction(r), iter_length(w)))
      edge <- (1 - w$residence@decay) * sum(per_step(w) * (1 - Matrix::rowSums(P)))
      dep <- terra::global(w$deposition, "sum")[[1]]
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
      ea <- quietly(random_walk(r, a, mode = "stream", half_life = 24))$residence
      eb <- quietly(random_walk(r, b, mode = "stream", half_life = 24))$residence
      eab <- quietly(random_walk(r, a + b, mode = "stream", half_life = 24))$residence
      expect_equal(vals(eab), vals(ea) + vals(eb), tolerance = 1e-10)
})

test_that("coordinate init matches a unit-mass raster init", {
      r <- noisy_rose()
      pts <- cbind(c(-95, -90), c(33, 38))
      cells <- terra::cellFromXY(r, pts)
      a <- quietly(random_walk(r, pts, mode = "stream", half_life = 24))$residence
      b <- quietly(random_walk(r, point_raster(r, cells), mode = "stream", half_life = 24))$residence
      expect_equal(vals(a), vals(b))
})

test_that("very short half-life leaves mass at the source", {
      r <- noisy_rose()
      n0 <- point_raster(r, 90)
      w <- quietly(random_walk(r, n0, mode = "stream", half_life = 1e-9, density = FALSE))
      expect_equal(vals(w$deposition), vals(n0), tolerance = 1e-8)
})

test_that("longer half-life carries mass farther from the source", {
      r <- noisy_rose()
      n0 <- point_raster(r, 90)
      short <- quietly(random_walk(r, n0, mode = "stream", half_life = 6))$residence
      long <- quietly(random_walk(r, n0, mode = "stream", half_life = 48))$residence
      share_at_source <- function(w) vals(w)[90] / sum(vals(w))
      expect_lt(share_at_source(long), share_at_source(short))
})

test_that("stream deposition is exactly invariant to timescale", {
      r <- noisy_rose()
      run <- function(ts) quietly(random_walk(r, cbind(-90, 35), mode = "stream",
                                              half_life = 24, timescale = ts))
      w1 <- run(1)
      w2 <- run(0.3)
      expect_equal(vals(w1$deposition), vals(w2$deposition), tolerance = 1e-10)
      expect_equal(vals(w1$residence), vals(w2$residence), tolerance = 1e-10)
      expect_lt(iter_length(w2), iter_length(w1))
})

test_that("stream deposition equals the continuous-time solution", {
      r <- noisy_rose()
      rc <- rw_latitude_correction(r)
      t <- rw_max_step(rc)
      Q <- (rw_matrix(rw_prob(rc, t)) - Matrix::Diagonal(terra::ncell(r))) / t
      k <- log(2) / 24
      n0 <- point_raster(r, c(20, 100), c(1, 2))
      ref <- as.vector(k * Matrix::solve(k * Matrix::Diagonal(terra::ncell(r)) - Matrix::t(Q), vals(n0)))
      w <- quietly(random_walk(r, n0, mode = "stream", half_life = 24, timescale = 0.5, density = FALSE))
      expect_equal(vals(w$deposition), ref, tolerance = 1e-10)
})

test_that("uniform wind displaces mass downwind; calm wind does not", {
      src <- 221 # center of a 21 x 21 grid
      calm <- uniform_rose()
      west <- uniform_rose(u = 3)
      south <- uniform_rose(v = -3)
      run <- function(r) quietly(random_walk(r, point_raster(r, src), mode = "stream", half_life = 24))$residence
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
                               mode = "stream", half_life = 24))$residence
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

test_that("pulse mode runs and conserves mass away from edges", {
      r <- uniform_rose(u = 1)
      w <- quietly(random_walk(r, point_raster(r, 221), iter = 3))$airborne
      expect_equal(terra::global(w, "sum")[[1]], 1, tolerance = 1e-12)
})

# pulse mode with decay ------------------------

test_that("pulse decay removes a fraction lambda of airborne mass per step", {
      r <- uniform_rose(u = 1)
      w <- quietly(random_walk(r, point_raster(r, 221), iter = 3, record = 0:3, half_life = 24))$airborne
      lam <- w@decay
      expect_equal(lam, lambda_24(w@iter_length))
      expect_equal(unname(terra::global(w, "sum")[[1]]), (1 - lam)^(0:3), tolerance = 1e-12)
})

test_that("pulse with no decay stores lambda = 0", {
      r <- noisy_rose()
      w <- quietly(random_walk(r, cbind(-90, 35), iter = 2))$airborne
      expect_equal(w@decay, 0)
})

test_that("record = 0 returns the initial state", {
      r <- noisy_rose()
      n0 <- point_raster(r, c(20, 90), c(1, 2))
      w <- quietly(random_walk(r, n0, iter = 3, record = c(3, 0), half_life = 24, density = FALSE))$airborne
      expect_equal(names(w), c("iter0", "iter3"))
      expect_equal(vals(w[[1]]), vals(n0))
})

test_that("pulse returns airborne mass and cumulative deposition", {
      r <- noisy_rose()
      w <- quietly(random_walk(r, cbind(-90, 35), iter = 4, record = c(0, 2, 4), half_life = 24))
      expect_named(w, c("airborne", "deposition"))
      for(x in w) expect_s4_class(x, "wind_walk")
      expect_equal(names(w$deposition), c("iter0", "iter2", "iter4"))
      expect_true(all(terra::values(w$deposition)[, 1] == 0))     # nothing deposited at release
      d <- terra::values(w$deposition)
      expect_true(all(d[, 3] >= d[, 2] & d[, 2] >= d[, 1]))   # cumulative
})

test_that("pulse deposition is lambda times the sum of earlier airborne mass", {
      r <- noisy_rose()
      w <- quietly(random_walk(r, cbind(-90, 35), iter = 5, record = 0:5, half_life = 24))
      lam <- w$airborne@decay
      a <- terra::values(w$airborne)
      d <- terra::values(w$deposition)
      for(j in 1:5) expect_equal(d[, j + 1], lam * rowSums(a[, 1:j, drop = FALSE]), tolerance = 1e-12)
})

test_that("pulse airborne + deposition + edge loss equals release", {
      r <- uniform_rose(u = 1)
      w <- quietly(random_walk(r, point_raster(r, 221), iter = 4, record = 0:4, half_life = 24))
      total <- terra::global(w$airborne, "sum")[[1]] + terra::global(w$deposition, "sum")[[1]]
      expect_equal(total, rep(1, 5), tolerance = 1e-12) # source far from edges: no loss yet
})

test_that("pulse with no decay deposits nothing", {
      r <- noisy_rose()
      w <- quietly(random_walk(r, cbind(-90, 35), iter = 2))
      expect_true(all(terra::values(w$deposition) == 0))
})


# pulse-stream equivalence ------------------------

test_that("stream steady state equals the sum of pulse snapshots", {
      r <- noisy_rose()
      n0 <- point_raster(r, c(20, 100), c(1, 2))
      s <- quietly(random_walk(r, n0, mode = "stream", half_life = 6))
      lam <- s$residence@decay
      K <- ceiling(log(1e-13) / log(1 - lam))
      p <- quietly(random_walk(r, n0, iter = K, record = 0:K, half_life = 6))
      expect_equal((1 - lam) * iter_length(s) * rowSums(terra::values(p$airborne)), vals(s$residence),
                   tolerance = 1e-10)
      # pulse cumulative deposition converges to stream deposition for a one-time release
      expect_equal(vals(p$deposition[[K + 1]]), vals(s$deposition), tolerance = 1e-10)
})


# stream mode without decay ------------------------

test_that("stream with infinite half-life is a fixed point of n0 + disperse(n)", {
      r <- noisy_rose()
      n0 <- point_raster(r, 90)
      w <- quietly(random_walk(r, n0, mode = "stream", density = FALSE))
      expect_equal(w$residence@decay, 0)
      expect_true(all(vals(w$deposition) == 0))
      n_star <- terra::rast(r, nlyrs = 1, vals = per_step(w))
      step <- quietly(random_walk(r, n_star, iter = 1, density = FALSE))$airborne
      expect_equal(vals(n0) + vals(step), per_step(w), tolerance = 1e-10)
})

test_that("stream with infinite half-life is the limit of long half-lives", {
      r <- noisy_rose()
      a <- quietly(random_walk(r, cbind(-90, 35), mode = "stream"))$residence
      b <- quietly(random_walk(r, cbind(-90, 35), mode = "stream", half_life = 1e9))$residence
      expect_equal(vals(a), vals(b), tolerance = 1e-6)
})

test_that("stream iteration converges without decay", {
      r <- noisy_rose()
      a <- quietly(random_walk(r, cbind(-90, 35), mode = "stream", method = "solve"))$residence
      b <- quietly(random_walk(r, cbind(-90, 35), mode = "stream", method = "iterate", tol = 1e-10))$residence
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
      w <- quietly(random_walk(r, cbind(-90, 35), mode = "stream", half_life = 24))$residence
      expect_true(all(is.finite(vals(w))))
})


# rw_self_retention() ------------------------

test_that("exact self-retention matches a unit-source stream run", {
      r <- noisy_rose()
      cells <- c(20, 95)
      g <- rw_self_retention(r, half_life = 24, cells = cells, exact = TRUE)
      brute <- vapply(cells, function(cc){
            w <- quietly(random_walk(r, point_raster(r, cc), mode = "stream", half_life = 24, density = FALSE))$residence
            vals(w)[cc]
      }, numeric(1))
      expect_equal(g, brute, tolerance = 1e-10)
})

test_that("diagonal approximation is a lower bound on the exact value", {
      r <- noisy_rose()
      cells <- c(20, 95, 150)
      approx <- rw_self_retention(r, half_life = 24, cells = cells)
      exact <- rw_self_retention(r, half_life = 24, cells = cells, exact = TRUE)
      expect_true(all(approx > 0))
      expect_true(all(approx < exact))
})

test_that("approximation raster equals (1 - lambda) * t / (1 - (1 - lambda) * p_stay)", {
      r <- noisy_rose()
      g <- rw_self_retention(r, half_life = 24)
      expect_s4_class(g, "SpatRaster")
      rc <- rw_latitude_correction(r)
      t <- rw_max_step(rc)
      lam <- lambda_24(t)
      p_stay <- vals(rw_prob(rc, t)[[1]])
      expect_equal(vals(g), (1 - lam) * t / (1 - (1 - lam) * p_stay))
})

test_that("leave-one-out correction removes exactly a source's own contribution", {
      r <- noisy_rose()
      src <- c(20, 95, 150)
      e <- c(1, 2, 0.5)
      all_src <- quietly(random_walk(r, point_raster(r, src, e), mode = "stream", half_life = 24, density = FALSE))$residence
      g <- rw_self_retention(r, half_life = 24, cells = src, exact = TRUE)
      loo <- vals(all_src)[src] - e * g
      others <- vapply(seq_along(src), function(i){
            w <- quietly(random_walk(r, point_raster(r, src[-i], e[-i]), mode = "stream", half_life = 24, density = FALSE))$residence
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
      w <- quietly(random_walk(r, point_raster(r, 95), mode = "stream", density = FALSE))$residence
      expect_equal(g, vals(w)[95], tolerance = 1e-10)
      expect_lt(rw_self_retention(r, cells = 95), g)
})


# latitude correction ------------------------

# one-row rose at latitude `lat` from winds u, v (single cell column repeated)
lat_rose <- function(lat, u, v, nc = 3){
      build_rose(1, nc, function(x, y) list(u = u, v = v), xmin = 0, ymin = lat - 0.5)
}

# drift and spread rates (km/h, km^2/h) of each cell's jump process
jump_moments <- function(rose){
      v <- terra::values(rose)
      lat <- terra::yFromRow(rose, terra::rowFromCell(rose, seq_len(nrow(v))))
      t(vapply(seq_len(nrow(v)), function(i){
            D <- rw_neighbor_displacements(lat[i], mean(terra::res(rose)))
            c(drift_x = sum(v[i, ] * D[, 1]), drift_y = sum(v[i, ] * D[, 2]),
              vx = sum(v[i, ] * D[, 1]^2), vy = sum(v[i, ] * D[, 2]^2))
      }, numeric(4)))
}

test_that("latitude correction makes isotropic spread isotropic", {
      th <- (1:720 - 0.5) / 720 * 2 * pi
      for(lat in c(30, 45, 60)){
            r <- lat_rose(lat, 5 * sin(th), 5 * cos(th))
            before <- jump_moments(r)[1, ]
            after <- jump_moments(rw_latitude_correction(r))[1, ]
            expect_lt(before[["vx"]] / before[["vy"]], 0.9)
            expect_equal(after[["vx"]] / after[["vy"]], 1, tolerance = 0.06)
      }
})

test_that("latitude correction matches square-cell spread across wind regimes", {
      set.seed(1)
      n <- 2000
      regimes <- list(westerly = list(u = rnorm(n, 6, 3), v = rnorm(n, 0, 3)),
                      ew_axis = list(u = rnorm(n, 0, 6), v = rnorm(n, 0, 1.5)),
                      ns_axis = list(u = rnorm(n, 0, 1.5), v = rnorm(n, 0, 6)))
      for(w in regimes){
            # reference: the same winds on ~square 1-degree cells at the equator, rescaled to
            # the north-south cell size at 60 degrees
            sq <- jump_moments(lat_rose(0, w$u, w$v))[1, ]
            k <- unname(rw_neighbor_displacements(60, 1)[4, "dy"] / rw_neighbor_displacements(0, 1)[4, "dy"])
            ll <- jump_moments(rw_latitude_correction(lat_rose(60, w$u, w$v)))[1, ]
            expect_equal(ll[["vx"]] / (sq[["vx"]] * k), 1, tolerance = 0.08)
      }
})

test_that("latitude correction preserves drift and only changes E and W", {
      r <- noisy_rose()
      rc <- rw_latitude_correction(r)
      expect_s4_class(rc, "wind_rose")
      a <- jump_moments(r)
      b <- jump_moments(rc)
      expect_equal(b[, c("drift_x", "drift_y")], a[, c("drift_x", "drift_y")], tolerance = 1e-10)
      d <- terra::values(rc) - terra::values(r)
      expect_true(all(d[, c("SW", "NW", "N", "NE", "SE", "S")] == 0))
      expect_equal(d[, "E"], d[, "W"])
      expect_true(all(d[, "E"] > 0))
})

test_that("latitude correction is a no-op on planar grids", {
      r <- uniform_rose(u = 3)
      expect_identical(terra::values(rw_latitude_correction(r)), terra::values(r))
})

test_that("latitude_correction = FALSE reproduces the uncorrected walk", {
      r <- noisy_rose()
      n <- point_raster(r, 90)
      P <- rw_matrix(rw_prob(r, rw_max_step(r)))
      w <- quietly(random_walk(r, n, iter = 1, latitude_correction = FALSE, density = FALSE))$airborne
      expect_equal(as.vector(Matrix::t(P) %*% vals(n)), vals(w), tolerance = 1e-12)
      wc <- quietly(random_walk(r, n, iter = 1, density = FALSE))$airborne
      expect_lt(wc@iter_length, w@iter_length) # added conductance shortens the time step
})

test_that("latitude correction widens east-west spread in a pulse walk", {
      th <- (1:72 - 0.5) / 72 * 2 * pi
      # domain wide enough (~2200 km E-W, ~1600 km N-S) that little mass reaches the edges
      r <- build_rose(15, 41, function(x, y) list(u = 3 * sin(th), v = 3 * cos(th)),
                      xmin = -20, ymin = 52)
      src <- terra::cellFromRowCol(r, 8, 21)
      # spread (km^2) of a pulse after the same simulated time, by axis
      spread <- function(correct){
            w <- quietly(random_walk(r, point_raster(r, src), iter = 1, latitude_correction = correct))$airborne
            hrs <- 5 * rw_max_step(r) # simulate a fixed duration
            k <- ceiling(hrs / w@iter_length)
            w <- quietly(random_walk(r, point_raster(r, src), iter = k, latitude_correction = correct))$airborne
            m <- vals(w) / sum(vals(w))
            xy <- terra::xyFromCell(r, seq_along(m))
            km_x <- (xy[, 1] - terra::xFromCell(r, src)) * 111.32 * cos(xy[, 2] * pi / 180)
            km_y <- (xy[, 2] - terra::yFromCell(r, src)) * 110.57
            c(x = sum(m * km_x^2), y = sum(m * km_y^2))
      }
      raw <- spread(FALSE)
      cor <- spread(TRUE)
      expect_lt(raw[["x"]] / raw[["y"]], 0.6)
      expect_gt(cor[["x"]] / cor[["y"]], 0.85)
      expect_equal(cor[["y"]], raw[["y"]], tolerance = 0.15)
})


# upwind walks ------------------------

test_that("upwind stream is the adjoint: upwind from r at s equals downwind from s at r", {
      r <- noisy_rose()
      for(cells in list(c(20, 140), c(90, 33), c(150, 7))){
            s <- cells[1]
            rc <- cells[2]
            down <- quietly(random_walk(r, point_raster(r, s), mode = "stream", half_life = 24))
            up <- quietly(random_walk(r, point_raster(r, rc), mode = "stream", direction = "upwind",
                                      half_life = 24))
            expect_equal(vals(up$deposition)[s], vals(down$deposition)[rc], tolerance = 1e-10)
            expect_equal(vals(up$residence)[s], vals(down$residence)[rc], tolerance = 1e-10)
      }
})

test_that("upwind pulse is the adjoint at each iteration", {
      r <- noisy_rose()
      s <- 20
      rc <- 60
      down <- quietly(random_walk(r, point_raster(r, s), iter = 6, record = c(2, 6), half_life = 24))
      up <- quietly(random_walk(r, point_raster(r, rc), iter = 6, record = c(2, 6), half_life = 24,
                                direction = "upwind"))
      expect_equal(terra::values(up$airborne)[s, ], terra::values(down$airborne)[rc, ], tolerance = 1e-12)
      expect_equal(terra::values(up$deposition)[s, ], terra::values(down$deposition)[rc, ], tolerance = 1e-12)
      expect_equal(up$airborne@direction, "upwind")
      expect_equal(down$airborne@direction, "downwind")
})

test_that("upwind deposition is a probability over origins", {
      r <- noisy_rose()
      up <- quietly(random_walk(r, point_raster(r, 90), mode = "stream", direction = "upwind",
                                half_life = 24, density = FALSE))
      d <- vals(up$deposition)
      expect_true(all(d >= 0 & d <= 1))
      # equals the direct calculation for every origin: deposition at 90 from a unit release
      P <- rw_matrix(rw_prob(rw_latitude_correction(r), iter_length(up)))
      lam <- up$residence@decay
      G <- solve(diag(terra::ncell(r)) - (1 - lam) * as.matrix(Matrix::t(P)))
      expect_equal(d, lam * G[90, ], tolerance = 1e-8)
})

test_that("upwind windsheds lie upwind", {
      r <- uniform_rose(u = 3)
      down <- quietly(random_walk(r, point_raster(r, 221), mode = "stream", half_life = 24))
      up <- quietly(random_walk(r, point_raster(r, 221), mode = "stream", half_life = 24,
                                direction = "upwind"))
      expect_gt(centroid(down$residence)[["col"]], 11)
      expect_lt(centroid(up$residence)[["col"]], 11)
      expect_equal(11 - centroid(up$residence)[["col"]], centroid(down$residence)[["col"]] - 11,
                   tolerance = 1e-8) # mirror images in a uniform field
})

test_that("upwind works with iteration and infinite half-life", {
      r <- noisy_rose()
      a <- quietly(random_walk(r, point_raster(r, 90), mode = "stream", direction = "upwind",
                               method = "solve"))
      b <- quietly(random_walk(r, point_raster(r, 90), mode = "stream", direction = "upwind",
                               method = "iterate", tol = 1e-10))
      expect_equal(vals(a$residence), vals(b$residence), tolerance = 1e-7)
      expect_true(all(vals(a$deposition) == 0))
      msgs <- capture_messages(random_walk(r, point_raster(r, 90), mode = "stream",
                                           direction = "upwind", half_life = 24))
      expect_false(any(grepl("lost across domain edges", msgs)))
      expect_error(random_walk(r, point_raster(r, 90), direction = "sideways"))
})


# net flux ------------------------

test_that("net flux satisfies the mass budget: net outflow = release - deposition", {
      r <- noisy_rose()
      init <- point_raster(r, c(30, 90), c(1, 2))
      w <- quietly(random_walk(r, init, mode = "stream", half_life = 24, flux = TRUE, density = FALSE))
      res <- vals(w$residence)
      J <- rw_edge_flows(res, rw_prob(rw_latitude_correction(r), iter_length(w)), iter_length(w))
      expect_equal(rowSums(J), vals(init) - log(2) / 24 * res, tolerance = 1e-10)
      # with no deposition, net outflow is the release itself
      w0 <- quietly(random_walk(r, init, mode = "stream", flux = TRUE, density = FALSE))
      J0 <- rw_edge_flows(vals(w0$residence), rw_prob(rw_latitude_correction(r), iter_length(w0)),
                          iter_length(w0))
      expect_equal(rowSums(J0), vals(init), tolerance = 1e-10)
})

test_that("net flux is a wind_field usable by the wind geoms", {
      r <- noisy_rose()
      w <- quietly(random_walk(r, cbind(-90, 35), mode = "stream", half_life = 24, flux = TRUE))
      expect_named(w, c("residence", "deposition", "flux"))
      expect_s4_class(w$flux, "wind_field")
      expect_equal(names(w$flux), c("u", "v"))
      d <- ggplot2::fortify(w$flux)
      expect_named(d, c("x", "y", "u", "v", "speed", "bearing"))
      p <- ggplot2::ggplot(d, ggplot2::aes(x, y)) + geom_wind_arrow(fixed_length = TRUE)
      expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
      expect_named(quietly(random_walk(r, cbind(-90, 35), mode = "stream", half_life = 24)),
                   c("residence", "deposition"))
})

test_that("in a uniform westerly, flux runs downwind, also across the plume's flanks", {
      th <- (1:72 - 0.5) / 72 * 2 * pi
      r <- build_rose(21, 41, function(x, y) list(u = 3 + 2 * sin(th), v = 2 * cos(th)),
                      xmin = -110, ymin = 30)
      src <- terra::cellFromRowCol(r, 11, 8)
      w <- quietly(random_walk(r, point_raster(r, src), mode = "stream", half_life = 48, flux = TRUE))
      f <- ggplot2::fortify(w$flux)
      rc <- terra::rowColFromCell(r, seq_len(nrow(f)))
      downwind <- rc[, 2] > 12 & f$speed > 1e-3 * max(f$speed)
      expect_true(all(abs(f$bearing[downwind] - 90) < 60))      # eastward-ish everywhere downwind
      axis <- downwind & rc[, 1] == 11
      expect_true(all(abs(f$bearing[axis] - 90) < 1))           # due east on the axis
})

test_that("flux is validated by mode and direction", {
      r <- noisy_rose()
      expect_error(quietly(random_walk(r, cbind(-90, 35), flux = TRUE)), "stream mode")
      expect_error(random_walk(r, cbind(-90, 35), mode = "stream", direction = "upwind", flux = TRUE),
                   "finite")
})


test_that("upwind flux satisfies the destined-material budget", {
      r <- noisy_rose()
      init <- point_raster(r, c(60, 140), c(1, 0.5))           # weighted receptors
      w <- quietly(random_walk(r, init, mode = "stream", direction = "upwind", half_life = 24,
                               flux = TRUE, density = FALSE))
      expect_s4_class(w$flux, "wind_field")
      t <- iter_length(w)
      rc <- rw_latitude_correction(r)
      area <- vals(terra::cellSize(methods::as(r, "SpatRaster")[[1]], unit = "km"))
      m <- vals(quietly(random_walk(r, terra::rast(r, nlyrs = 1, vals = area), mode = "stream",
                                    half_life = 24, density = FALSE))$residence)
      h <- vals(w$deposition)
      J <- rw_edge_flows(m, rw_prob(rc, t), t, h = h)
      # net outflow of destined material = destined release (area x h, at 1 unit per km^2 per
      # hour) minus destined deposition (k * mass, at the receptors, weighted)
      k <- log(2) / 24
      expect_equal(rowSums(J), area * h - k * m * vals(init), tolerance = 1e-10)
})

test_that("upwind flux converges on the receptor along the wind", {
      th <- (1:72 - 0.5) / 72 * 2 * pi
      r <- build_rose(21, 41, function(x, y) list(u = 3 + 2 * sin(th), v = 2 * cos(th)),
                      xmin = -110, ymin = 30)
      rec <- terra::cellFromRowCol(r, 11, 30)
      w <- quietly(random_walk(r, point_raster(r, rec), mode = "stream", direction = "upwind",
                               half_life = 48, flux = TRUE))
      f <- ggplot2::fortify(w$flux)
      rcl <- terra::rowColFromCell(r, seq_len(nrow(f)))
      upwind_axis <- rcl[, 1] == 11 & rcl[, 2] < 26 & rcl[, 2] > 5 & f$speed > 1e-3 * max(f$speed)
      expect_true(any(upwind_axis))
      expect_true(all(abs(f$bearing[upwind_axis] - 90) < 1)) # due east, toward the receptor
})



# density ------------------------

area_of <- function(r) vals(terra::cellSize(methods::as(r, "SpatRaster")[[1]], unit = "km"))

test_that("downwind density divides each cell's values by its area", {
      r <- noisy_rose()
      a <- area_of(r)
      s <- quietly(random_walk(r, cbind(-90, 35), mode = "stream", half_life = 24, flux = TRUE, density = FALSE))
      sd <- quietly(random_walk(r, cbind(-90, 35), mode = "stream", half_life = 24, flux = TRUE,
                                density = TRUE))
      expect_equal(vals(sd$residence), vals(s$residence) / a)
      expect_equal(vals(sd$deposition), vals(s$deposition) / a)
      expect_equal(terra::values(sd$flux), terra::values(s$flux) / a, ignore_attr = TRUE)
      expect_s3_class(sd, "random_walk")
      expect_s4_class(sd$flux, "wind_field")
      p <- quietly(random_walk(r, cbind(-90, 35), iter = 3, record = c(1, 3), half_life = 24, density = FALSE))
      pd <- quietly(random_walk(r, cbind(-90, 35), iter = 3, record = c(1, 3), half_life = 24,
                                density = TRUE))
      expect_equal(terra::values(pd$airborne), terra::values(p$airborne) / a, ignore_attr = TRUE)
      expect_equal(terra::values(pd$deposition), terra::values(p$deposition) / a, ignore_attr = TRUE)
})

test_that("upwind density divides by receptor area, each receptor's own", {
      r <- noisy_rose()
      a <- area_of(r)
      rec <- c(20, 160)                                   # receptors at different latitudes
      init <- point_raster(r, rec, c(1, 2))
      u <- quietly(random_walk(r, init, mode = "stream", direction = "upwind", half_life = 24, density = FALSE))
      ud <- quietly(random_walk(r, init, mode = "stream", direction = "upwind", half_life = 24,
                                density = TRUE))
      u1 <- quietly(random_walk(r, point_raster(r, rec[1]), mode = "stream", direction = "upwind",
                                half_life = 24, density = FALSE))
      u2 <- quietly(random_walk(r, point_raster(r, rec[2]), mode = "stream", direction = "upwind",
                                half_life = 24, density = FALSE))
      expected <- vals(u1$deposition) / a[rec[1]] + 2 * vals(u2$deposition) / a[rec[2]]
      expect_equal(vals(ud$deposition), expected, tolerance = 1e-10)
      expect_false(isTRUE(all.equal(vals(ud$deposition), vals(u$deposition) / a))) # not per origin area
})

test_that("results are per km^2 by default", {
      r <- noisy_rose()
      a <- area_of(r)
      init <- point_raster(r, c(20, 160), c(1, 2))
      for(direction in c("downwind", "upwind")){
            s <- quietly(random_walk(r, init, mode = "stream", direction = direction, half_life = 24,
                                     flux = TRUE))
            sd <- quietly(random_walk(r, init, mode = "stream", direction = direction, half_life = 24,
                                      flux = TRUE, density = TRUE))
            expect_identical(vals(s$deposition), vals(sd$deposition))
            expect_identical(terra::values(s$flux), terra::values(sd$flux))
            p <- quietly(random_walk(r, init, direction = direction, iter = 3, half_life = 24))
            pd <- quietly(random_walk(r, init, direction = direction, iter = 3, half_life = 24,
                                      density = TRUE))
            expect_identical(vals(p$airborne), vals(pd$airborne))
      }
      # downwind default differs from per-cell values by cell area
      s0 <- quietly(random_walk(r, init, mode = "stream", half_life = 24, density = FALSE))
      s <- quietly(random_walk(r, init, mode = "stream", half_life = 24))
      expect_equal(vals(s$deposition), vals(s0$deposition) / a)
})

test_that("defaults match pairwise_random_walk, with and without density", {
      r <- noisy_rose()
      cells <- c(20, 75, 140)
      sites <- terra::xyFromCell(r, cells)
      m <- pairwise_random_walk(r, sites, half_life = 24)
      w <- quietly(random_walk(r, sites[2, , drop = FALSE], mode = "stream", half_life = 24))
      expect_equal(m[2, ], vals(w$deposition)[cells], tolerance = 1e-10)
      m0 <- pairwise_random_walk(r, sites, half_life = 24, density = FALSE)
      w0 <- quietly(random_walk(r, sites[2, , drop = FALSE], mode = "stream", half_life = 24,
                                density = FALSE))
      expect_equal(m0[2, ], vals(w0$deposition)[cells], tolerance = 1e-10)
})


# origin and source ------------------------

test_that("origin follows Bayes' rule: downwind deposition from each source times its release", {
      r <- noisy_rose()
      rec <- 90
      srcs <- c(20, 75, 140)
      q <- c(1, 3, 0.5)
      source <- point_raster(r, srcs, q)
      up <- quietly(random_walk(r, point_raster(r, rec), mode = "stream", direction = "upwind",
                                half_life = 24, source = source, density = FALSE))
      expect_equal(sum(vals(up$origin)), 1)
      expect_true(all(vals(up$origin)[-srcs] == 0))          # no release, no attribution
      # independently: deposition at the receptor from each source, by downwind walks
      d <- vapply(srcs, function(s){
            w <- quietly(random_walk(r, point_raster(r, s), mode = "stream", half_life = 24, density = FALSE))
            vals(w$deposition)[rec]
      }, numeric(1))
      expect_equal(vals(up$origin)[srcs], d * q / sum(d * q), tolerance = 1e-8)
})

test_that("default origin uses uniform release per km^2", {
      r <- noisy_rose()
      a <- area_of(r)
      up <- quietly(random_walk(r, point_raster(r, 90), mode = "stream", direction = "upwind",
                                half_life = 24, density = FALSE))
      h <- vals(up$deposition)
      expect_equal(vals(up$origin), h * a / sum(h * a), tolerance = 1e-10)
      upd <- quietly(random_walk(r, point_raster(r, 90), mode = "stream", direction = "upwind",
                                 half_life = 24, density = TRUE))
      expect_equal(sum(vals(upd$origin) * a), 1, tolerance = 1e-10)       # integrates to 1 per km^2
      expect_equal(vals(upd$origin) / sum(vals(upd$origin)), h / sum(h),  # proportional to reach
                   tolerance = 1e-10)
})

test_that("origin for a multi-cell receptor ignores receptor-area normalization", {
      r <- noisy_rose()
      a <- area_of(r)
      patch <- point_raster(r, c(20, 160))                     # cells at different latitudes
      up <- quietly(random_walk(r, patch, mode = "stream", direction = "upwind", half_life = 24, density = FALSE))
      upd <- quietly(random_walk(r, patch, mode = "stream", direction = "upwind", half_life = 24,
                                 density = TRUE))
      expect_equal(vals(upd$origin) * a, vals(up$origin), tolerance = 1e-10)
})

test_that("source drives upwind flux", {
      r <- noisy_rose()
      init <- point_raster(r, 90)
      q <- terra::rast(r, nlyrs = 1, vals = rep(c(0, 2), length.out = terra::ncell(r)))
      w <- quietly(random_walk(r, init, mode = "stream", direction = "upwind", half_life = 24,
                               flux = TRUE, source = q, density = FALSE))
      t <- iter_length(w)
      m <- vals(quietly(random_walk(r, q, mode = "stream", half_life = 24, density = FALSE))$residence)
      h <- vals(w$deposition)
      J <- rw_edge_flows(m, rw_prob(rw_latitude_correction(r), t), t, h = h)
      expect_equal(rowSums(J), vals(q) * h - log(2) / 24 * m * vals(init), tolerance = 1e-10)
      expect_named(w, c("residence", "deposition", "origin", "flux"))
})

test_that("origin and source scope and validation", {
      r <- noisy_rose()
      init <- point_raster(r, 90)
      expect_false("origin" %in% names(quietly(random_walk(r, init, mode = "stream", direction = "upwind"))))
      expect_false("origin" %in% names(quietly(random_walk(r, init, mode = "stream", half_life = 24))))
      expect_error(random_walk(r, init, mode = "stream", half_life = 24, source = init), "upwind stream")
      expect_error(random_walk(r, init, direction = "upwind", half_life = 24, source = init), "upwind stream")
      expect_error(random_walk(r, init, mode = "stream", direction = "upwind", half_life = 24,
                               source = point_raster(r, 5, -1)), "non-negative")
      expect_warning(utils::capture.output(suppressMessages(
            random_walk(r, init, mode = "stream", direction = "upwind", half_life = 24,
                        source = terra::rast(r, nlyrs = 1, vals = 0)))), "undefined")
      d <- ggplot2::fortify(quietly(random_walk(r, init, mode = "stream", direction = "upwind",
                                                half_life = 24)))
      expect_named(d, c("x", "y", "residence", "deposition", "origin"))
})
