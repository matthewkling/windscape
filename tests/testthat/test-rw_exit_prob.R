# exit probability of a unit source, from a random_walk() mass balance
lost_stream <- function(r, cell, ...){
      w <- quietly(random_walk(r, point_raster(r, cell), mode = "stream", density = FALSE, ...))
      1 - sum(vals(w$deposition), na.rm = TRUE)
}

test_that("exit probability equals stream-mode edge loss from a unit source", {
      r <- noisy_rose()
      h <- vals(rw_exit_prob(r, half_life = 24))
      cells <- c(1, 20, 95, 150, terra::ncell(r))
      lost <- vapply(cells, function(cc) lost_stream(r, cc, half_life = 24), numeric(1))
      expect_equal(h[cells], lost, tolerance = 1e-8)
})

test_that("with iter, exit probability equals pulse-mode edge loss at that iteration", {
      r <- noisy_rose()
      for(ts in c(1, 0.4)){
            h <- vals(rw_exit_prob(r, half_life = 24, iter = 30, timescale = ts))
            for(cc in c(20, 95)){
                  w <- quietly(random_walk(r, point_raster(r, cc), half_life = 24, iter = 30,
                                           timescale = ts, density = FALSE))
                  lost <- 1 - sum(vals(w$airborne)) - sum(vals(w$deposition))
                  expect_equal(h[cc], lost, tolerance = 1e-10)
            }
      }
})

test_that("all-time exit probability is the limit of long pulse horizons", {
      r <- noisy_rose()
      h <- vals(rw_exit_prob(r, half_life = 6))
      hk <- vals(rw_exit_prob(r, half_life = 6, iter = 5000))
      expect_equal(hk, h, tolerance = 1e-8)
      expect_true(all(vals(rw_exit_prob(r, half_life = 6, iter = 0)) == 0))
})

test_that("exit probability solves h = (1 - lambda) (l + P h)", {
      r <- noisy_rose()
      h <- vals(rw_exit_prob(r, half_life = 24))
      rc <- rw_latitude_correction(r)
      t <- rw_max_step(rc)
      lam <- rw_decay(24, t)
      p <- rw_prob(rc, t)
      P <- rw_matrix(p)
      l <- rowSums(rw_leak(p))
      expect_equal(h, (1 - lam) * (l + as.vector(P %*% h)), tolerance = 1e-10)
})

test_that("all-time exit probability does not depend on timescale", {
      r <- noisy_rose()
      expect_equal(vals(rw_exit_prob(r, 24, timescale = 0.3)), vals(rw_exit_prob(r, 24)))
})

test_that("leak accounts for every transition dropped from the matrix", {
      r <- noisy_rose()
      r[c(40, 41, 55, 100)] <- NA
      p <- rw_prob(rw_latitude_correction(r), rw_max_step(rw_latitude_correction(r)))
      L <- rw_leak(p)
      P <- rw_matrix(p)
      valid <- attr(P, "valid")
      expect_equal(rowSums(L)[valid], 1 - Matrix::rowSums(P)[valid], tolerance = 1e-12)
      expect_true(all(L[!valid, ] == 0))
      # only neighbors of NA cells leak into them
      nbr <- unique(as.vector(terra::adjacent(r, c(40, 41, 55, 100), directions = 8)))
      far <- setdiff(which(valid), nbr)
      expect_true(all(L[far, "na"] == 0))
      expect_gt(sum(L[intersect(which(valid), nbr), "na"]), 0)
})

test_that("edge and NA exits sum to all exits", {
      r <- noisy_rose()
      r[terra::cellFromRowCol(r, 4:9, 8)] <- NA
      h <- vals(rw_exit_prob(r, half_life = 24))
      he <- vals(rw_exit_prob(r, half_life = 24, exits = "edges"))
      hn <- vals(rw_exit_prob(r, half_life = 24, exits = "na"))
      expect_equal(he + hn, h, tolerance = 1e-10)
      expect_gt(max(hn, na.rm = TRUE), 0)
      expect_true(all(is.na(h[terra::cellFromRowCol(r, 4:9, 8)])))
      # without NA cells, every exit is across an edge
      r0 <- noisy_rose()
      expect_true(all(vals(rw_exit_prob(r0, 24, exits = "na")) == 0))
      expect_equal(vals(rw_exit_prob(r0, 24, exits = "edges")), vals(rw_exit_prob(r0, 24)))
})

test_that("in a uniform westerly, exit probability is highest near the downwind edge", {
      r <- uniform_rose(nr = 15, nc = 15, u = 3)
      h <- matrix(vals(rw_exit_prob(r, half_life = 24)), 15, 15, byrow = TRUE)
      expect_gt(mean(h[, 15]), mean(h[, 1]))
      expect_true(all(diff(h[8, 4:15]) > 0)) # rises toward the east edge
})

test_that("longer half-life raises exit probability", {
      r <- noisy_rose()
      expect_true(all(vals(rw_exit_prob(r, 48)) >= vals(rw_exit_prob(r, 12))))
})

test_that("without decay, exit probability is 1 except where no exit can be reached", {
      r <- noisy_rose()
      r[100] <- 0 # zero conductance: mass that arrives never leaves
      h <- vals(rw_exit_prob(r))
      expect_equal(h[100], 0)
      expect_true(all(h[-100] > 0 & h[-100] <= 1))
      h0 <- vals(rw_exit_prob(noisy_rose()))
      expect_equal(h0, rep(1, length(h0)), tolerance = 1e-8)
})

test_that("exit probability bounds the effect of enlarging the domain", {
      big <- noisy_rose(nr = 18, nc = 21, seed = 3)
      small <- as_wind_rose(terra::crop(big, terra::ext(big) - 3), trans = 1, n_steps = 50)
      h <- vals(rw_exit_prob(small, half_life = 24))
      idx <- terra::cellFromXY(big, terra::crds(small))
      receptor <- terra::xyFromCell(small, 60)

      # upwind deposition: error at each origin is at most h there
      up_small <- quietly(random_walk(small, receptor, mode = "stream", direction = "upwind",
                                      half_life = 24, density = FALSE))$deposition
      up_big <- quietly(random_walk(big, receptor, mode = "stream", direction = "upwind",
                                    half_life = 24, density = FALSE))$deposition
      err <- abs(vals(up_big)[idx] - vals(up_small))
      expect_true(all(err <= h + 1e-10))
      expect_gt(max(err), 0)

      # downwind deposition from a unit source: total error over the small domain is at most h
      for(s in c(1, 60, 130)){
            dn_small <- quietly(random_walk(small, terra::xyFromCell(small, s), mode = "stream",
                                            half_life = 24, density = FALSE))$deposition
            dn_big <- quietly(random_walk(big, terra::xyFromCell(small, s), mode = "stream",
                                          half_life = 24, density = FALSE))$deposition
            expect_lte(sum(abs(vals(dn_big)[idx] - vals(dn_small))), h[s] + 1e-10)
      }
})

test_that("invalid inputs are rejected", {
      r <- noisy_rose()
      expect_error(rw_exit_prob(r, -1), "half_life")
      expect_error(rw_exit_prob(r, 24, exits = "coast"))
      expect_error(rw_exit_prob(r, 24, iter = 2.5), "iter")
      expect_error(rw_exit_prob(r, 24, iter = c(1, 2)), "iter")
      expect_error(rw_exit_prob(r, 24, timescale = 0), "timescale")
})
