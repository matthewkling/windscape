# site handling in pairwise_least_cost() ------------------------

test_that("lc_leg_time solves the continuum least-cost problem", {
      # flow vectors: 10 km/h east, 5 km/h north, nothing else
      Vx <- matrix(c(0, 0, 0, 0, 0, 10, 0, 0), 1)
      Vy <- matrix(c(0, 0, 0, 5, 0, 0, 0, 0), 1)
      D <- rbind(c(20, 0), c(0, 10), c(20, 10), c(-1, 0), c(0, 0))
      t <- lc_leg_time(D, Vx[rep(1, 5), ], Vy[rep(1, 5), ])
      expect_equal(t, c(2, 2, 4, Inf, 0))
})

test_that("sites in the same cell get exact travel times in uniform wind", {
      r <- uv_rose(nr = 3, nc = 3, u = 5, v = 0) # 5 m/s westerly, 1-degree cells
      g <- wind_graph(r)
      ctr <- terra::xyFromCell(r, 5)
      xy <- rbind(ctr + c(-0.3, 0.1), ctr + c(0.2, 0.1)) # same cell, due east
      d <- pairwise_least_cost(g, xy)
      km <- geosphere::distGeo(xy[1, ], xy[2, ]) / 1000
      expect_equal(d[1, 2], km / 18, tolerance = 1e-3) # 5 m/s = 18 km/h
      expect_equal(d[2, 1], Inf)
      expect_equal(pairwise_least_cost(g, xy, snap = TRUE)[1, 2], 0) # snapped: same cell
})

test_that("travel times change continuously as a site crosses a cell boundary", {
      r <- noisy_rose()
      g <- wind_graph(r)
      edge <- terra::xmax(r) - 5 # a column boundary
      a <- cbind(-98.6, 34.3)
      b <- function(x) rbind(a, cbind(x, 36.7))
      d1 <- pairwise_least_cost(g, b(edge - 1e-6))
      d2 <- pairwise_least_cost(g, b(edge + 1e-6))
      expect_equal(d1, d2, tolerance = 1e-4)
      s1 <- pairwise_least_cost(g, b(edge - 1e-6), snap = TRUE)
      s2 <- pairwise_least_cost(g, b(edge + 1e-6), snap = TRUE)
      expect_false(isTRUE(all.equal(s1, s2, tolerance = 1e-4))) # snapping jumps
})

test_that("for sites at cell centers, results are close to and no greater than snapped results", {
      r <- noisy_rose()
      g <- wind_graph(r)
      xy <- terra::xyFromCell(r, c(5, 40, 77, 160))
      s <- pairwise_least_cost(g, xy)
      grid <- pairwise_least_cost(g, xy, snap = TRUE)
      off <- row(s) != col(s)
      expect_true(all(s[off] <= grid[off] * (1 + 1e-10)))
      expect_true(all(s[off] / grid[off] > 0.95))
})

test_that("upwind results are the transpose of downwind", {
      r <- noisy_rose()
      set.seed(3)
      xy <- cbind(runif(6, -99.5, -86), runif(6, 30.5, 41.5))
      xy <- rbind(xy, xy[1, ] + c(0.1, 0.05)) # a pair within one cell
      down <- pairwise_least_cost(wind_graph(r), xy)
      up <- pairwise_least_cost(wind_graph(r, direction = "upwind"), xy)
      expect_equal(up, t(down), tolerance = 1e-10)
})

test_that("results are finite and positive for distinct sites, zero on the diagonal", {
      r <- noisy_rose()
      g <- wind_graph(r)
      ctr <- terra::xyFromCell(r, 80)
      xy <- rbind(ctr + c(-0.2, -0.2), ctr + c(0.2, 0.2), ctr + c(1.3, 0.1))
      d <- pairwise_least_cost(g, xy)
      expect_equal(unname(diag(d)), c(0, 0, 0))
      expect_true(all(d[row(d) != col(d)] > 0))
      expect_true(all(is.finite(d)))
      expect_equal(pairwise_least_cost(g, xy, rate = TRUE), 1 / d)
})

test_that("sites outside the grid get NA with a warning", {
      r <- noisy_rose()
      xy <- rbind(terra::xyFromCell(r, c(5, 40)), c(0, 0))
      expect_warning(d <- pairwise_least_cost(wind_graph(r), xy), "outside")
      expect_true(all(is.na(d[3, ])) && all(is.na(d[, 3])))
      expect_true(is.finite(d[1, 2]) || d[1, 2] == Inf)
})
