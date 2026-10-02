# rose() and the edge_loadings() C++ routine ------------------------

layers <- c("SW", "W", "NW", "N", "NE", "E", "SE", "S")

test_that("wind toward a neighbor loads only that neighbor, at speed / distance", {
      nd <- neighbor_distances(0)
      for(k in c(4, 6, 8, 2)){ # N, E, S, W
            toward <- c(N = 0, E = 90, S = 180, W = 270)[layers[k]]
            l <- rose(rose_input(0, toward, speed = 5))
            expect_equal(l[k], 5 * 3600 / nd[k], tolerance = 1e-10)
            expect_true(all(l[-k] == 0))
      }
})

test_that("wind toward a diagonal neighbor's bearing loads only that neighbor", {
      lat <- 45
      ne <- geosphere::bearingRhumb(c(0, lat), c(1, lat + 1))
      l <- rose(rose_input(lat, ne))
      expect_equal(l[5], 5 * 3600 / neighbor_distances(lat)[5], tolerance = 1e-8)
      expect_lt(max(l[-5]), 1e-10 * l[5])
})

test_that("wind between two neighbors is split between only those two", {
      l <- rose(rose_input(30, 20))   # between N (0) and NE (~41 deg at 30 N)
      expect_true(all(l[c("N", "NE") == layers] > 0))
      expect_true(all(l[!layers %in% c("N", "NE")] == 0))
})

test_that("loadings conserve total transformed speed", {
      set.seed(1)
      for(lat in c(-60, 0, 35, 70)){
            toward <- runif(50, 0, 360)
            speed <- runif(50, 0, 10)
            l <- rose(rose_input(lat, toward, speed))
            expect_equal(sum(l * neighbor_distances(lat) / 3600) * 50, sum(speed), tolerance = 1e-8)
      }
})

test_that("rose averages over time steps", {
      a <- rose(rose_input(10, 0, 4))
      b <- rose(rose_input(10, 90, 6))
      ab <- rose(c(10, 1, c(0, 6), c(4, 0)))   # u = (0, 6), v = (4, 0)
      expect_equal(ab, (a + b) / 2, tolerance = 1e-12)
})

test_that("trans is applied to speed", {
      x <- rose_input(20, 135, speed = 3)
      expect_equal(rose(x, trans = function(s) s^2), rose(x) * 3, tolerance = 1e-12)
      expect_true(all(rose(x, trans = function(s) s * 0) == 0))
})

test_that("reversing the wind reverses the rose at the equator", {
      set.seed(2)
      toward <- runif(20, 0, 360)
      speed <- runif(20, 1, 8)
      fwd <- rose(rose_input(0, toward, speed))
      bwd <- rose(rose_input(0, toward + 180, speed))
      expect_equal(bwd, fwd[c(5:8, 1:4)], tolerance = 1e-8)
})

test_that("calm wind gives zero conductance", {
      expect_true(all(rose(c(40, 1, rep(0, 5), rep(0, 5))) == 0))
})
