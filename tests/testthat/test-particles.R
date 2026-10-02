# wind_field(), particle_flow(), generate_particles() ------------------------

field <- function(u, v, ymin = 30, nr = 10, nc = 10){
      r <- terra::rast(nrows = nr, ncols = nc, xmin = -100, xmax = -100 + nc, ymin = ymin,
                       ymax = ymin + nr, crs = "EPSG:4326", nlyrs = 2)
      terra::values(r) <- cbind(rep(u, nr * nc), rep(v, nr * nc))
      wind_field(r)
}

test_that("wind_field validates input", {
      expect_s4_class(field(1, 1), "wind_field")
      r <- terra::rast(nrows = 3, ncols = 3, nlyrs = 3, vals = 1)
      expect_error(wind_field(r), "two layers")
      p <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3e5, ymin = 0, ymax = 3e5, nlyrs = 2, vals = 1)
      expect_error(wind_field(p), "lon-lat")
})

test_that("particles move downwind and upwind with the wind", {
      f <- field(5, 0)
      p0 <- cbind(x = -95, y = 35)
      tr <- particle_flow(f, p0, n_iter = 2, scale = 0.01)
      expect_s3_class(tr, "data.frame")
      expect_true(all(c("p", "x", "y", "v", "t") %in% names(tr)))
      expect_equal(sort(unique(tr$t)), -2:2)
      expect_true(all(diff(tr$x[order(tr$t)]) > 0))       # eastward over time
      expect_equal(unique(round(tr$y, 10)), 35)            # no north-south motion
      expect_equal(unique(tr$v), 5)
})

test_that("east-west displacement accounts for latitude", {
      # equal u and v components should move a particle along a 45-degree rhumb line
      f <- field(5, 5, ymin = 55)
      tr <- particle_flow(f, cbind(-95, 60), n_iter = 1, scale = 0.01, direction = "downwind")
      p <- as.matrix(tr[order(tr$t), c("x", "y")])
      expect_lt(abs(geosphere::bearingRhumb(p[1, ], p[2, ]) - 45), 0.5)
})

test_that("ignore_speed moves particles a fixed distance", {
      f1 <- field(2, 0)
      f2 <- field(8, 0)
      p0 <- cbind(-95, 35)
      step <- function(f){
            tr <- particle_flow(f, p0, n_iter = 1, scale = 0.05, direction = "downwind", ignore_speed = TRUE)
            diff(tr$x[order(tr$t)])
      }
      expect_equal(step(f1), step(f2))
      expect_equal(step(f1), 0.05 * aspect(35), tolerance = 1e-6)
})

test_that("particles leaving the domain are dropped unless wrapped", {
      f <- field(5, 0)
      tr <- particle_flow(f, cbind(-90.5, 35), n_iter = 10, scale = 0.05, direction = "downwind")
      expect_lt(max(tr$t), 10)
      trw <- particle_flow(f, cbind(-90.5, 35), n_iter = 10, scale = 0.05, direction = "downwind",
                           wrap = "horizontal")
      expect_equal(max(trw$t), 10)
      expect_true(all(trw$x >= -100 & trw$x <= -90))
})

test_that("generate_particles samples within the domain", {
      f <- field(1, 1, ymin = 20, nr = 50, nc = 40)
      set.seed(1)
      a <- generate_particles(f, n = 500, equalarea = FALSE)
      b <- generate_particles(f, n = 500)
      g <- generate_particles(f, n = 400, sample = "grid", equalarea = FALSE)
      for(p in list(a, b, g)){
            expect_true(all(p[, 1] >= -100 & p[, 1] <= -60 & p[, 2] >= 20 & p[, 2] <= 70))
      }
      expect_equal(nrow(a), 500)
      expect_equal(nrow(b), 500)
      expect_equal(nrow(g), 400, tolerance = 0.1)
      # equal-area sampling puts fewer points at high latitude
      expect_lt(mean(b[, 2]), mean(a[, 2]))
})

test_that("equal-area grid sampling works", {
      skip_if_not_installed("sf")
      f <- field(1, 1, ymin = 20, nr = 50, nc = 40)
      g <- generate_particles(f, n = 400, sample = "grid")
      expect_equal(ncol(g), 2)
      expect_gt(nrow(g), 100)
})

test_that("trails can be returned as sf linestrings", {
      skip_if_not_installed("sf")
      f <- field(5, 1)
      s <- particle_flow(f, cbind(c(-95, -94), c(35, 36)), n_iter = 3, scale = 0.01, sf = TRUE)
      expect_s3_class(s, "sf")
      expect_equal(nrow(s), 2)
      expect_true(all(sf::st_geometry_type(s) == "LINESTRING"))
})
