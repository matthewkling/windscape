# wind_field(), wind_trails(), generate_particles() ------------------------

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

test_that("trails move downwind and upwind with the wind", {
      f <- field(5, 0)
      tr <- wind_trails(f, cbind(x = -95, y = 35), hours = 2, steps = 4)
      expect_s3_class(tr, "data.frame")
      expect_named(tr, c("trail", "particle", "step", "hours", "x", "y", "speed"))
      expect_equal(tr$step, -2:2)
      expect_equal(tr$hours, -2:2 * 0.5)
      expect_true(all(diff(tr$x) > 0))                     # eastward, ordered upwind to downwind
      expect_equal(unique(round(tr$y, 10)), 35)            # no north-south motion
      expect_equal(unique(tr$speed), 5)
})

test_that("hours gives transport distance at the wind speed", {
      f <- field(5, 0)
      tr <- wind_trails(f, cbind(-95, 35), hours = 3, steps = 3, direction = "downwind")
      km <- diff(range(tr$x)) * 111.32 * cos(35 * pi / 180)
      expect_equal(km, 5 * 3.6 * 3, tolerance = 0.01) # 5 m/s for 3 h = 54 km
})

test_that("distance gives fixed-length trails regardless of speed", {
      step <- function(f){
            tr <- wind_trails(f, cbind(-95, 35), distance = 50, steps = 1, direction = "downwind")
            diff(tr$x) * 111.32 * cos(35 * pi / 180)
      }
      expect_equal(step(field(2, 0)), step(field(8, 0)))
      expect_equal(step(field(2, 0)), 50, tolerance = 0.01)
      expect_named(wind_trails(field(2, 0), cbind(-95, 35), distance = 50, steps = 2),
                   c("trail", "particle", "step", "km", "x", "y", "speed"))
})

test_that("east-west displacement accounts for latitude", {
      # equal u and v components should move a particle along a 45-degree rhumb line
      f <- field(5, 5, ymin = 55)
      tr <- wind_trails(f, cbind(-95, 60), hours = 1, steps = 1, direction = "downwind")
      p <- as.matrix(tr[, c("x", "y")])
      expect_lt(abs(geosphere::bearingRhumb(p[1, ], p[2, ]) - 45), 0.5)
})

test_that("trails leaving the domain end, unless wrapped into new trails", {
      f <- field(5, 0)
      tr <- wind_trails(f, cbind(-90.5, 35), hours = 20, steps = 10, direction = "downwind")
      expect_lt(max(tr$step), 10)
      trw <- wind_trails(f, cbind(-90.5, 35), hours = 20, steps = 10, direction = "downwind",
                         wrap = "horizontal")
      expect_equal(max(trw$step), 10)
      expect_true(all(trw$x >= -100 & trw$x <= -90))
      expect_gt(length(unique(trw$trail)), 1)  # wrapping starts a new trail
      expect_equal(unique(trw$particle), 1)
})

test_that("global fields wrap horizontally by default", {
      g <- terra::rast(nrows = 18, ncols = 36, xmin = -180, xmax = 180, ymin = -90, ymax = 90,
                       crs = "EPSG:4326", nlyrs = 2, vals = 0)
      g[[1]] <- 5
      f <- wind_field(g)
      seed <- cbind(175, 5)
      auto <- wind_trails(f, seed, hours = 48, steps = 20, direction = "downwind")
      expect_gt(length(unique(auto$trail)), 1)
      expect_equal(max(auto$step), 20)
      expect_equal(auto, wind_trails(f, seed, hours = 48, steps = 20, direction = "downwind",
                                     wrap = TRUE))
      off <- wind_trails(f, seed, hours = 48, steps = 20, direction = "downwind", wrap = FALSE)
      expect_lt(max(off$step), 20)
      # regional fields don't wrap by default
      r <- field(5, 0)
      expect_lt(max(wind_trails(r, cbind(-90.5, 35), hours = 20, steps = 10,
                                direction = "downwind")$step), 10)
})

test_that("wind_trails validates input", {
      f <- field(5, 0)
      expect_error(wind_trails(f, cbind(-95, 35)), "exactly one")
      expect_error(wind_trails(f, cbind(-95, 35), hours = 1, distance = 1), "exactly one")
      expect_error(wind_trails(f, cbind(-95, 35), hours = -1), "positive")
      expect_error(wind_trails(f, cbind(-95, 35), hours = 1, steps = 0), "steps")
      expect_error(wind_trails(f, cbind(-95, 35, 1), hours = 1), "two-column")
      expect_error(wind_trails(methods::as(f, "SpatRaster"), cbind(-95, 35), hours = 1), "wind_field")
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
      s <- wind_trails(f, cbind(c(-95, -94), c(35, 36)), hours = 3, steps = 6, sf = TRUE)
      expect_s3_class(s, "sf")
      expect_equal(nrow(s), 2)
      expect_named(s, c("trail", "particle", "geometry"))
      expect_true(all(sf::st_geometry_type(s) == "LINESTRING"))
})

test_that("wind_field accepts a time step taken from a wind_series", {
      ws <- windscape_example("wind_series")
      n <- ws@n_steps
      f <- wind_field(c(ws[[5]], ws[[n + 5]]))
      expect_s4_class(f, "wind_field")
      expect_equal(terra::values(f), terra::values(c(ws[[5]], ws[[n + 5]])), ignore_attr = TRUE)
})
