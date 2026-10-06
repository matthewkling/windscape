# windscape_example() ------------------------

test_that("example wind_series loads as a 96-step wind_series", {
      ws <- windscape_example("wind_series")
      expect_s4_class(ws, "wind_series")
      expect_equal(ws@n_steps, 96)
      expect_equal(terra::nlyr(ws), 192)
      expect_true(isTRUE(terra::is.lonlat(ws)))
      expect_equal(unname(as.vector(terra::ext(ws))), c(-120.1579, -89.84211, 29.84127, 50.15873),
                   tolerance = 1e-5)
      expect_false(anyNA(terra::values(ws)))
})

test_that("example wind_rose loads and agrees with a rose built from the example series", {
      rose <- windscape_example("wind_rose")
      expect_s4_class(rose, "wind_rose")
      expect_equal(names(rose), c("SW", "W", "NW", "N", "NE", "E", "SE", "S"))
      expect_equal(rose@n_steps, 576)
      expect_true(all(terra::values(rose) >= 0))
      expect_equal(rose@trans(2), 2)
      # the 96-step series is a subsample of the 576 hours behind the rose
      sub <- wind_rose(windscape_example("wind_series"), trans = 1)
      expect_gt(cor(as.vector(terra::values(sub)), as.vector(terra::values(rose))), 0.98)
})

test_that("example wind_field is Hurricane Katrina", {
      f <- windscape_example("wind_field")
      expect_s4_class(f, "wind_field")
      expect_equal(terra::nlyr(f), 2)
      speed <- sqrt(sum(f^2))
      expect_gt(terra::global(speed, "max")[[1]], 33) # hurricane-force winds (64 knots)
})

test_that("example rose works in the connectivity functions", {
      rose <- windscape_example("wind_rose")
      xy <- cbind(c(-110, -95), c(40, 40))
      d <- pairwise_least_cost(rose, xy)
      expect_true(all(is.finite(d)))
      w <- quietly(random_walk(rose, xy[1, , drop = FALSE], mode = "stream", half_life = 24))
      expect_gt(terra::global(w[["deposition"]], "sum", na.rm = TRUE)[[1]], 0)
})

test_that("windscape_example rejects unknown names", {
      expect_error(windscape_example("nope"))
})
