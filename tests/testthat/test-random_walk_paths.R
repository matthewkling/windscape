# random_walk_paths() ---------------------------------------------------------------------

rose <- windscape_example("wind_rose")
site <- cbind(-105, 40)
first <- function(p) p[!duplicated(p$trail), ]
last <- function(p) p[!duplicated(p$trail, fromLast = TRUE), ]

test_that("downwind paths run from the source, deterministically", {
      p <- quietly(random_walk_paths(rose, site, n = 20, half_life = 48))
      expect_named(p, c("trail", "step", "x", "y"))
      expect_equal(length(unique(p$trail)), 20)
      expect_true(all(first(p)$x == site[1] & first(p)$y == site[2]))
      expect_true(all(tapply(p$step, p$trail, function(s) all(diff(s) == 1))))
      expect_identical(p, quietly(random_walk_paths(rose, site, n = 20, half_life = 48)))
})

test_that("path ends are placed in proportion to deposition", {
      n <- 40
      p <- quietly(random_walk_paths(rose, site, n = n, half_life = 48))
      w <- quietly(random_walk(rose, site, mode = "stream", half_life = 48, density = FALSE))
      dep <- terra::values(w$deposition)[, 1]
      near <- terra::adjacent(rose, terra::cellFromXY(rose, site), directions = "queen", include = TRUE)
      share <- sum(dep[near], na.rm = TRUE) / sum(dep, na.rm = TRUE)
      ends <- terra::cellFromXY(rose, as.matrix(last(p)[, c("x", "y")]))
      expect_lt(abs(mean(ends %in% near) - share), 1.5 / n) # within sampling resolution
})

test_that("with no deposition, paths end at the domain edges", {
      p <- quietly(random_walk_paths(rose, site, n = 15))
      rc <- terra::rowColFromCell(rose, terra::cellFromXY(rose, as.matrix(last(p)[, c("x", "y")])))
      expect_true(all(rc[, 1] %in% c(1, terra::nrow(rose)) | rc[, 2] %in% c(1, terra::ncol(rose))))
})

test_that("upwind paths run to the receptor", {
      p <- quietly(random_walk_paths(rose, site, n = 15, direction = "upwind", half_life = 48))
      expect_equal(length(unique(p$trail)), 15)
      expect_true(all(last(p)$x == site[1] & last(p)$y == site[2]))
})

test_that("multiple sources and raster sources work", {
      sites <- cbind(c(-110, -100), c(38, 44))
      p <- quietly(random_walk_paths(rose, sites, n = 20, half_life = 48))
      s <- first(p)
      expect_true(all(paste(s$x, s$y) %in% paste(sites[, 1], sites[, 2])))
      expect_true(all(paste(sites[, 1], sites[, 2]) %in% paste(s$x, s$y))) # both sources used

      r <- terra::rast(rose[[1]])
      terra::values(r) <- 0
      r[terra::cellFromXY(r, site)] <- 1
      pr <- quietly(random_walk_paths(rose, r, n = 10, half_life = 48))
      expect_true(all(terra::cellFromXY(rose, as.matrix(first(pr)[, c("x", "y")])) ==
                            terra::cellFromXY(rose, site)))
})

test_that("random_walk_paths rejects arguments it controls", {
      expect_error(random_walk_paths(rose, site, mode = "pulse"), "can't be supplied")
      expect_error(random_walk_paths(rose, site, flux = FALSE), "can't be supplied")
      expect_error(random_walk_paths(rose, site, n = 0), "positive integer")
      expect_error(random_walk_paths(as(rose, "SpatRaster"), site), "wind_rose")
})

test_that("paths can be traced to chosen points", {
      to <- cbind(c(-95, -100, -98), c(45, 32, 41))
      p <- quietly(random_walk_paths(rose, site, to = to, half_life = 48))
      expect_equal(sort(unique(p$trail)), 1:3)
      expect_true(all(first(p)$x == site[1] & first(p)$y == site[2]))
      e <- last(p)
      expect_equal(unname(as.matrix(e[order(e$trail), c("x", "y")])), to)
      expect_error(random_walk_paths(rose, site, to = cbind(0, 0)), "outside the extent")
      expect_error(random_walk_paths(rose, site, to = 1:3), "two-column")
})
