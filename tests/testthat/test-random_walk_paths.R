# random_walk_paths() ---------------------------------------------------------------------

rose <- windscape_example("wind_rose")
site <- cbind(-105, 40)
first <- function(p) p[!duplicated(p$path), ]
last <- function(p) p[!duplicated(p$path, fromLast = TRUE), ]

test_that("downwind paths run from the source, deterministically", {
      p <- quietly(random_walk_paths(rose, site, n = 20, half_life = 48))
      expect_named(p, c("trail", "path", "step", "x", "y"))
      expect_equal(p$trail, p$path) # no seam on a regional grid
      expect_equal(length(unique(p$path)), 20)
      expect_true(all(first(p)$x == site[1] & first(p)$y == site[2]))
      expect_true(all(tapply(p$step, p$path, function(s) all(diff(s) == 1))))
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
      expect_equal(length(unique(p$path)), 15)
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
      expect_equal(sort(unique(p$path)), 1:3)
      expect_true(all(first(p)$x == site[1] & first(p)$y == site[2]))
      e <- last(p)
      expect_equal(unname(as.matrix(e[order(e$path), c("x", "y")])), to)
      expect_error(random_walk_paths(rose, site, to = cbind(0, 0)), "outside the extent")
      expect_error(random_walk_paths(rose, site, to = 1:3), "two-column")
})

test_that("on global grids, paths cross the seam, splitting into trails there", {
      g <- global_rose()
      src <- cbind(155, 5)
      p0 <- quietly(random_walk_paths(g, src, n = 30, half_life = 400, wrap = FALSE))
      p <- quietly(random_walk_paths(g, src, n = 30, half_life = 400))
      expect_equal(length(unique(p$path)), 30)
      expect_true(all(last(p0)$x > 0)) # without wrap, material can't get past the east edge
      expect_gt(sum(last(p)$x < 0), 0)  # with wrap (the default here), it continues across it
      expect_true(all(first(p)$x == 155))
      # trails never jump across the map; paths that cross the seam have two or more trails
      expect_true(all(tapply(p$x, p$trail, function(x) all(abs(diff(x)) < 180))))
      expect_gt(length(unique(p$trail)), 30)
      expect_equal(p$trail, seam_trails(p$path, p$x, 360))
      # a path to a point across the seam runs east from the source in two pieces
      to <- quietly(random_walk_paths(g, src, to = cbind(-165, 5), half_life = 400))
      expect_equal(to$x[1], 155)
      expect_equal(to$x[nrow(to)], -165)
      expect_equal(length(unique(to$trail)), 2)
      expect_equal(unique(to$path), 1)
      # upwind paths to a receptor west of the seam start east of it
      up <- quietly(random_walk_paths(g, cbind(-165, 5), n = 10, half_life = 400,
                                      direction = "upwind"))
      expect_gt(sum(first(up)$x > 0), 0)
      expect_true(all(last(up)$x == -165))
})
