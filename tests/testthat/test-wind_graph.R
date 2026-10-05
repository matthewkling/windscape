# wind_graph(), pairwise_least_cost(), least_cost_surface() ------------------------

# cell-center coordinates of row `row`, columns `cols`, of raster r
centers <- function(r, row, cols) terra::xyFromCell(r, terra::cellFromRowCol(r, row, cols))

test_that("wind_graph returns a downwind or upwind wind_graph", {
      r <- uv_rose()
      g <- wind_graph(r)
      expect_s4_class(g, "wind_graph")
      expect_s4_class(g, "TransitionLayer")
      expect_equal(g@direction, "downwind")
      expect_equal(wind_graph(r, direction = "upwind")@direction, "upwind")
})

test_that("downwind travel time along a uniform westerly is exact", {
      r <- uv_rose(nr = 3, nc = 6, u = 5, v = 0)
      xy <- centers(r, 2, c(1, 6))
      d <- pairwise_least_cost(wind_graph(r), xy, snap = TRUE)
      dE <- geosphere::distGeo(c(0, xy[1, 2]), c(1, xy[1, 2])) # m per hop
      expect_equal(d[1, 2], 5 * dE / (5 * 3600), tolerance = 1e-8) # hours, 5 hops at 5 m/s
      expect_equal(d[2, 1], Inf)                                    # no westward wind
      expect_equal(unname(diag(d)), c(0, 0))
})

test_that("upwind graph is the transpose of the downwind graph", {
      r <- noisy_rose()
      xy <- terra::xyFromCell(r, c(5, 40, 77, 160))
      down <- pairwise_least_cost(wind_graph(r), xy, snap = TRUE)
      up <- pairwise_least_cost(wind_graph(r, direction = "upwind"), xy, snap = TRUE)
      expect_equal(up, t(down), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("rate = TRUE returns inverse cost distances", {
      r <- noisy_rose()
      xy <- terra::xyFromCell(r, c(5, 40, 77))
      g <- wind_graph(r)
      expect_equal(pairwise_least_cost(g, xy, snap = TRUE, rate = TRUE),
                   1 / pairwise_least_cost(g, xy, snap = TRUE))
})

test_that("least_cost_surface agrees with pairwise_least_cost", {
      r <- noisy_rose()
      xy <- terra::xyFromCell(r, c(5, 40, 77))
      g <- wind_graph(r)
      s <- least_cost_surface(g, xy[1, , drop = FALSE])
      expect_s4_class(s, "SpatRaster")
      expect_equal(names(s), "hours")
      d <- pairwise_least_cost(g, xy, snap = TRUE)
      expect_equal(terra::extract(s, xy)[, 1], d[1, ], tolerance = 1e-8)
      sr <- least_cost_surface(g, xy[1, , drop = FALSE], rate = TRUE)
      expect_equal(names(sr), "rate")
      expect_equal(terra::values(sr), 1 / terra::values(s), ignore_attr = TRUE)
})

test_that("wrap connects the eastern and western edges", {
      r <- uv_rose(nr = 3, nc = 4, u = 5, v = 0)
      xy <- centers(r, 2, c(4, 1))
      open <- pairwise_least_cost(wind_graph(r), xy, snap = TRUE)
      wrapped <- pairwise_least_cost(wind_graph(r, wrap = TRUE), xy, snap = TRUE)
      dE <- geosphere::distGeo(c(0, xy[1, 2]), c(1, xy[1, 2]))
      expect_equal(open[1, 2], Inf)
      expect_equal(wrapped[1, 2], dE / (5 * 3600), tolerance = 1e-8) # one hop across the suture
})

# least_cost_paths() ------------------------

test_that("path costs equal least-cost distances, with vertices from origin to destination", {
      r <- noisy_rose()
      g <- wind_graph(r)
      from <- terra::xyFromCell(r, c(5, 77))
      to <- terra::xyFromCell(r, c(40, 120, 160))
      p <- least_cost_paths(g, from, to)
      expect_named(p, c("trail", "from", "to", "step", "hours", "x", "y"))
      expect_equal(length(unique(p$trail)), 6)                    # all 2 x 3 pairs
      cost <- matrix(gdistance::costDistance(g, from, to), 2, 3)
      fin <- p[!duplicated(p$trail, fromLast = TRUE), ]
      expect_equal(fin$hours, cost[cbind(fin$from, fin$to)], tolerance = 1e-8)
      for(d in split(p, p$trail)){
            expect_equal(terra::cellFromXY(r, as.matrix(d[1, c("x", "y")])), terra::cellFromXY(r, from[d$from[1], , drop = FALSE]))
            expect_equal(terra::cellFromXY(r, as.matrix(d[nrow(d), c("x", "y")])), terra::cellFromXY(r, to[d$to[1], , drop = FALSE]))
            expect_false(is.unsorted(d$hours))
            expect_equal(d$step, seq_len(nrow(d)) - 1)
      }
})

test_that("a uniform westerly gives straight eastward paths; westward pairs are omitted", {
      r <- uv_rose(nr = 3, nc = 6, u = 5, v = 0)
      g <- wind_graph(r)
      west <- terra::xyFromCell(r, terra::cellFromRowCol(r, 2, 1))
      east <- terra::xyFromCell(r, terra::cellFromRowCol(r, 2, 6))
      p <- least_cost_paths(g, west, east)
      expect_equal(nrow(p), 6)
      expect_true(all(diff(p$x) > 0))
      expect_equal(unique(p$y), unname(west[, 2]))
      expect_warning(q <- least_cost_paths(g, east, west), "no path")
      expect_equal(nrow(q), 0)
})

test_that("paths are directional", {
      r <- noisy_rose()
      g <- wind_graph(r)
      a <- terra::xyFromCell(r, 5)
      b <- terra::xyFromCell(r, 160)
      ab <- least_cost_paths(g, a, b)
      ba <- least_cost_paths(g, b, a)
      expect_false(isTRUE(all.equal(max(ab$hours), max(ba$hours))))
})

test_that("nearest and matched pair options", {
      r <- noisy_rose()
      g <- wind_graph(r)
      from <- terra::xyFromCell(r, c(5, 77, 100))
      to <- terra::xyFromCell(r, c(40, 120, 160))
      cost <- matrix(gdistance::costDistance(g, from, to), 3, 3)
      n <- least_cost_paths(g, from, to, pairs = "nearest")
      chosen <- unique(n[, c("from", "to")])
      expect_equal(nrow(chosen), 3)
      expect_equal(chosen$to, unname(apply(cost, 1, which.min)))
      m <- least_cost_paths(g, from, to, pairs = "matched")
      expect_equal(unique(m[, c("from", "to")]), data.frame(from = 1:3, to = 1:3), ignore_attr = TRUE)
      expect_error(least_cost_paths(g, from, to[1:2, ], pairs = "matched"), "same number")
})

test_that("input handling: data frames, SpatVectors, same-cell pairs, validation", {
      r <- noisy_rose()
      g <- wind_graph(r)
      from <- terra::xyFromCell(r, 5)
      to <- terra::xyFromCell(r, c(5, 160))     # first destination is the origin's cell
      p <- least_cost_paths(g, as.data.frame(from), terra::vect(to, crs = "EPSG:4326"))
      expect_equal(unique(p$to), 2)
      expect_error(least_cost_paths(r, from, to), "wind_graph")
      expect_error(least_cost_paths(g, cbind(1, 2, 3), to), "two-column")
      expect_error(least_cost_paths(g, from, to, pairs = "some"))
})

test_that("paths draw with geom_wind_trail", {
      r <- noisy_rose()
      g <- wind_graph(r)
      p <- least_cost_paths(g, terra::xyFromCell(r, 5), terra::xyFromCell(r, c(40, 160)))
      plt <- ggplot2::ggplot(p, ggplot2::aes(x, y)) + geom_wind_trail(ggplot2::aes(color = hours))
      ld <- ggplot2::layer_data(plt)
      expect_equal(length(unique(ld$group)), 2)
      expect_s3_class(ggplot2::ggplotGrob(plt), "gtable")
})
