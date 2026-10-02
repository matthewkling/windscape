# wind_graph(), least_cost_distance(), least_cost_surface() ------------------------

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
      d <- least_cost_distance(wind_graph(r), xy, adjust = FALSE)
      dE <- geosphere::distGeo(c(0, xy[1, 2]), c(1, xy[1, 2])) # m per hop
      expect_equal(d[1, 2], 5 * dE / (5 * 3600), tolerance = 1e-8) # hours, 5 hops at 5 m/s
      expect_equal(d[2, 1], Inf)                                    # no westward wind
      expect_equal(unname(diag(d)), c(0, 0))
})

test_that("upwind graph is the transpose of the downwind graph", {
      r <- noisy_rose()
      xy <- terra::xyFromCell(r, c(5, 40, 77, 160))
      down <- least_cost_distance(wind_graph(r), xy, adjust = FALSE)
      up <- least_cost_distance(wind_graph(r, direction = "upwind"), xy, adjust = FALSE)
      expect_equal(up, t(down), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("rate = TRUE returns inverse cost distances", {
      r <- noisy_rose()
      xy <- terra::xyFromCell(r, c(5, 40, 77))
      g <- wind_graph(r)
      expect_equal(least_cost_distance(g, xy, adjust = FALSE, rate = TRUE),
                   1 / least_cost_distance(g, xy, adjust = FALSE))
})

test_that("adjust = FALSE and TRUE agree for sites at cell centers", {
      r <- noisy_rose()
      xy <- terra::xyFromCell(r, c(5, 40, 77))
      g <- wind_graph(r)
      expect_equal(least_cost_distance(g, xy, adjust = TRUE),
                   least_cost_distance(g, xy, adjust = FALSE), tolerance = 1e-8)
})

test_that("least_cost_surface agrees with least_cost_distance", {
      r <- noisy_rose()
      xy <- terra::xyFromCell(r, c(5, 40, 77))
      g <- wind_graph(r)
      s <- least_cost_surface(g, xy[1, , drop = FALSE])
      expect_s4_class(s, "SpatRaster")
      d <- least_cost_distance(g, xy, adjust = FALSE)
      expect_equal(terra::extract(s, xy)[, 1], d[1, ], tolerance = 1e-8)
      sr <- least_cost_surface(g, xy[1, , drop = FALSE], rate = TRUE)
      expect_equal(terra::values(sr), 1 / terra::values(s))
})

test_that("wrap connects the eastern and western edges", {
      r <- uv_rose(nr = 3, nc = 4, u = 5, v = 0)
      xy <- centers(r, 2, c(4, 1))
      open <- least_cost_distance(wind_graph(r), xy, adjust = FALSE)
      wrapped <- least_cost_distance(wind_graph(r, wrap = TRUE), xy, adjust = FALSE)
      dE <- geosphere::distGeo(c(0, xy[1, 2]), c(1, xy[1, 2]))
      expect_equal(open[1, 2], Inf)
      expect_equal(wrapped[1, 2], dE / (5 * 3600), tolerance = 1e-8) # one hop across the suture
})
