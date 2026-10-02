# ws_summarize() ------------------------

ws_grid <- function(){
      terra::rast(nrows = 21, ncols = 21, xmin = -110.5, xmax = -89.5, ymin = 29.5, ymax = 50.5,
                  crs = "EPSG:4326", vals = 0)
}
origin <- matrix(c(-100, 40), 1)

test_that("a single accessible cell defines the centroid and has zero dispersion", {
      x <- ws_grid()
      x[terra::cellFromXY(x, cbind(-95, 40))] <- 1
      s <- ws_summarize(x, origin)
      expect_equal(unname(s[c("centroid_x", "centroid_y")]), c(-95, 40))
      expect_equal(unname(s["centroid_distance"]),
                   geosphere::distGeo(origin, c(-95, 40)) / 1000, tolerance = 1e-6)
      expect_equal(unname(s["centroid_bearing"]), geosphere::bearing(origin, c(-95, 40)), tolerance = 1e-6)
      expect_lt(abs(unname(s["windshed_bearing"]) - 90), 2)
      expect_equal(unname(s["windshed_isotropy"]), 0, tolerance = 1e-6)
})

test_that("mass to the north gives a northward windshed", {
      x <- ws_grid()
      x[terra::cellFromXY(x, cbind(c(-101, -100, -99), 46))] <- 1
      s <- ws_summarize(x, origin)
      expect_lt(abs(unname(s["centroid_bearing"])), 0.5)
      expect_gt(unname(s["windshed_isotropy"]), 0)
})

test_that("radius excludes distant cells", {
      x <- ws_grid()
      x[terra::cellFromXY(x, cbind(c(-97, -91), 40))] <- 1
      s_all <- ws_summarize(x, origin)
      s_near <- ws_summarize(x, origin, radius = 500)
      expect_gt(unname(s_all["windshed_distance"]), unname(s_near["windshed_distance"]))
      expect_equal(unname(s_near["centroid_x"]), -97)
})

test_that("summary statistics are named", {
      x <- ws_grid()
      x[] <- 1
      expect_named(ws_summarize(x, origin),
                   c("centroid_x", "centroid_y", "centroid_distance", "centroid_bearing",
                     "windshed_distance", "windshed_bearing", "windshed_isotropy",
                     "windshed_size", "windshed_landarea"))
})
