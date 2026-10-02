# point_distance(), cell_distance(), check_cell_distance() ------------------------

test_that("point_distance returns symmetric geodesic distances in km", {
      ll <- cbind(c(-100, -99, -95), c(40, 41, 40))
      d <- point_distance(ll)
      expect_equal(dim(d), c(3, 3))
      expect_equal(d, t(d))
      expect_equal(diag(d), rep(0, 3))
      expect_equal(d[1, 3], geosphere::distGeo(ll[1, ], ll[3, ]) / 1000)
})

test_that("cell_distance measures between the centers of the cells containing points", {
      r <- noisy_rose()
      ll <- cbind(c(-99.8, -91.3, -88.6), c(31.2, 36.9, 33.4))
      cc <- terra::xyFromCell(r, terra::cellFromXY(r, ll))
      expect_equal(cell_distance(r, ll), point_distance(cc), ignore_attr = TRUE)
})

test_that("check_cell_distance reports and optionally returns distance ratios", {
      r <- noisy_rose()
      ll <- cbind(c(-99.8, -91.3, -88.6, -88.5), c(31.2, 36.9, 33.4, 33.45))
      msgs <- capture_messages(check_cell_distance(r, ll))
      expect_true(any(grepl("Total point pairs: 6", msgs)))
      expect_true(any(grepl("same grid cell: 1", msgs)))
      out <- suppressMessages(check_cell_distance(r, ll, return = TRUE))
      expect_equal(out, cell_distance(r, ll) / point_distance(ll))
})
