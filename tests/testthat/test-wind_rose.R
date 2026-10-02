# Grid checks and geometry in wind_series() / wind_rose() ------------------------

test_that("wind_series accepts longitude/latitude grids, including CRS-less ones", {
      expect_s4_class(wind_series(uv_raster()), "wind_series")
      r <- uv_raster()
      terra::crs(r) <- ""
      expect_s4_class(wind_series(r), "wind_series")
})

test_that("projected grids are rejected", {
      r <- uv_raster(xmin = 0, ymin = 0, xres = 10000, crs = "EPSG:5070")
      expect_error(wind_series(r), "longitude/latitude")
      expect_error(wind_series(uv_raster(crs = "local", xmin = 0, ymin = 0)), "longitude/latitude")
})

test_that("non-square cells are rejected", {
      expect_error(wind_series(uv_raster(xres = 1, yres = 0.5)), "square")
})

test_that("wind_rose re-checks grids that bypass wind_series", {
      r <- uv_raster(xmin = 0, ymin = 0, xres = 10000, crs = "EPSG:5070")
      x <- methods::as(r, "wind_series")
      x@n_steps <- 3
      expect_error(wind_rose(x), "longitude/latitude")
})

test_that("wind_series collates alternating u and v layers", {
      r <- uv_raster(u = 1, v = 2, n_steps = 2)       # u u v v
      alt <- r[[c(1, 3, 2, 4)]]                       # u v u v
      a <- wind_series(r, order = "uuvv")
      b <- wind_series(alt, order = "uvuv")
      expect_equal(terra::values(a), terra::values(b), ignore_attr = TRUE)
      expect_equal(b@n_steps, 2)
      expect_error(wind_series(r, order = "vvuu"))
})

test_that("conductance accounts for east-west cell width at each latitude", {
      # uniform 5 m/s westerly at latitudes 0.5 to 59.5
      r <- uv_raster(nr = 60, nc = 2, u = 5, v = 0, ymin = 0)
      w <- wind_rose(wind_series(r), trans = 1)
      e <- terra::values(w)[, "E"]
      lat <- terra::crds(w)[, 2]
      expect_true(all(abs(terra::values(w)[, c("SW", "W", "NW", "N", "NE", "SE", "S")]) < 1e-12))
      # E conductance = speed / east-west distance, so e * cos(lat) is ~constant
      # (up to ellipsoid effects)
      ec <- e * cos(lat * pi / 180)
      expect_lt(diff(range(ec)) / mean(ec), 0.01)
      expect_gt(e[which.max(lat)] / e[which.min(lat)], 1.9)
})

test_that("nearly square cells (within 1%) are accepted", {
      # CFSR's native Gaussian grid is about 0.316 x 0.317 degrees
      expect_s4_class(wind_series(uv_raster(xres = 0.3158, yres = 0.3175)), "wind_series")
})

test_that("numeric trans is stored as a working power function", {
      r <- as_wind_rose(methods::as(noisy_rose(), "SpatRaster"), trans = 2)
      expect_equal(r@trans(3), 9)
      expect_equal(as_wind_rose(methods::as(noisy_rose(), "SpatRaster"), trans = sqrt)@trans(9), 3)
})
