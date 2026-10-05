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

test_that("wind_rose builds in chunks with the same result as all at once", {
      series <- windscape_example("wind_series")
      whole <- wind_rose(series, trans = 2)
      chunked <- withr::with_options(list(windscape.chunk_values = 2 * terra::ncell(series) * 7),
                                     wind_rose(series, trans = 2)) # 7 steps per chunk
      expect_equal(chunked@n_steps, whole@n_steps)
      expect_equal(terra::values(chunked), terra::values(whole), tolerance = 1e-10)
      expect_equal(chunked@trans(3), 9)
})

test_that("wind_rose loads saved roses from rasters and files", {
      r <- wind_rose(windscape_example("wind_series"), trans = 2)
      f <- withr::local_tempfile(fileext = ".tif")
      terra::writeRaster(r, f)
      loaded <- wind_rose(f, trans = 2, n_steps = r@n_steps)
      expect_s4_class(loaded, "wind_rose")
      expect_equal(names(loaded), c("SW", "W", "NW", "N", "NE", "E", "SE", "S"))
      expect_equal(terra::values(loaded), terra::values(r), tolerance = 1e-6)
      expect_equal(loaded@n_steps, r@n_steps)
      expect_equal(loaded@trans(3), 9)
      expect_s4_class(combine_roses(r, loaded), "wind_rose")

      plain <- terra::rast(f)
      names(plain) <- paste0("lyr", 1:8)
      expect_equal(terra::values(wind_rose(plain)), terra::values(loaded), ignore_attr = TRUE)
      expect_identical(wind_rose(r), r)
})

test_that("wind_rose writes to a file when given a filename", {
      f <- withr::local_tempfile(fileext = ".tif")
      r <- wind_rose(windscape_example("wind_series"), filename = f)
      expect_true(file.exists(f))
      expect_s4_class(r, "wind_rose")
      expect_equal(terra::values(wind_rose(f)), terra::values(r), tolerance = 1e-6)
})

test_that("wind_rose rejects inputs that aren't series or roses", {
      series <- windscape_example("wind_series")
      four <- as(subset_series(series, steps = 1:4), "SpatRaster") # 8 layers, but wind data
      expect_error(wind_rose(four), "use wind_series\\(\\) first")
      expect_error(wind_rose(as(series, "SpatRaster")[[1:6]]), "must have 8 layers|use wind_series")
      expect_error(wind_rose(c("a.tif", "b.tif")), "single file path")
      expect_error(wind_rose("no/such/file.tif"), "not found")
      expect_error(wind_rose(1:8), "must be a wind_series")
})
