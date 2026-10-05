# wind_series() ---------------------------------------------------------------------------

series <- windscape_example("wind_series")
n <- series@n_steps
r <- as(series, "SpatRaster")
uv <- function(x) terra::values(as(x, "SpatRaster"))
halves <- list(steps = list(1:48, 49:n))

# a piece of the series as a raster, in the given layer order
piece <- function(steps, order = "uuvv"){
      lyr <- if(order == "uuvv") c(steps, n + steps) else as.vector(rbind(steps, n + steps))
      r[[lyr]]
}

test_that("wind_series builds a series from one raster in either order", {
      expect_equal(uv(wind_series(r)), uv(series))
      w <- wind_series(piece(1:n, "uvuv"), order = "uvuv")
      expect_equal(uv(w), uv(series), ignore_attr = TRUE)
      expect_equal(names(w), names(series))
      expect_equal(w@n_steps, n)
})

test_that("wind_series combines a list of rasters, or files, into one series", {
      a <- 1:48
      b <- 49:n
      expect_equal(uv(wind_series(list(piece(a), piece(b)))), uv(series), ignore_attr = TRUE)
      expect_equal(uv(wind_series(list(piece(a, "uvuv"), piece(b, "uvuv")), order = "uvuv")),
                   uv(series), ignore_attr = TRUE)

      dir <- withr::local_tempdir()
      files <- file.path(dir, c("a.tif", "b.tif"))
      terra::writeRaster(piece(a), files[1])
      terra::writeRaster(piece(b), files[2])
      w <- wind_series(files)
      expect_s4_class(w, "wind_series")
      expect_equal(w@n_steps, n)
      expect_equal(names(w), names(series))
      expect_equal(uv(w), uv(series), ignore_attr = TRUE, tolerance = 1e-6)
})

test_that("wind_series rejects bad input", {
      expect_error(wind_series(r[[1:3]]), "even number")
      expect_error(wind_series(list(piece(1:10), piece(11:13)[[1:5]])), "input 2")
      shifted <- terra::shift(piece(11:20), dx = 1)
      expect_error(wind_series(list(piece(1:10), shifted)), "different grids")
      expect_error(wind_series(list(piece(1:10), "not a raster")), "SpatRasters")
      expect_error(wind_series(character(0)), "at least one")
})

test_that("wind_series accepts named inputs", {
      w <- wind_series(list(first = piece(1:48), second = piece(49:n)))
      expect_equal(w@n_steps, n)
})

test_that("malformed series are rejected, not silently misused", {
      expect_false(inherits(series[[1:10]], "wind_series")) # layer selection drops the class
      bad <- series
      bad@n_steps <- 10 # e.g. a series modified by other means
      expect_error(wind_rose(bad), "malformed wind_series")
      expect_error(mean(bad), "malformed wind_series")
      expect_error(subset_series(bad, steps = 1), "malformed wind_series")
      expect_error(wind_times(bad), "malformed wind_series")
      expect_error(ggplot2::fortify(bad), "malformed wind_series")
})
