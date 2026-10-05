# [[, subset(), and c() return plain SpatRasters for windscape classes ---------------------------

series <- windscape_example("wind_series")
rose <- windscape_example("wind_rose")
field <- wind_field(series, step = 1)
is_plain <- function(x) expect_identical(class(x), structure("SpatRaster", package = "terra"))

test_that("layer selection and combination drop the windscape class", {
      for(x in list(series, rose, field)){
            is_plain(x[[1]])
            is_plain(terra::subset(x, 1))
            is_plain(c(x, x))
      }
      is_plain(c(as(series, "SpatRaster"), series))
})

test_that("values are unchanged by the plain-raster methods", {
      v <- terra::values(as(series, "SpatRaster"))
      expect_equal(terra::values(series[[3:5]]), v[, 3:5])
      expect_equal(terra::values(terra::subset(series, 3:5)), v[, 3:5])
      expect_equal(terra::nlyr(c(series, rose)), terra::nlyr(series) + 8)
      expect_equal(terra::values(rose[["N"]]), terra::values(as(rose, "SpatRaster"))[, "N", drop = FALSE])
})

test_that("operations that keep the layers keep the class", {
      expect_s4_class(series * 2, "wind_series")
      expect_s4_class(terra::crop(series, terra::ext(-110, -100, 35, 45)), "wind_series")
      expect_s4_class(rose / 2, "wind_rose")
})

test_that("grid-changing operations are refused for wind roses", {
      expect_error(terra::aggregate(rose, 2), "downscale")
      expect_error(terra::disagg(rose, 2), "downscale")
      expect_error(terra::resample(rose, terra::rast(rose)), "downscale")
      expect_error(terra::project(rose, "EPSG:3857"), "downscale")
      # wind series can still be aggregated, then made into a rose
      coarse <- terra::aggregate(series, 2)
      expect_s4_class(coarse, "wind_series")
      expect_s4_class(wind_rose(coarse), "wind_rose")
})

test_that("downscale still works, keeping the rose's metadata", {
      d <- downscale(rose, 2)
      expect_s4_class(d, "wind_rose")
      expect_equal(terra::ncell(d), 4 * terra::ncell(rose))
      expect_equal(d@n_steps, rose@n_steps)
      expect_equal(d@trans(3), rose@trans(3))
})
