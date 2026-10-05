# wind_field() and mean() for wind_series ----------------------------------------------------

series <- windscape_example("wind_series")
n <- series@n_steps
r <- as(series, "SpatRaster")
vals <- function(x) terra::values(as(x, "SpatRaster"))

test_that("wind_field takes one time step of a wind_series", {
      f <- wind_field(series, step = 3)
      expect_s4_class(f, "wind_field")
      expect_equal(vals(f), vals(r[[c(3, n + 3)]]), ignore_attr = TRUE)
      one <- subset_series(series, steps = 5)
      expect_equal(vals(wind_field(one)), vals(r[[c(5, n + 5)]]), ignore_attr = TRUE)
})

test_that("wind_field reads two-layer rasters and files", {
      expect_s4_class(wind_field(r[[c(1, n + 1)]]), "wind_field")
      f <- withr::local_tempfile(fileext = ".tif")
      terra::writeRaster(r[[c(1, n + 1)]], f)
      expect_equal(vals(wind_field(f)), vals(r[[c(1, n + 1)]]), ignore_attr = TRUE, tolerance = 1e-6)
})

test_that("wind_field rejects ambiguous or invalid input", {
      expect_error(wind_field(series), "choose one with `step`")
      expect_error(wind_field(series, step = n + 1), "between 1 and")
      expect_error(wind_field(series, step = 1:2), "single integer")
      expect_error(wind_field(r[[c(1, n + 1)]], step = 1), "only to a wind_series")
      expect_error(wind_field(r[[1:4]]), "use wind_series")
      expect_error(wind_field("no/such/file.tif"), "not found")
})

test_that("mean() of a wind_series is the time-mean wind field", {
      m <- mean(series)
      expect_s4_class(m, "wind_field")
      expect_equal(names(m), c("u", "v"))
      v <- vals(r)
      expect_equal(vals(m)[, 1], rowMeans(v[, seq_len(n)]), tolerance = 1e-6)
      expect_equal(vals(m)[, 2], rowMeans(v[, n + seq_len(n)]), tolerance = 1e-6)

      # over selected steps, and with missing values
      s <- subset_series(series, months = 7)
      expect_equal(vals(mean(s))[, 1], rowMeans(vals(s)[, seq_len(s@n_steps)]), tolerance = 1e-6)
      gap <- r
      gap[[1]][1] <- NA
      expect_true(is.na(vals(mean(wind_series(gap)))[1, 1]))
      expect_false(is.na(vals(mean(wind_series(gap), na.rm = TRUE))[1, 1]))
})

test_that("mean wind and net flow point in nearly the same direction with trans = 1", {
      m <- ggplot2::fortify(mean(series))
      nf <- ggplot2::fortify(net_flow(wind_rose(series)))
      strong <- m$speed > stats::quantile(m$speed, 0.5)
      dif <- abs(((m$bearing - nf$bearing + 180) %% 360) - 180)[strong]
      expect_lt(stats::median(dif), 10)
})
