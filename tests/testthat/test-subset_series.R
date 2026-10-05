# wind_times() and subset_series() ---------------------------------------------------------

series <- windscape_example("wind_series")
n <- series@n_steps
times <- wind_times(series)
uv <- function(x) terra::values(as(x, "SpatRaster"))

test_that("wind_times parses layer names, including midnight written without a time", {
      expect_s3_class(times, "POSIXct")
      expect_length(times, n)
      expect_false(anyNA(times))
      expect_equal(times[1:2], as.POSIXct(c("2000-01-01 00:00:00", "2000-01-01 06:00:00"), tz = "UTC"))
})

test_that("wind_times falls back to terra::time, and errors without times", {
      r <- as(series, "SpatRaster")
      names(r) <- paste0("layer", seq_len(terra::nlyr(r)))
      expect_error(wind_times(wind_series(r)), "no times found")
      terra::time(r) <- rep(times, 2)
      expect_equal(wind_times(wind_series(r)), times)
})

test_that("subset_series keeps u and v of the selected steps together, with their values", {
      s <- subset_series(series, months = 6:8)
      keep <- which(as.integer(format(times, "%m")) %in% 6:8)
      expect_s4_class(s, "wind_series")
      expect_equal(s@n_steps, length(keep))
      expect_equal(names(s), names(series)[c(keep, n + keep)])
      expect_equal(uv(s), uv(series)[, c(keep, n + keep)], ignore_attr = TRUE)
})

test_that("criteria select the right steps and combine", {
      h <- as.integer(format(times, "%H"))
      m <- as.integer(format(times, "%m"))
      expect_equal(wind_times(subset_series(series, hours = c(0, 12))), times[h %in% c(0, 12)])
      expect_equal(wind_times(subset_series(series, months = 1, hours = 6)), times[m == 1 & h == 6])
      expect_equal(wind_times(subset_series(series, steps = 3:5)), times[3:5])
      expect_equal(wind_times(subset_series(series, steps = h == 18)), times[h == 18])
      expect_equal(wind_times(subset_series(series, steps = 1:10, hours = 0)), times[1:10][h[1:10] == 0])
})

test_that("start and end are inclusive, and a bare end date includes the whole day", {
      s <- subset_series(series, start = "2000-01-15", end = "2000-02-01")
      expect_equal(wind_times(s), times[times >= as.POSIXct("2000-01-15", tz = "UTC") &
                                              times < as.POSIXct("2000-02-02", tz = "UTC")])
      expect_true(any(format(wind_times(s), "%H") == "18")) # the last day's evening is included
      s2 <- subset_series(series, end = as.POSIXct("2000-01-01 06:00:00", tz = "UTC"))
      expect_equal(wind_times(s2), times[1:2])
      s3 <- subset_series(series, start = as.Date("2000-12-15"))
      expect_equal(wind_times(s3), times[times >= as.POSIXct("2000-12-15", tz = "UTC")])
})

test_that("subset_series rejects bad input", {
      expect_error(subset_series(series, months = 13), "months")
      expect_error(subset_series(series, hours = 24), "hours")
      expect_error(subset_series(series, steps = c(TRUE, FALSE)), "one element per time step")
      expect_error(subset_series(series, steps = n + 1), "steps")
      expect_error(subset_series(series, start = "not a date"), "could not be read")
      expect_error(subset_series(series, start = "2001-01-01"), "no time steps")
      expect_error(subset_series(as(series, "SpatRaster")), "wind_series")
})

# write the example series as two files, January-June and July-December
split_files <- function(){
      dir <- withr::local_tempdir(.local_envir = parent.frame())
      m <- as.integer(format(times, "%m"))
      vapply(list(h1 = m <= 6, h2 = m > 6), function(k){
            f <- file.path(dir, paste0(if(k[1]) "h1" else "h2", ".tif"))
            terra::writeRaster(as(subset_series(series, steps = k), "SpatRaster"), f)
            f
      }, character(1))
}

test_that("roses from selected time steps of files match roses from the whole series", {
      files <- split_files()
      expect_equal(wind_series(files)@n_steps, n)
      for(sel in list(list(months = 7:8), list(hours = c(0, 6)))){
            r <- wind_rose(do.call(subset_series, c(list(wind_series(files)), sel)), trans = 2)
            ref <- wind_rose(do.call(subset_series, c(list(series), sel)), trans = 2)
            expect_equal(r@n_steps, ref@n_steps)
            expect_equal(terra::values(r), terra::values(ref), tolerance = 1e-6)
      }
})
