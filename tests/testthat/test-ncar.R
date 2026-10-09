# download_wind_data(), wind_series(), download_land_mask(), combine_roses() --------------------------
# Downloads are tested against a fake server (helper-ncar.R); no network access is needed.

# value of a raster layer at a lon/lat point
at <- function(r, lon, lat) terra::extract(r, cbind(lon, lat))[1, 1]

test_that("bounding boxes are split at the 0/360 seam", {
      expect_equal(check_bbox(c(-120, -90), c(30, 50))$pieces, list(c(240, 270)))
      expect_equal(check_bbox(c(-10, 10), c(30, 50))$pieces, list(c(350, 360), c(0, 10)))
      expect_equal(check_bbox(c(170, 200), c(30, 50))$pieces, list(c(170, 200)))
      expect_equal(check_bbox(c(-180, 180), c(30, 50))$pieces, list(c(0, 360)))
      expect_equal(check_bbox(c(0, 10), c(30, 50))$pieces, list(c(0, 10)))
      expect_equal(check_bbox(c(-10, 0), c(30, 50))$pieces, list(c(350, 360)))
      expect_equal(check_bbox(c(-120, -90), c(30, 50))$convention, "-180_180")
      expect_equal(check_bbox(c(240, 270), c(30, 50))$convention, "0_360")
      expect_error(check_bbox(c(-10, 200), c(30, 50)), "mixes")
      expect_error(check_bbox(c(0, 10), c(30, 95)), "ylim")
})

test_that("NetCDF time units are parsed in both ERA5 and CFSR styles", {
      expect_equal(parse_nc_time(c(0, 6), "hours since 1900-01-01 00:00:00"),
                   as.POSIXct(c("1900-01-01 00:00:00", "1900-01-01 06:00:00"), tz = "UTC"))
      expect_equal(parse_nc_time(1, "Hour since 2000-04-01T00:00:00Z"),
                   as.POSIXct("2000-04-01 01:00:00", tz = "UTC"))
      expect_equal(parse_nc_time(1, "days since 2000-01-01"), as.POSIXct("2000-01-02", tz = "UTC"))
      expect_error(parse_nc_time(1, "fortnights since 2000-01-01"), "unrecognized")
})

test_that("ERA5 download produces a monthly wind_series with correct values and times", {
      skip_if_not_installed("ncdf4")
      local_mocked_bindings(ncss_fetch = fake_ncss())
      dir <- withr::local_tempdir()
      f <- download_wind_data("era5", xlim = c(-120, -100), ylim = c(30, 40), years = 2005,
                         months = 1:2, dir = dir, quiet = TRUE)
      expect_length(f, 2)
      expect_true(all(file.exists(f)))
      expect_match(basename(f[1]), "^era5_10m_200501_")

      ws <- wind_series(f)
      expect_s4_class(ws, "wind_series")
      expect_equal(ws@n_steps, 16)
      expect_equal(as.vector(terra::ext(ws)), c(xmin = -121, xmax = -99, ymin = 29, ymax = 41))
      expect_equal(names(ws)[c(1, 9, 17)],
                   c("u 2005-01-01 00:00:00", "u 2005-02-01 00:00:00", "v 2005-01-01 00:00:00"))
      t <- as.POSIXct("2005-02-01 05:00:00", tz = "UTC")
      expect_equal(at(ws[["u 2005-02-01 05:00:00"]], -110, 36), truth_u(250, 36, t), tolerance = 0.006)
      expect_equal(at(ws[["v 2005-02-01 05:00:00"]], -110, 36), truth_v(250, 36, t), tolerance = 0.006)
})

test_that("boxes crossing the prime meridian and the antimeridian are assembled", {
      skip_if_not_installed("ncdf4")
      log <- new.env()
      local_mocked_bindings(ncss_fetch = fake_ncss(log))
      dir <- withr::local_tempdir()
      t <- as.POSIXct("2005-01-01 03:00:00", tz = "UTC")

      ws <- wind_series(download_wind_data("era5", xlim = c(-10, 10), ylim = c(40, 50), years = 2005,
                                           months = 1, dir = dir, quiet = TRUE))
      expect_length(log$urls, 4) # two pieces x two variables
      expect_equal(as.vector(terra::ext(ws))[1:2], c(xmin = -11, xmax = 11))
      expect_equal(terra::ncol(ws), 11)
      for(lon in c(-8, 0, 8))
            expect_equal(at(ws[["u 2005-01-01 03:00:00"]], lon, 44), truth_u(lon %% 360, 44, t), tolerance = 0.006)

      ws <- wind_series(download_wind_data("era5", xlim = c(170, 200), ylim = c(40, 50), years = 2005,
                                           months = 1, dir = dir, quiet = TRUE))
      expect_equal(as.vector(terra::ext(ws))[1:2], c(xmin = 169, xmax = 201))
      expect_equal(at(ws[["v 2005-01-01 03:00:00"]], 190, 44), truth_v(190, 44, t), tolerance = 0.006)
})

test_that("CFSR and CFSv2 downloads read their layout, levels, and time units", {
      skip_if_not_installed("ncdf4")
      log <- new.env()
      local_mocked_bindings(ncss_fetch = fake_ncss(log))
      dir <- withr::local_tempdir()
      f <- download_wind_data("cfsr", xlim = c(250, 270), ylim = c(30, 40), years = 2000, months = 4,
                         dir = dir, quiet = TRUE)
      expect_match(log$urls[1], "files/g/d093001/2000/wnd10m.gdas.200004.grb2")
      ws <- wind_series(f)
      expect_equal(names(ws)[1], "u 2000-04-01 01:00:00") # CFSR files start at hour 1
      expect_equal(as.vector(terra::ext(ws))[1:2], c(xmin = 249, xmax = 271))
      t <- as.POSIXct("2000-04-01 08:00:00", tz = "UTC")
      expect_equal(at(ws[["u 2000-04-01 08:00:00"]], 260, 34), truth_u(260, 34, t), tolerance = 0.006)
      expect_equal(at(ws[["v 2000-04-01 08:00:00"]], 260, 34), truth_v(260, 34, t), tolerance = 0.006)

      download_wind_data("cfsv2", level = "850hPa", xlim = c(250, 270), ylim = c(30, 40), years = 2015,
                    months = 4, dir = dir, quiet = TRUE)
      expect_match(utils::tail(log$urls, 1), "files/g/d094001/2015/wnd850.cdas1.201504.grb2")
      expect_match(utils::tail(log$urls, 1), "var=u-component_of_wind_isobaric")
})

test_that("days and hours select time steps on the server", {
      skip_if_not_installed("ncdf4")
      log <- new.env()
      local_mocked_bindings(ncss_fetch = fake_ncss(log, n_hours = 72)) # three days per month
      dir <- withr::local_tempdir()
      args <- list(source = "era5", xlim = c(-120, -100), ylim = c(30, 40), years = 2005,
                   months = 1, dir = dir, quiet = TRUE)

      # regular hours on one day: a single strided request per variable, after the time probe
      f <- do.call(download_wind_data, c(args, list(days = 2, hours = c(0, 6, 12, 18))))
      expect_match(basename(f), "_d2_h0\\.6\\.12\\.18\\.tif$")
      expect_length(log$urls, 3)
      expect_match(log$urls[2], "time_start=2005-01-02T00:00:00Z&time_end=2005-01-02T18:00:00Z&timeStride=6")
      ws <- wind_series(f)
      expect_equal(names(ws)[1:4], paste("u 2005-01-02", c("00:00:00", "06:00:00", "12:00:00", "18:00:00")))
      t <- as.POSIXct("2005-01-02 12:00:00", tz = "UTC")
      expect_equal(at(ws[["v 2005-01-02 12:00:00"]], -110, 36), truth_v(250, 36, t), tolerance = 0.006)

      # a single hour
      ws <- wind_series(do.call(download_wind_data, c(args, list(days = 3, hours = 23))))
      expect_equal(names(ws), c("u 2005-01-03 23:00:00", "v 2005-01-03 23:00:00"))

      # irregular hours (evenings, across midnight)
      ws <- wind_series(do.call(download_wind_data, c(args, list(hours = c(20:23, 0)))))
      tm <- wind_times(ws)
      expect_equal(ws@n_steps, 15)
      expect_true(all(as.integer(format(tm, "%H", tz = "UTC")) %in% c(20:23, 0)))
      expect_equal(at(ws[["u 2005-01-02 21:00:00"]], -110, 36),
                   truth_u(250, 36, as.POSIXct("2005-01-02 21:00:00", tz = "UTC")), tolerance = 0.006)

      # whole months and selections are cached separately
      full <- do.call(download_wind_data, args)
      expect_equal(wind_series(full)@n_steps, 72)
      expect_length(list.files(dir, "^era5_10m_200501_"), 4)

      expect_error(do.call(download_wind_data, c(args, list(days = 4))), "no time steps")
})

test_that("CFSR days and hours follow the calendar, across monthly files", {
      skip_if_not_installed("ncdf4")
      log <- new.env()
      local_mocked_bindings(ncss_fetch = fake_ncss(log, full_months = TRUE))
      dir <- withr::local_tempdir()
      args <- list(source = "cfsr", xlim = c(250, 270), ylim = c(30, 40), dir = dir, quiet = TRUE)
      hrs <- function(f) format(wind_times(wind_series(f)), "%Y-%m-%d %H", tz = "UTC")

      # 1 January: 00:00 from December's file, 01:00-23:00 from January's; not 1 February 00:00
      f <- do.call(download_wind_data, c(args, list(years = 2000, months = 1, days = 1)))
      expect_equal(hrs(f), sprintf("2000-01-01 %02d", 0:23))
      expect_true(any(grepl("199912.grb2", log$urls)))
      ws <- wind_series(f)
      t <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
      expect_equal(at(ws[["u 2000-01-01 00:00:00"]], 260, 34), truth_u(260, 34, t), tolerance = 0.006)

      # midnight on the 1st only
      f <- do.call(download_wind_data, c(args, list(years = 2000, months = 3, days = 1, hours = 0)))
      expect_equal(hrs(f), "2000-03-01 00")

      # last day of a month, without the next month's first step
      f <- do.call(download_wind_data, c(args, list(years = 2000, months = 2, days = 29)))
      expect_equal(hrs(f), sprintf("2000-02-29 %02d", 0:23))

      # the first month of the data set has no earlier file
      f <- do.call(download_wind_data, c(args, list(years = 1979, months = 1, days = 1)))
      expect_equal(hrs(f), sprintf("1979-01-01 %02d", 1:23))

      # whole months keep the file as is
      f <- do.call(download_wind_data, c(args, list(years = 2000, months = 4)))
      expect_equal(range(hrs(f)), c("2000-04-01 01", "2000-05-01 00"))
})

test_that("large requests are split by time to fit the server's size limit", {
      skip_if_not_installed("ncdf4")
      # the fake grid has 11 x 6 cells in this box; allow requests of up to 3 time steps
      cap <- 11 * 6 * 4 * 3
      dir <- withr::local_tempdir()
      args <- list(source = "era5", xlim = c(-120, -100), ylim = c(30, 40), years = 2005,
                   months = 1, dir = dir, quiet = TRUE)
      local_mocked_bindings(ncss_fetch = fake_ncss(max_bytes = cap))
      expect_error(do.call(download_wind_data, args), "Request URL")

      log <- new.env()
      local_mocked_bindings(ncss_fetch = fake_ncss(log, max_bytes = cap))
      step_bytes <- 81 * 41 * 4 # the package's estimate, from ERA5's 0.25 degree grid
      withr::local_options(windscape.ncss_max_bytes = 3 * step_bytes)
      f <- do.call(download_wind_data, args)
      expect_length(log$urls, 1 + 2 * 3) # time probe, then 3 chunks per variable
      expect_true(all(grepl("time_start", log$urls[-1])))
      ws <- wind_series(f)
      expect_equal(ws@n_steps, 8)
      expect_equal(names(ws)[1:8], paste("u 2005-01-01", sprintf("%02d:00:00", 0:7)))
      t <- as.POSIXct("2005-01-01 07:00:00", tz = "UTC")
      expect_equal(at(ws[["u 2005-01-01 07:00:00"]], -110, 36), truth_u(250, 36, t), tolerance = 0.006)
})

test_that("selections are spelled out in file names, or hashed when long", {
      expect_equal(selection_tag(NULL, NULL), "")
      expect_equal(selection_tag(28, 23), "_d28_h23")
      expect_equal(selection_tag(1:15, NULL), "_d1-15")
      expect_equal(selection_tag(NULL, c(0, 20:23)), "_h0.20-23")
      long <- selection_tag(c(1, 3, 5, 7, 9, 11), c(1, 3, 5, 7, 9, 11))
      expect_match(long, "^_s[0-9a-f]{8}$")
      expect_false(identical(long, selection_tag(c(1, 3, 5, 7, 9, 11), c(1, 3, 5, 7, 9, 13))))
})

test_that("time steps are grouped into evenly spaced runs", {
      expect_equal(time_runs(5), list(list(idx = 5, by = 1)))
      expect_equal(time_runs(c(1, 7, 13, 19)), list(list(idx = c(1, 7, 13, 19), by = 6)))
      # nearly regular: one run, with a few extra steps dropped after download
      expect_equal(time_runs(c(1, 2, 4)), list(list(idx = 1:4, by = 1)))
      # irregular: separate runs, covering exactly the selected steps
      idx <- c(1, 21:25, 45:49, 69:72)
      runs <- time_runs(idx)
      expect_gt(length(runs), 1)
      expect_equal(sort(unlist(lapply(runs, `[[`, "idx"))), idx)
      for(r in runs) expect_true(length(r$idx) == 1 || all(diff(r$idx) == r$by))
})

test_that("downloaded months are cached unless overwrite = TRUE", {
      skip_if_not_installed("ncdf4")
      log <- new.env()
      local_mocked_bindings(ncss_fetch = fake_ncss(log))
      dir <- withr::local_tempdir()
      args <- list(source = "era5", xlim = c(-120, -100), ylim = c(30, 40), years = 2005, dir = dir)
      do.call(download_wind_data, c(args, list(months = 1, quiet = TRUE)))
      n <- length(log$urls)
      msg <- testthat::capture_messages(do.call(download_wind_data, c(args, list(months = 1:2))))
      expect_true(any(grepl("1 month\\(s\\) already downloaded", msg)))
      expect_length(log$urls, 2 * n) # only month 2 was fetched
      do.call(download_wind_data, c(args, list(months = 1:2, overwrite = TRUE, quiet = TRUE)))
      expect_length(log$urls, 4 * n)
})

test_that("download falls back to netCDF-3 and reports server errors", {
      skip_if_not_installed("ncdf4")
      log <- new.env()
      local_mocked_bindings(ncss_fetch = fake_ncss(log, reject_netcdf4 = TRUE))
      dir <- withr::local_tempdir()
      f <- download_wind_data("era5", xlim = c(-120, -100), ylim = c(30, 40), years = 2005, months = 1,
                         dir = dir, quiet = TRUE)
      expect_true(file.exists(f))
      expect_match(log$urls[2], "accept=netcdf$")

      # CFSR is requested as netCDF-3 directly
      log2 <- new.env()
      local_mocked_bindings(ncss_fetch = fake_ncss(log2, reject_netcdf4 = TRUE))
      download_wind_data("cfsr", xlim = c(250, 270), ylim = c(30, 40), years = 2000, months = 4,
                         dir = dir, quiet = TRUE)
      expect_length(log2$urls, 1)
      expect_match(log2$urls[1], "accept=netcdf$")

      local_mocked_bindings(ncss_fetch = fake_ncss(html = TRUE))
      expect_error(download_wind_data("era5", xlim = c(-120, -100), ylim = c(30, 40), years = 2005,
                                 months = 2, dir = dir, quiet = TRUE), "Request URL")
      expect_false(any(grepl("200502", list.files(dir)))) # no partial file left behind
})

test_that("invalid requests are rejected before downloading", {
      expect_error(download_wind_data("cfsr", xlim = c(0, 10), ylim = c(0, 10), years = 2015), "1979-2010")
      expect_error(download_wind_data("era5", level = "850hPa", xlim = c(0, 10), ylim = c(0, 10),
                                 years = 2000), "\"10m\", \"100m\"")
      expect_error(download_wind_data("cfsv2", level = "100m", xlim = c(0, 10), ylim = c(0, 10),
                                 years = 2015), "level")
      expect_error(download_wind_data("era5", xlim = c(0, 10), ylim = c(0, 10), years = 2000, months = 13),
                   "months")
      expect_error(download_wind_data("era5", xlim = c(0, 10), ylim = c(0, 10), years = 2000, days = 32),
                   "days")
      expect_error(download_wind_data("era5", xlim = c(0, 10), ylim = c(0, 10), years = 2000, hours = 24),
                   "hours")
      expect_error(download_wind_data("era5", xlim = c(0, 10), ylim = c(0, 10), years = 2000, hours = 1.5),
                   "hours")
})

test_that("ERA5 land layer is a land fraction on the wind grid", {
      skip_if_not_installed("ncdf4")
      local_mocked_bindings(ncss_fetch = fake_ncss())
      land <- download_land_mask("era5", xlim = c(170, 200), ylim = c(40, 50))
      expect_equal(names(land), "land")
      expect_equal(at(land, 176, 44), 0)
      expect_equal(at(land, 190, 44), 0.75)
      expect_error(download_land_mask("cfsv2", xlim = c(0, 10), ylim = c(0, 10)), "not yet available")
})

test_that("a rose from downloaded files matches one built in chunks, and combines", {
      skip_if_not_installed("ncdf4")
      local_mocked_bindings(ncss_fetch = fake_ncss())
      f <- download_wind_data("era5", xlim = c(-120, -100), ylim = c(30, 40), years = 2005, months = 1:3,
                         dir = withr::local_tempdir(), quiet = TRUE)
      whole <- wind_rose(wind_series(f), trans = 2)
      chunked <- withr::with_options(list(windscape.chunk_values = 1),
                                     wind_rose(wind_series(f), trans = 2))
      expect_s4_class(chunked, "wind_rose")
      expect_equal(chunked@n_steps, 24)
      expect_equal(terra::values(chunked), terra::values(whole), tolerance = 1e-6)
      expect_equal(chunked@trans(3), 9)
      expect_error(wind_rose(f), "use wind_series\\(\\) first")

      r1 <- wind_rose(wind_series(f[1]), trans = 1)
      r2 <- wind_rose(wind_series(f[2]), trans = 2)
      expect_error(combine_roses(r1, r2), "different `trans`")
      expect_equal(terra::values(combine_roses(list(r1, wind_rose(wind_series(f[2]), trans = 1)))),
                   terra::values(wind_rose(wind_series(f[1:2]), trans = 1)), tolerance = 1e-6)
      r1@n_steps <- NA_real_
      expect_error(combine_roses(r1, r1), "n_steps")
})

test_that("wind_series() checks its inputs", {
      expect_error(wind_series("no/such/file.tif"), "not found")
      r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4, nlyrs = 3, vals = 1)
      f <- withr::local_tempfile(fileext = ".tif")
      terra::writeRaster(r, f)
      expect_error(wind_series(f), "even number")
})
