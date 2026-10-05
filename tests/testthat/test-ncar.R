# ncar_download(), read_wind_series(), ncar_land(), combine_roses() --------------------------
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
      f <- ncar_download("era5", xlim = c(-120, -100), ylim = c(30, 40), years = 2005,
                         months = 1:2, dir = dir, quiet = TRUE)
      expect_length(f, 2)
      expect_true(all(file.exists(f)))
      expect_match(basename(f[1]), "^era5_10m_200501_")

      ws <- read_wind_series(f)
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

      ws <- read_wind_series(ncar_download("era5", xlim = c(-10, 10), ylim = c(40, 50), years = 2005,
                                           months = 1, dir = dir, quiet = TRUE))
      expect_length(log$urls, 4) # two pieces x two variables
      expect_equal(as.vector(terra::ext(ws))[1:2], c(xmin = -11, xmax = 11))
      expect_equal(terra::ncol(ws), 11)
      for(lon in c(-8, 0, 8))
            expect_equal(at(ws[["u 2005-01-01 03:00:00"]], lon, 44), truth_u(lon %% 360, 44, t), tolerance = 0.006)

      ws <- read_wind_series(ncar_download("era5", xlim = c(170, 200), ylim = c(40, 50), years = 2005,
                                           months = 1, dir = dir, quiet = TRUE))
      expect_equal(as.vector(terra::ext(ws))[1:2], c(xmin = 169, xmax = 201))
      expect_equal(at(ws[["v 2005-01-01 03:00:00"]], 190, 44), truth_v(190, 44, t), tolerance = 0.006)
})

test_that("CFSR and CFSv2 downloads read their layout, levels, and time units", {
      skip_if_not_installed("ncdf4")
      log <- new.env()
      local_mocked_bindings(ncss_fetch = fake_ncss(log))
      dir <- withr::local_tempdir()
      f <- ncar_download("cfsr", xlim = c(250, 270), ylim = c(30, 40), years = 2000, months = 4,
                         dir = dir, quiet = TRUE)
      expect_match(log$urls[1], "files/g/d093001/2000/wnd10m.gdas.200004.grb2")
      ws <- read_wind_series(f)
      expect_equal(names(ws)[1], "u 2000-04-01 01:00:00") # CFSR files start at hour 1
      expect_equal(as.vector(terra::ext(ws))[1:2], c(xmin = 249, xmax = 271))
      t <- as.POSIXct("2000-04-01 08:00:00", tz = "UTC")
      expect_equal(at(ws[["u 2000-04-01 08:00:00"]], 260, 34), truth_u(260, 34, t), tolerance = 0.006)
      expect_equal(at(ws[["v 2000-04-01 08:00:00"]], 260, 34), truth_v(260, 34, t), tolerance = 0.006)

      ncar_download("cfsv2", level = "850hPa", xlim = c(250, 270), ylim = c(30, 40), years = 2015,
                    months = 4, dir = dir, quiet = TRUE)
      expect_match(utils::tail(log$urls, 1), "files/g/d094001/2015/wnd850.cdas1.201504.grb2")
      expect_match(utils::tail(log$urls, 1), "var=u-component_of_wind_isobaric")
})

test_that("time_stride thins time steps on the server", {
      skip_if_not_installed("ncdf4")
      log <- new.env()
      local_mocked_bindings(ncss_fetch = fake_ncss(log))
      f <- ncar_download("era5", xlim = c(-120, -100), ylim = c(30, 40), years = 2005, months = 1,
                         time_stride = 3, dir = withr::local_tempdir(), quiet = TRUE)
      expect_match(log$urls[1], "timeStride=3")
      expect_match(basename(f), "_t3\\.tif$")
      ws <- read_wind_series(f)
      expect_equal(ws@n_steps, 3)
      expect_equal(names(ws)[1:3], paste("u 2005-01-01", c("00:00:00", "03:00:00", "06:00:00")))
})

test_that("downloaded months are cached unless overwrite = TRUE", {
      skip_if_not_installed("ncdf4")
      log <- new.env()
      local_mocked_bindings(ncss_fetch = fake_ncss(log))
      dir <- withr::local_tempdir()
      args <- list(source = "era5", xlim = c(-120, -100), ylim = c(30, 40), years = 2005, dir = dir)
      do.call(ncar_download, c(args, list(months = 1, quiet = TRUE)))
      n <- length(log$urls)
      msg <- testthat::capture_messages(do.call(ncar_download, c(args, list(months = 1:2))))
      expect_true(any(grepl("1 month\\(s\\) already downloaded", msg)))
      expect_length(log$urls, 2 * n) # only month 2 was fetched
      do.call(ncar_download, c(args, list(months = 1:2, overwrite = TRUE, quiet = TRUE)))
      expect_length(log$urls, 4 * n)
})

test_that("download falls back to netCDF-3 and reports server errors", {
      skip_if_not_installed("ncdf4")
      log <- new.env()
      local_mocked_bindings(ncss_fetch = fake_ncss(log, reject_netcdf4 = TRUE))
      dir <- withr::local_tempdir()
      f <- ncar_download("era5", xlim = c(-120, -100), ylim = c(30, 40), years = 2005, months = 1,
                         dir = dir, quiet = TRUE)
      expect_true(file.exists(f))
      expect_match(log$urls[2], "accept=netcdf$")

      local_mocked_bindings(ncss_fetch = fake_ncss(html = TRUE))
      expect_error(ncar_download("era5", xlim = c(-120, -100), ylim = c(30, 40), years = 2005,
                                 months = 2, dir = dir, quiet = TRUE), "Request URL")
      expect_false(any(grepl("200502", list.files(dir)))) # no partial file left behind
})

test_that("invalid requests are rejected before downloading", {
      expect_error(ncar_download("cfsr", xlim = c(0, 10), ylim = c(0, 10), years = 2015), "1979-2010")
      expect_error(ncar_download("era5", level = "850hPa", xlim = c(0, 10), ylim = c(0, 10),
                                 years = 2000), "\"10m\", \"100m\"")
      expect_error(ncar_download("cfsv2", level = "100m", xlim = c(0, 10), ylim = c(0, 10),
                                 years = 2015), "level")
      expect_error(ncar_download("era5", xlim = c(0, 10), ylim = c(0, 10), years = 2000, months = 13),
                   "months")
      expect_error(ncar_download("era5", xlim = c(0, 10), ylim = c(0, 10), years = 2000, time_stride = 0),
                   "time_stride")
})

test_that("ERA5 land layer is a land fraction on the wind grid", {
      skip_if_not_installed("ncdf4")
      local_mocked_bindings(ncss_fetch = fake_ncss())
      land <- ncar_land("era5", xlim = c(170, 200), ylim = c(40, 50))
      expect_equal(names(land), "land")
      expect_equal(at(land, 176, 44), 0)
      expect_equal(at(land, 190, 44), 0.75)
      expect_error(ncar_land("cfsv2", xlim = c(0, 10), ylim = c(0, 10)), "not yet available")
})

test_that("wind_rose() on files matches a rose built from all time steps at once", {
      skip_if_not_installed("ncdf4")
      local_mocked_bindings(ncss_fetch = fake_ncss())
      f <- ncar_download("era5", xlim = c(-120, -100), ylim = c(30, 40), years = 2005, months = 1:3,
                         dir = withr::local_tempdir(), quiet = TRUE)
      chunked <- wind_rose(f, trans = 2)
      whole <- wind_rose(read_wind_series(f), trans = 2)
      expect_s4_class(chunked, "wind_rose")
      expect_equal(chunked@n_steps, 24)
      expect_equal(terra::values(chunked), terra::values(whole), tolerance = 1e-6)
      expect_equal(chunked@trans(3), 9)
      expect_error(wind_rose(f, filename = tempfile(fileext = ".tif")), "filename")

      r1 <- wind_rose(read_wind_series(f[1]), trans = 1)
      r2 <- wind_rose(read_wind_series(f[2]), trans = 2)
      expect_error(combine_roses(r1, r2), "different `trans`")
      expect_equal(terra::values(combine_roses(list(r1, wind_rose(read_wind_series(f[2]), trans = 1)))),
                   terra::values(wind_rose(read_wind_series(f[1:2]), trans = 1)), tolerance = 1e-6)
      r1@n_steps <- NA_real_
      expect_error(combine_roses(r1, r1), "n_steps")
})

test_that("read_wind_series() checks its inputs", {
      expect_error(read_wind_series("no/such/file.tif"), "not found")
      r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4, nlyrs = 3, vals = 1)
      f <- withr::local_tempfile(fileext = ".tif")
      terra::writeRaster(r, f)
      expect_error(read_wind_series(f), "odd number")
})
