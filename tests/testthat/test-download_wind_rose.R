# download_wind_rose(), wind_rose_catalog(), wind_rose_cache() ---------------------------------

test_that("the shipped catalog is complete and consistent", {
      catalog <- wind_rose_catalog()
      expect_equal(as.vector(table(catalog$tier)[c("all", "month", "month_of_year", "year")]),
                   c(1, 384, 12, 32))
      m <- catalog[catalog$tier == "month", ]
      expect_equal(m$n_steps, days_in_month(m$year, m$month) * 24)
      expect_equal(catalog$n_steps[catalog$tier == "all"], sum(m$n_steps))
      expect_equal(anyDuplicated(catalog$file), 0)
      expect_true(all(nchar(catalog$md5) == 32))
      expect_true(all(startsWith(catalog$url, paste0(rose_host(), "/cfsr-roses-v1/"))))
})

test_that("rose_files picks the fewest files covering a request", {
      catalog <- wind_rose_catalog()
      expect_equal(rose_files(catalog, NULL, NULL)$tier, "all")
      expect_equal(rose_files(catalog, 1979:2010, 1:12)$tier, "all")
      f <- rose_files(catalog, NULL, c(8, 6, 7))
      expect_equal(f$tier, rep("month_of_year", 3))
      expect_equal(f$month, 6:8)
      f <- rose_files(catalog, 1990:1991, NULL)
      expect_equal(f$tier, rep("year", 2))
      f <- rose_files(catalog, 1990:1991, 6:8)
      expect_equal(f$tier, rep("month", 6))
      expect_equal(sum(f$n_steps), 2 * (30 + 31 + 31) * 24)
})

test_that("rose_files validates years and months", {
      catalog <- wind_rose_catalog()
      expect_error(rose_files(catalog, 2015, NULL), "1979-2010.*2015")
      expect_error(rose_files(catalog, 1990.5, NULL), "whole years")
      expect_error(rose_files(catalog, NULL, 13), "1 to 12")
      expect_error(rose_files(catalog, NULL, integer()), "1 to 12")
      expect_error(rose_files(catalog[catalog$tier != "all", ], NULL, NULL), "missing files")
})

test_that("download_wind_rose returns single files and hour-weighted combinations", {
      h <- fake_rose_host()
      r <- download_wind_rose(year = 2001, month = 3, quiet = TRUE)
      expect_s4_class(r, "wind_rose")
      expect_equal(r@n_steps, 31 * 24)
      expect_equal(r@trans(3), 3)
      expect_equal(terra::values(r), terra::values(h$monthly[[3]]), ignore_attr = TRUE)

      check <- function(year, month){
            i <- which(h$ym$year %in% (if(is.null(year)) 2001:2002 else year) &
                             h$ym$month %in% (if(is.null(month)) 1:12 else month))
            x <- download_wind_rose(year = year, month = month, quiet = TRUE)
            expect_equal(x@n_steps, sum(h$ym$n_steps[i]))
            expect_equal(terra::values(x), terra::values(combine_roses(h$monthly[i])),
                         tolerance = 1e-6, ignore_attr = TRUE)
      }
      check(NULL, NULL)       # full-period file
      check(NULL, 6:8)        # month-of-year files
      check(2002, NULL)       # annual file
      check(2002, c(1, 12))   # monthly files
})

test_that("download_wind_rose caches whole files, and wind_rose_cache manages them", {
      h <- fake_rose_host()
      expect_equal(nrow(wind_rose_cache()), 0)
      expect_message(expect_message(download_wind_rose(year = 2001), "getting 1 wind rose file"),
                     "downloading cfsr_10m_year_2001.tif")
      cached <- wind_rose_cache()
      expect_equal(cached$file, "cfsr_10m_year_2001.tif")
      expect_equal(cached$release, h$release)

      # once cached, the host isn't consulted
      writeBin(as.raw(1:10), file.path(h$host, h$release, "cfsr_10m_year_2001.tif"))
      expect_silent(r <- download_wind_rose(year = 2001, quiet = TRUE))
      expect_equal(r@n_steps, 365 * 24)

      expect_message(wind_rose_cache(clear = TRUE), "deleted 1 file")
      expect_equal(nrow(wind_rose_cache()), 0)

      # a corrupt download fails its checksum, leaving nothing behind
      expect_error(download_wind_rose(year = 2001, quiet = TRUE), "checksum")
      expect_equal(nrow(wind_rose_cache()), 0)
})

test_that("download_wind_rose crops to ext, reading remotely or from the cache", {
      h <- fake_rose_host()
      e <- c(-50, 20, -10, 35)
      ref <- terra::crop(h$monthly[[1]], terra::ext(e), snap = "out")

      rem <- suppressMessages(download_wind_rose(year = 2001, month = 1, ext = e))
      expect_equal(nrow(wind_rose_cache()), 0)
      expect_true(terra::compareGeom(rem, ref))
      expect_equal(terra::values(rem), terra::values(ref), ignore_attr = TRUE)
      expect_equal(rem@n_steps, 31 * 24)

      loc <- download_wind_rose(year = 2001, month = 1, ext = e, cache = TRUE, quiet = TRUE)
      expect_equal(nrow(wind_rose_cache()), 1)
      expect_equal(terra::values(loc), terra::values(ref), ignore_attr = TRUE)

      # extent of another object
      x <- download_wind_rose(year = 2001, month = 1, ext = ref, cache = TRUE, quiet = TRUE)
      expect_true(terra::compareGeom(x, ref))
})

test_that("download_wind_rose validates its inputs", {
      fake_rose_host()
      expect_error(download_wind_rose(source = "era5"), "no pre-built roses")
      expect_error(download_wind_rose(ext = c(170, 200, 0, 10)), "-180 to 180")
      expect_error(download_wind_rose(ext = 1:3), "c\\(xmin, xmax, ymin, ymax\\)")
      expect_error(download_wind_rose(year = 1990), "available for 2001-2002")
})

test_that("download_wind_rose reads hosted roses", {
      skip_on_cran()
      skip_if_offline("github.com")
      url <- wind_rose_catalog()$url[1]
      ok <- tryCatch(attr(curlGetHeaders(url), "status") < 400, error = function(e) FALSE)
      skip_if_not(ok, "hosted roses not reachable")
      r <- download_wind_rose(year = 1990, month = 1, ext = c(-125, -120, 40, 45), quiet = TRUE)
      expect_s4_class(r, "wind_rose")
      expect_equal(r@n_steps, 744)
      expect_false(anyNA(terra::values(r)))
})


# Reading windscape metadata in wind_rose() ------------------------------------------------------

test_that("wind_rose reads n_steps and trans recorded in a file's metadata", {
      skip_if_not_installed("sf")
      plain <- withr::local_tempfile(fileext = ".tif")
      terra::writeRaster(windscape_example("wind_rose"), plain)
      tagged <- function(...){
            f <- withr::local_tempfile(fileext = ".tif", .local_envir = parent.frame(2))
            mo <- as.vector(rbind("-mo", c(...)))
            sf::gdal_utils("translate", plain, f, options = mo)
            f
      }
      f <- tagged("windscape_rose_format=1", "windscape_n_steps=576", "windscape_trans=2")
      x <- wind_rose(f)
      expect_equal(x@n_steps, 576)
      expect_equal(x@trans(3), 9)
      expect_silent(wind_rose(f, trans = 2, n_steps = 576))
      expect_warning(y <- wind_rose(f, trans = 1), "differs")
      expect_equal(y@trans(3), 3)
      expect_warning(wind_rose(f, n_steps = 10), "differs")

      z <- wind_rose(plain)
      expect_true(is.na(z@n_steps))
      expect_equal(z@trans(3), 3)

      expect_error(wind_rose(tagged("windscape_rose_format=99")), "newer format")
})
