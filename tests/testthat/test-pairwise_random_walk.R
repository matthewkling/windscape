# pairwise_random_walk() ------------------------

test_that("each origin's row equals a stream walk from that origin", {
      r <- noisy_rose()
      cells <- c(20, 75, 140)
      sites <- terra::xyFromCell(r, cells)
      area <- terra::values(terra::cellSize(methods::as(r, "SpatRaster")[[1]], unit = "km"))[cells, 1]
      for(value in c("deposition", "residence")){
            m <- pairwise_random_walk(r, sites, half_life = 24, value = value, density = FALSE)
            md <- pairwise_random_walk(r, sites, half_life = 24, value = value)
            for(i in seq_along(cells)){
                  w <- quietly(random_walk(r, point_raster(r, cells[i]), mode = "stream", half_life = 24))
                  expect_equal(m[i, ], vals(w[[value]])[cells], tolerance = 1e-10)
                  expect_equal(md[i, ], vals(w[[value]])[cells] / area, tolerance = 1e-10)
            }
      }
})

test_that("connectivity is asymmetric, following the wind", {
      th <- (1:72 - 0.5) / 72 * 2 * pi
      r <- build_rose(11, 21, function(x, y) list(u = 3 + 2 * sin(th), v = 2 * cos(th)),
                      xmin = -110, ymin = 35)
      sites <- terra::xyFromCell(r, terra::cellFromRowCol(r, 6, c(6, 12)))   # west, east
      m <- pairwise_random_walk(r, sites, half_life = 48)
      expect_gt(m[1, 2], 10 * m[2, 1])
})

test_that("results don't depend on timescale", {
      r <- noisy_rose()
      sites <- terra::xyFromCell(r, c(20, 75, 140))
      expect_equal(pairwise_random_walk(r, sites, half_life = 24),
                   pairwise_random_walk(r, sites, half_life = 24, timescale = 0.4), tolerance = 1e-10)
})

test_that("half-life handling", {
      r <- noisy_rose()
      sites <- terra::xyFromCell(r, c(20, 75))
      expect_error(pairwise_random_walk(r, sites), "half_life")
      expect_error(pairwise_random_walk(r, sites, half_life = Inf), "finite")
      m <- pairwise_random_walk(r, sites, value = "residence")          # no deposition
      expect_true(all(is.finite(m)) && all(m > 0))
})

test_that("site handling: shared cells, names, SpatVectors, chunks, validation", {
      r <- noisy_rose()
      sites <- rbind(a = terra::xyFromCell(r, 20)[1, ], b = terra::xyFromCell(r, 20)[1, ] + 0.1,
                     c = terra::xyFromCell(r, 140)[1, ])
      m <- pairwise_random_walk(r, sites, half_life = 24)
      expect_equal(dimnames(m), list(c("a", "b", "c"), c("a", "b", "c")))
      expect_equal(m["a", ], m["b", ])                   # same grid cell
      expect_equal(unname(pairwise_random_walk(r, terra::vect(unname(sites), crs = "EPSG:4326"), half_life = 24)),
                   unname(m))
      expect_equal(pairwise_random_walk(r, sites, half_life = 24, chunk = 1), m)
      expect_error(pairwise_random_walk(methods::as(r, "SpatRaster"), sites, half_life = 24), "wind_rose")
      expect_error(pairwise_random_walk(r, cbind(0, 0), half_life = 24), "outside")
      expect_error(pairwise_random_walk(r, cbind(1, 2, 3), half_life = 24), "two-column")
      r2 <- r
      r2[20] <- NA
      expect_error(pairwise_random_walk(r2, sites, half_life = 24), "NA")
})
