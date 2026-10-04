test_that("net_flow returns a wind_field of net flow on the rose's grid", {
      r <- noisy_rose()
      f <- net_flow(r)
      expect_s4_class(f, "wind_field")
      expect_equal(names(f), c("u", "v"))
      expect_true(terra::compareGeom(f, as(r, "SpatRaster")))
})

test_that("net_flow matches the net flow statistics from fortify()", {
      r <- noisy_rose()
      f <- terra::values(net_flow(r))
      d <- ggplot2::fortify(r, na.rm = FALSE)
      expect_equal(sqrt(f[, "u"]^2 + f[, "v"]^2), d$net, tolerance = 1e-12)
      expect_equal((atan2(f[, "u"], f[, "v"]) * 180 / pi) %% 360, d$bearing, tolerance = 1e-10)
})

test_that("steady wind toward a neighbor gives net flow equal to wind speed, in km/h", {
      for(lat in c(0, 30, 60)){
            east <- build_rose(3, 3, function(x, y) list(u = rep(5, 10), v = rep(0, 10)), ymin = lat)
            fe <- terra::values(net_flow(east))
            expect_equal(fe[, "u"], rep(5 * 3.6, 9), tolerance = 1e-10)
            expect_equal(fe[, "v"], rep(0, 9), tolerance = 1e-10)
            north <- build_rose(3, 3, function(x, y) list(u = rep(0, 10), v = rep(5, 10)), ymin = lat)
            fn <- terra::values(net_flow(north))
            expect_equal(fn[, "u"], rep(0, 9), tolerance = 1e-10)
            expect_equal(fn[, "v"], rep(5 * 3.6, 9), tolerance = 1e-10)
      }
})

test_that("opposing winds cancel", {
      r <- build_rose(3, 3, function(x, y) list(u = c(4, -4), v = c(0, 0)))
      f <- terra::values(net_flow(r))
      expect_equal(f[, "u"], rep(0, 9), tolerance = 1e-10)
})

test_that("net_flow keeps NA cells and rejects bad input", {
      r <- noisy_rose()
      r[5] <- NA
      f <- net_flow(r)
      expect_true(all(is.na(terra::values(f)[5, ])))
      expect_false(anyNA(terra::values(f)[-5, ]))
      expect_error(net_flow(as(r, "SpatRaster")), "must be a wind_rose")
      expect_error(net_flow(uniform_rose()), "longitude/latitude")
})

test_that("net flow is the drift velocity of a random walk on the rose", {
      r <- noisy_rose()
      res <- terra::res(r)
      xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
      f <- terra::values(net_flow(r))
      for(cell in c(50, 95)){
            d <- rw_neighbor_displacements(terra::yFromCell(r, cell), mean(res))
            k <- match(paste(round((xy[, 1] - xy[cell, 1]) / res[1]), round((xy[, 2] - xy[cell, 2]) / res[2])),
                       paste(c(-1, -1, -1, 0, 1, 1, 1, 0), c(-1, 0, 1, 1, 1, 0, -1, -1)))
            nb <- !is.na(k)
            for(lc in c(TRUE, FALSE)){
                  w <- quietly(random_walk(r, point_raster(r, cell), iter = 1, density = FALSE,
                                           latitude_correction = lc))
                  a <- vals(w$airborne)
                  t <- rw_max_step(if(lc) rw_latitude_correction(r) else r)
                  drift <- c(sum(a[nb] * d[k[nb], 1]), sum(a[nb] * d[k[nb], 2])) / sum(a) / t
                  expect_equal(drift, unname(f[cell, ]), tolerance = 1e-10)
            }
      }
})
