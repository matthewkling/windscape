# fortify methods, block aggregation, bearing scales ------------------------

test_that("fortify converts roses and fields to tidy data frames", {
      rose <- windscape_example("wind_rose")
      d <- ggplot2::fortify(rose)
      expect_named(d, c("x", "y", "SW", "W", "NW", "N", "NE", "E", "SE", "S", "total",
                        "speed", "net", "bearing", "consistency"))
      expect_equal(nrow(d), terra::ncell(rose))
      expect_equal(d$E, terra::values(rose)[, "E"])
      expect_equal(d$total, rowSums(terra::values(rose)))
      f <- ggplot2::fortify(windscape_example("wind_field"))
      expect_named(f, c("x", "y", "u", "v", "speed", "bearing"))
      expect_equal(f$speed, sqrt(f$u^2 + f$v^2))
      expect_equal(f$bearing, (atan2(f$u, f$v) * 180 / pi) %% 360)
})

test_that("fortify converts a wind_series to long format with times", {
      ws <- windscape_example("wind_series")
      d <- ggplot2::fortify(ws)
      expect_named(d, c("x", "y", "step", "time", "u", "v"))
      expect_equal(nrow(d), terra::ncell(ws) * ws@n_steps)
      expect_equal(d$u[d$step == 3], terra::values(ws[[3]])[, 1])
      expect_equal(d$v[d$step == 3], terra::values(ws[[ws@n_steps + 3]])[, 1])
      expect_s3_class(d$time, "POSIXct")
      expect_false(anyNA(d$time))
      expect_equal(format(d$time[d$step == 2][1], "%H"), "06")
})

test_that("random_walk results have a class and fortify by mode", {
      r <- noisy_rose()
      p <- quietly(random_walk(r, cbind(-90, 35), iter = 4, record = c(2, 4), half_life = 24))
      expect_s3_class(p, "random_walk")
      d <- ggplot2::fortify(p)
      expect_named(d, c("x", "y", "iteration", "hours", "airborne", "deposition"))
      expect_equal(sort(unique(d$iteration)), c(2, 4))
      expect_equal(unique(d$hours), c(2, 4) * iter_length(p))
      expect_equal(d$airborne[d$iteration == 4], terra::values(p$airborne)[, "iter4"])
      s <- quietly(random_walk(r, cbind(-90, 35), mode = "stream", half_life = 24))
      expect_named(ggplot2::fortify(s), c("x", "y", "residence", "deposition"))
      expect_equal(nrow(ggplot2::fortify(s)), terra::ncell(r))
})

test_that("ggplot accepts windscape objects directly", {
      rose <- windscape_example("wind_rose")
      p <- ggplot2::ggplot(rose, ggplot2::aes(x, y, fill = N)) + ggplot2::geom_raster()
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

test_that("block_means matches raster aggregation of the same blocks", {
      rose <- windscape_example("wind_rose")
      d <- ggplot2::fortify(rose)
      b <- block_means(d, names(rose), res = 15)
      blk <- attr(b, "block")
      a <- terra::aggregate(as(rose, "SpatRaster"), c(blk[["fy"]], blk[["fx"]]), mean, na.rm = TRUE)
      expect_equal(nrow(b), terra::ncell(a))
      # the block containing each centroid has the same mean as the raster aggregate
      cells <- terra::cellFromXY(a, as.matrix(b[, c("x", "y")]))
      expect_equal(as.matrix(b[, names(rose)]), terra::values(a)[cells, ], ignore_attr = TRUE)
      expect_equal(sum(b$n), nrow(d))
})

test_that("block_means: square blocks, NA handling, and res", {
      rose <- windscape_example("wind_rose")
      d <- ggplot2::fortify(rose)
      blk <- attr(block_means(d, "N", res = 15), "block")
      w_km <- blk[["fx"]] * blk[["dx"]] * 111.32 * cos(40 * pi / 180)
      h_km <- blk[["fy"]] * blk[["dy"]] * 110.57
      expect_lt(abs(w_km / h_km - 1), 0.2)                         # roughly square in km
      expect_gt(nrow(block_means(d, "N", 25)), nrow(block_means(d, "N", 15)))
      d$N[1:50] <- NA
      expect_equal(sum(block_means(d, "N", 15)$n), nrow(d) - 50)
      expect_error(block_means(d, "N", 0), "res")
})

test_that("bearing scales are cyclic and labeled by compass direction", {
      d <- data.frame(x = 1:5, bearing = c(0, 90, 180, 270, 360))
      p <- ggplot2::ggplot(d, ggplot2::aes(x, 1, fill = bearing)) + ggplot2::geom_tile() +
            scale_fill_bearing()
      fills <- ggplot2::layer_data(p)$fill
      expect_equal(fills[1], fills[5])                 # 0 and 360 match
      expect_equal(length(unique(fills[1:4])), 4)
      expect_equal(fills[1:4], grDevices::hcl(c(0, 90, 180, 270), 90, 65))
      d2 <- data.frame(x = 1:5, bearing = c(-90, 90, 180, 270, 450)) # wrapped
      p2 <- ggplot2::ggplot(d2, ggplot2::aes(x, 1, fill = bearing)) + ggplot2::geom_tile() +
            scale_fill_bearing()
      fills2 <- ggplot2::layer_data(p2)$fill
      expect_equal(fills2[1], fills[4])
      expect_equal(fills2[5], fills[2])
      pc <- ggplot2::ggplot(d2, ggplot2::aes(x, 1, colour = bearing)) + ggplot2::geom_point() +
            scale_colour_bearing()
      expect_equal(length(unique(ggplot2::layer_data(pc)$colour)), 3) # -90 = 270 and 450 = 90
})

test_that("bearing scale legend shows eight directions, and breaks can be overridden", {
      d <- data.frame(x = 1:8, bearing = seq(0, 315, 45))
      p <- ggplot2::ggplot(d, ggplot2::aes(x, 1, fill = bearing)) + ggplot2::geom_tile()
      sc <- ggplot2::ggplot_build(p + scale_fill_bearing())$plot$scales$get_scales("fill")
      expect_equal(sc$get_labels(), c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))
      sc4 <- ggplot2::ggplot_build(p + scale_fill_bearing(breaks = c(0, 90, 180, 270),
                                                          labels = c("N", "E", "S", "W")))$plot$scales$get_scales("fill")
      expect_equal(sc4$get_labels(), c("N", "E", "S", "W"))
})

test_that("rose fortify flow summaries match the glyph stat at one glyph per cell", {
      rose <- windscape_example("wind_rose")
      d <- ggplot2::fortify(rose)
      p <- ggplot2::ggplot(rose, ggplot2::aes(x, y)) + geom_wind_rose(res = 1000)
      g <- ggplot2::layer_data(p)
      g <- g[!duplicated(g$group), ]
      m <- match(paste(round(d$x, 6), round(d$y, 6)), paste(round(g$x0, 6), round(g$y0, 6)))
      expect_false(anyNA(m))
      for(v in c("speed", "net", "bearing", "consistency")) expect_equal(d[[v]], g[[v]][m], tolerance = 1e-10)
      expect_true(all(d$consistency >= 0 & d$consistency <= 1))
})

test_that("rose fortify speed is the mean wind speed in km/h", {
      ws <- windscape_example("wind_series")
      d <- ggplot2::fortify(wind_rose(ws, trans = 1))
      n <- ws@n_steps
      v <- terra::values(ws)
      expect_equal(d$speed, rowMeans(sqrt(v[, 1:n]^2 + v[, n + 1:n]^2)) * 3.6, tolerance = 1e-8)
})

test_that("projected roses get total conductance but no flow summaries", {
      d <- ggplot2::fortify(uniform_rose(u = 3))
      expect_true("total" %in% names(d))
      expect_false(any(c("speed", "net", "bearing", "consistency") %in% names(d)))
})

test_that("fortified variables work as background layers", {
      f <- windscape_example("wind_field")
      p <- ggplot2::ggplot(f, ggplot2::aes(x, y)) + ggplot2::geom_raster(ggplot2::aes(fill = speed)) +
            geom_wind_arrow()
      expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
      rose <- windscape_example("wind_rose")
      p2 <- ggplot2::ggplot(rose, ggplot2::aes(x, y)) +
            ggplot2::geom_raster(ggplot2::aes(alpha = consistency), fill = "gray30") +
            geom_wind_rose()
      expect_s3_class(ggplot2::ggplotGrob(p2), "gtable")
})

test_that("British and American spellings of colour both work", {
      kat <- windscape_example("wind_field")
      rose <- windscape_example("wind_rose")
      a <- function(...) ggplot2::layer_data(ggplot2::ggplot(kat, ggplot2::aes(x, y)) + ...)
      # fixed colours, as layer arguments
      expect_equal(a(geom_wind_arrow(color = "red"))$colour, a(geom_wind_arrow(colour = "red"))$colour)
      expect_true(all(a(stat_wind_trail(color = "red"))$colour == "red"))
      r <- ggplot2::layer_data(ggplot2::ggplot(rose, ggplot2::aes(x, y)) + geom_wind_rose(color = "red"))
      expect_true(all(r$colour == "red"))
      # mapped colours, with either scale spelling
      p1 <- ggplot2::ggplot(kat, ggplot2::aes(x, y)) +
            geom_wind_arrow(ggplot2::aes(color = ggplot2::after_stat(bearing))) + scale_color_bearing()
      p2 <- ggplot2::ggplot(kat, ggplot2::aes(x, y)) +
            geom_wind_arrow(ggplot2::aes(colour = ggplot2::after_stat(bearing))) + scale_colour_bearing()
      expect_equal(ggplot2::layer_data(p1)$colour, ggplot2::layer_data(p2)$colour)
      expect_equal(ggplot2::layer_data(p1)$colour, grDevices::hcl(ggplot2::layer_data(p1)$bearing, 90, 65))
})
