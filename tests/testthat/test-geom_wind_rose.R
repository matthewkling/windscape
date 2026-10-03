# stat_wind_rose() / geom_wind_rose() ------------------------

rose_plot <- function(rose = windscape_example("wind_rose"), ...){
      ggplot2::ggplot(rose, ggplot2::aes(x, y)) + geom_wind_rose(...) + ggplot2::coord_quickmap()
}
# one row per glyph from a layer's data
glyph_rows <- function(ld) ld[!duplicated(ld$group), ]

test_that("with one glyph per cell, rays are conductance times neighbor distance", {
      r <- noisy_rose(nr = 4, nc = 5)
      ld <- ggplot2::layer_data(rose_plot(r, res = 100))
      g <- glyph_rows(ld)
      expect_equal(nrow(g), terra::ncell(r))
      expect_true(all(table(ld$group) == 8))
      # ray lengths in km, recovered from vertex offsets, equal flow x the size factor;
      # check the ratios between rays of one glyph against conductance x distance
      i <- which.min(abs(g$x0 - terra::xFromCell(r, 7)) + abs(g$y0 - terra::yFromCell(r, 7)))
      v <- ld[ld$group == g$group[i], ]
      km <- sqrt(((v$x - v$x0) * 111.32 * cos(v$y0 * pi / 180))^2 + ((v$y - v$y0) * 110.57)^2)
      flow <- terra::values(r)[7, ] * sqrt(rowSums(rw_neighbor_displacements(terra::yFromCell(r, 7), 1)^2))
      flow <- flow[match(c("N", "NE", "E", "SE", "S", "SW", "W", "NW"), names(flow))]
      expect_equal(km / sum(km), unname(flow / sum(flow)), tolerance = 1e-3)
      expect_equal(g$speed[i], sum(flow), tolerance = 1e-8)
})

test_that("glyph speed is the mean wind speed in km/h", {
      ws <- windscape_example("wind_series")
      rose <- wind_rose(ws, trans = 1)
      n <- ws@n_steps
      v <- terra::values(ws)
      speed_kmh <- rowMeans(sqrt(v[, 1:n]^2 + v[, n + 1:n]^2)) * 3.6
      g <- glyph_rows(ggplot2::layer_data(rose_plot(rose, res = 1000))) # one glyph per cell
      cells <- terra::cellFromXY(rose, as.matrix(g[, c("x0", "y0")]))
      expect_equal(g$speed, speed_kmh[cells], tolerance = 1e-8)
})

test_that("a uniform westerly gives eastward, fully consistent glyphs", {
      g <- glyph_rows(ggplot2::layer_data(rose_plot(uv_rose(u = 5, v = 0))))
      expect_equal(g$bearing, rep(90, nrow(g)), tolerance = 1e-6)
      expect_equal(g$consistency, rep(1, nrow(g)), tolerance = 1e-6)
      expect_equal(g$speed, rep(18, nrow(g)), tolerance = 1e-6) # 5 m/s = 18 km/h
})

test_that("block averaging matches raster aggregation", {
      rose <- windscape_example("wind_rose")
      g <- glyph_rows(ggplot2::layer_data(rose_plot(rose, res = 15)))
      expect_equal(sum(g$n), terra::ncell(rose))
      expect_true(all(g$consistency >= 0 & g$consistency <= 1))
      expect_true(all(g$bearing >= 0 & g$bearing < 360))
})

test_that("default fill is bearing, with the bearing scale added", {
      p <- rose_plot()
      ld <- ggplot2::layer_data(p)
      expect_equal(ld$fill, grDevices::hcl(ld$bearing, 90, 65))
      expect_true(any(vapply(p$scales$scales, function(s) "fill" %in% s$aesthetics, logical(1))))
})

test_that("facets share blocks and glyph scale", {
      d <- ggplot2::fortify(windscape_example("wind_rose"))
      dd <- rbind(cbind(d, version = "a"), cbind(d, version = "b"))
      p <- ggplot2::ggplot(dd, ggplot2::aes(x, y)) + geom_wind_rose(res = 10) +
            ggplot2::facet_wrap(~version)
      ld <- ggplot2::layer_data(p)
      a <- ld[ld$PANEL == 1, c("x", "y", "bearing")]
      b <- ld[ld$PANEL == 2, c("x", "y", "bearing")]
      expect_equal(a, b, ignore_attr = TRUE)
})

test_that("user fill mappings and constant columns carry through", {
      d <- ggplot2::fortify(windscape_example("wind_rose"))
      dd <- rbind(cbind(d, season = "summer"), cbind(d, season = "winter"))
      dd$N[dd$season == "winter"] <- dd$N[dd$season == "winter"] * 3
      p <- ggplot2::ggplot(dd, ggplot2::aes(x, y)) +
            geom_wind_rose(ggplot2::aes(fill = season), res = 8)
      expect_false(any(vapply(p$scales$scales, function(s) "fill" %in% s$aesthetics, logical(1))))
      ld <- ggplot2::layer_data(p)
      expect_equal(length(unique(ld$fill)), 2)       # group-level fill survives aggregation
      g <- ld[!duplicated(ld$group), ]
      expect_equal(nrow(g) %% 2, 0)                   # glyphs computed separately per season
      s <- split(g$bearing, g$fill)
      expect_false(isTRUE(all.equal(s[[1]], s[[2]]))) # winter's tripled northward flow differs
})

test_that("res and scale control glyph number and size", {
      n <- function(p) length(unique(ggplot2::layer_data(p)$group))
      expect_gt(n(rose_plot(res = 25)), n(rose_plot(res = 10)))
      a <- ggplot2::layer_data(rose_plot(scale = 1))
      b <- ggplot2::layer_data(rose_plot(scale = 2))
      expect_equal(b$x - b$x0, 2 * (a$x - a$x0), tolerance = 1e-10)
      expect_equal(b$bearing, a$bearing)
      expect_error(rose_plot(scale = 0), "scale")
      expect_error(rose_plot(res = 0), "res")
})

test_that("the geom draws, with and without centers and blending", {
      for(p in list(rose_plot(), rose_plot(center = FALSE), rose_plot(saturation = 0))){
            expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
      }
})

test_that("blend_gray blends toward gray as consistency falls", {
      cols <- c("#FF0000", "#FF0000", "#FF0000", NA)
      out <- blend_gray(cols, c(1, 0.25, 0, 0.5), saturation = 0.5)
      expect_equal(out[1], "#FF0000")
      expect_equal(out[3], "#A0A0A0")
      expect_equal(out[2], grDevices::rgb(t((grDevices::col2rgb("#FF0000") + grDevices::col2rgb("#A0A0A0")) / 2),
                                          maxColorValue = 255))
      expect_true(is.na(out[4]))
      expect_equal(blend_gray(cols, c(0, 0, 0, 0), saturation = 0), cols)
})

test_that("a wind_rose can be passed as layer data", {
      rose <- windscape_example("wind_rose")
      p <- ggplot2::ggplot() + geom_wind_rose(ggplot2::aes(x, y), data = rose)
      expect_gt(nrow(ggplot2::layer_data(p)), 0)
})

test_that("x and y must be mapped", {
      p <- ggplot2::ggplot(windscape_example("wind_rose")) + geom_wind_rose()
      expect_error(ggplot2::layer_data(p), "x|y")
})
