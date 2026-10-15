# stat_wind_rose() / geom_wind_rose() ------------------------

rose_plot <- function(rose = windscape_example("wind_rose"), ...){
      ggplot2::ggplot(rose, ggplot2::aes(x, y)) + geom_wind_rose(...) + ggplot2::coord_quickmap()
}
# one row per glyph from a layer's data
glyph_rows <- function(ld) ld[!duplicated(ld$group), ]
# the eight ray-tip vertices of one glyph (N first, clockwise): every tenth outline vertex
ray_tips <- function(v) v[seq(1, nrow(v), by = nrow(v) / 8), ]
# ray lengths (km) of one glyph's vertices
vertex_km <- function(v) sqrt(((v$x - v$x0) * 111.32 * cos(v$y0 * pi / 180))^2 + ((v$y - v$y0) * 110.57)^2)

test_that("with one glyph per cell, rays are conductance times neighbor distance", {
      r <- noisy_rose(nr = 4, nc = 5)
      ld <- ggplot2::layer_data(rose_plot(r, res = 100))
      g <- glyph_rows(ld)
      expect_equal(nrow(g), terra::ncell(r))
      expect_true(all(table(ld$group) == 80))
      # ray lengths in km, recovered from the ray-tip vertices, equal flow x the size factor;
      # check the ratios between rays of one glyph against conductance x distance
      i <- which.min(abs(g$x0 - terra::xFromCell(r, 7)) + abs(g$y0 - terra::yFromCell(r, 7)))
      v <- ray_tips(ld[ld$group == g$group[i], ])
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

test_that("the outline is smooth and stays within adjacent ray lengths", {
      r <- noisy_rose(nr = 4, nc = 5)
      ld <- ggplot2::layer_data(rose_plot(r, res = 100))
      for(grp in unique(ld$group)[1:5]){
            v <- ld[ld$group == grp, ]
            km <- vertex_km(v)
            tips <- km[seq(1, 80, by = 10)]
            lo <- pmin(tips, tips[c(2:8, 1)])
            hi <- pmax(tips, tips[c(2:8, 1)])
            gap <- rep(1:8, each = 10)
            expect_true(all(km >= lo[gap] - 1e-6 & km <= hi[gap] + 1e-6))
      }
})

test_that("rose_outline() passes through the ray tips and keeps a constant radius constant", {
      b <- c(0, 30, 90, 150, 180, 210, 270, 330)
      o <- rose_outline(b, rep(2, 8))
      expect_equal(sqrt(o$x^2 + o$y^2), rep(2, 80))
      r <- c(1, 5, 0, 0.5, 3, 0, 0, 2)
      o <- rose_outline(b, r)
      tips <- o[seq(1, 80, by = 10), ]
      expect_equal(sqrt(tips$x^2 + tips$y^2), r)
      expect_equal((atan2(tips$x, tips$y) * 180 / pi)[r > 0] %% 360, b[r > 0])
})

test_that("hat widths sum to 360 and are equal near the equator", {
      w <- rose_hat_widths(c(0, 45, 75), 0.5)
      expect_equal(unname(rowSums(w)), rep(360, 3))
      expect_equal(unname(w[1, ]), rep(45, 8), tolerance = 0.01)
      expect_lt(w[3, "N"], 20)  # at 75 N the N neighbor sits close to NE and NW
      expect_gt(w[3, "E"], 70)
})

test_that("normalize = TRUE draws an isotropic wind as a circle at high latitude", {
      set.seed(1)
      th <- seq(0, 2 * pi, length.out = 3601)[-1]
      iso <- build_rose(3, 3, function(x, y) list(u = 5 * sin(th), v = 5 * cos(th)),
                        ymin = 70, res = 0.5, n_steps = length(th))
      tip_km <- function(normalize){
            ld <- ggplot2::layer_data(rose_plot(iso, res = 100, normalize = normalize))
            v <- ld[ld$group == ld$group[1], ]
            vertex_km(v)
      }
      a <- tip_km(TRUE)
      expect_lt(diff(range(a)) / mean(a), 0.01)
      b <- tip_km(FALSE)
      expect_gt(diff(range(b)) / mean(b), 0.5)  # the peanut
      # computed variables describe unnormalized flow either way
      g1 <- glyph_rows(ggplot2::layer_data(rose_plot(iso, normalize = TRUE)))
      g0 <- glyph_rows(ggplot2::layer_data(rose_plot(iso, normalize = FALSE)))
      expect_equal(g1$speed, g0$speed)
      expect_equal(g1$net, g0$net)
      expect_error(rose_plot(normalize = NA), "normalize")
})

test_that("smooth = FALSE draws a star through the same ray tips", {
      r <- noisy_rose(nr = 4, nc = 5)
      a <- ggplot2::layer_data(rose_plot(r, res = 100))
      b <- ggplot2::layer_data(rose_plot(r, res = 100, smooth = FALSE))
      expect_true(all(table(b$group) == 8))
      for(grp in unique(a$group)[1:3]){
            tips <- ray_tips(a[a$group == grp, ])
            star <- b[b$group == grp, ]
            expect_equal(star$x, tips$x, ignore_attr = TRUE)
            expect_equal(star$y, tips$y, ignore_attr = TRUE)
      }
      expect_error(rose_plot(smooth = "yes"), "smooth")
})
