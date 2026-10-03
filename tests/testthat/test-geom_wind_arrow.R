# stat_wind_arrow() / geom_wind_arrow() ------------------------

arrow_plot <- function(field = windscape_example("wind_field"), ...){
      ggplot2::ggplot(field, ggplot2::aes(x, y)) + geom_wind_arrow(...) + ggplot2::coord_quickmap()
}
uniform_field <- function(u, v){
      r <- terra::rast(nrows = 10, ncols = 10, xmin = -100, xmax = -90, ymin = 30, ymax = 40,
                       crs = "EPSG:4326", nlyrs = 2)
      terra::values(r) <- cbind(rep(u, 100), rep(v, 100))
      wind_field(r)
}

test_that("a uniform westerly gives eastward, fully consistent, equal arrows", {
      ld <- ggplot2::layer_data(arrow_plot(uniform_field(5, 0), res = 5))
      expect_equal(ld$bearing, rep(90, nrow(ld)), tolerance = 1e-8)
      expect_equal(ld$consistency, rep(1, nrow(ld)), tolerance = 1e-8)
      expect_equal(ld$speed, rep(5, nrow(ld)), tolerance = 1e-8)
      expect_equal(ld$y, ld$yend)                                     # horizontal
      expect_true(all(ld$xend > ld$x))                                # pointing east
      expect_equal((ld$x + ld$xend) / 2, ld$x0, tolerance = 1e-10)    # centered (pivot 0.5)
})

test_that("block means are vector means of u and v", {
      f <- windscape_example("wind_field")
      ld <- ggplot2::layer_data(arrow_plot(f, res = 20))
      d <- ggplot2::fortify(f)
      b <- block_means(d, c("u", "v"), res = 20)
      m <- match(paste(round(b$x, 6), round(b$y, 6)), paste(round(ld$x0, 6), round(ld$y0, 6)))
      expect_false(anyNA(m))
      expect_equal(ld$speed[m], sqrt(b$u^2 + b$v^2), tolerance = 1e-10)
      expect_equal(ld$bearing[m], (atan2(b$u, b$v) * 180 / pi) %% 360, tolerance = 1e-10)
      expect_true(all(ld$consistency >= 0 & ld$consistency <= 1 + 1e-12))
      expect_equal(sum(ld$n), nrow(d))
})

test_that("rotation within blocks lowers consistency", {
      f <- windscape_example("wind_field")
      ld <- ggplot2::layer_data(arrow_plot(f, res = 6))
      single <- ggplot2::layer_data(arrow_plot(f, res = 1000))    # one arrow per cell
      expect_equal(single$consistency, rep(1, nrow(single)), tolerance = 1e-8)
      expect_lt(min(ld$consistency), 0.9)                         # e.g. the block holding the eye
})

test_that("pivot, scale, and res control arrow geometry", {
      a <- ggplot2::layer_data(arrow_plot(pivot = 0))
      expect_equal(a$x, a$x0)
      expect_equal(a$y, a$y0)
      b1 <- ggplot2::layer_data(arrow_plot(scale = 1))
      b2 <- ggplot2::layer_data(arrow_plot(scale = 2))
      expect_equal(b2$xend - b2$x, 2 * (b1$xend - b1$x), tolerance = 1e-10)
      expect_gt(nrow(ggplot2::layer_data(arrow_plot(res = 30))), nrow(b1))
      expect_error(arrow_plot(scale = 0), "scale")
      expect_error(arrow_plot(res = 0), "res")
      expect_error(arrow_plot(pivot = NA), "pivot")
})

test_that("facets share blocks and arrow scale", {
      d <- ggplot2::fortify(windscape_example("wind_field"))
      dd <- rbind(cbind(d, version = "a"), cbind(d, version = "b"))
      p <- ggplot2::ggplot(dd, ggplot2::aes(x, y)) + geom_wind_arrow() + ggplot2::facet_wrap(~version)
      ld <- ggplot2::layer_data(p)
      expect_equal(ld[ld$PANEL == 1, c("x", "y", "xend", "yend")],
                   ld[ld$PANEL == 2, c("x", "y", "xend", "yend")], ignore_attr = TRUE)
})

test_that("colour can be mapped to bearing, and the geom draws in all styles", {
      p <- arrow_plot(mapping = ggplot2::aes(colour = ggplot2::after_stat(bearing),
                                             consistency = ggplot2::after_stat(consistency))) +
            scale_colour_bearing()
      ld <- ggplot2::layer_data(p)
      expect_equal(ld$colour, grDevices::hcl(ld$bearing, 90, 65))
      for(q in list(p, arrow_plot(), arrow_plot(arrow = NULL, pivot = 0, center = TRUE))){
            expect_s3_class(ggplot2::ggplotGrob(q), "gtable")
      }
})

test_that("u and v map automatically; other column names can be mapped in the layer", {
      d <- ggplot2::fortify(windscape_example("wind_field"))
      names(d)[3:4] <- c("east", "north")
      p <- ggplot2::ggplot(d, ggplot2::aes(x, y)) + geom_wind_arrow(ggplot2::aes(u = east, v = north))
      expect_equal(ggplot2::layer_data(p)$speed, ggplot2::layer_data(arrow_plot())$speed)
})

test_that("fixed_length gives equal arrow lengths in km, with zero for calm blocks", {
      f <- windscape_example("wind_field")
      ld <- ggplot2::layer_data(arrow_plot(f, fixed_length = TRUE))
      km <- sqrt(((ld$xend - ld$x) * 111.32 * cos(ld$y0 * pi / 180))^2 + ((ld$yend - ld$y) * 110.57)^2)
      expect_equal(km, rep(km[1], length(km)), tolerance = 1e-8)
      # direction is unchanged from speed-scaled arrows
      ls <- ggplot2::layer_data(arrow_plot(f))
      expect_equal(atan2(ld$xend - ld$x, ld$yend - ld$y), atan2(ls$xend - ls$x, ls$yend - ls$y), tolerance = 1e-10)
      calm <- ggplot2::layer_data(arrow_plot(uniform_field(0, 0), fixed_length = TRUE, res = 5))
      expect_equal(calm$x, calm$xend)
      expect_equal(calm$y, calm$yend)
})

test_that("speed can be mapped to colour or linewidth", {
      p <- arrow_plot(mapping = ggplot2::aes(colour = ggplot2::after_stat(speed)), fixed_length = TRUE) +
            ggplot2::scale_colour_viridis_c()
      expect_gt(length(unique(ggplot2::layer_data(p)$colour)), 10)
      p2 <- arrow_plot(mapping = ggplot2::aes(linewidth = ggplot2::after_stat(speed)), fixed_length = TRUE)
      lw <- ggplot2::layer_data(p2)
      expect_true(all(diff(lw$linewidth[order(lw$speed)]) >= -1e-12)) # width increases with speed
      for(q in list(p, p2)) expect_s3_class(ggplot2::ggplotGrob(q), "gtable")
})
