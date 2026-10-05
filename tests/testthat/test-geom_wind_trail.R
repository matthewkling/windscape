# geom_wind_trail() / stat_wind_trail() ------------------------

trail_plot <- function(field = windscape_example("wind_field"), ...){
      ggplot2::ggplot(field, ggplot2::aes(x, y)) + stat_wind_trail(...) + ggplot2::coord_quickmap()
}
field_from <- function(u, v, nr = 10, nc = 10){
      r <- terra::rast(nrows = nr, ncols = nc, xmin = -100, xmax = -100 + nc, ymin = 30,
                       ymax = 30 + nr, crs = "EPSG:4326", nlyrs = 2)
      terra::values(r) <- cbind(rep_len(u, nr * nc), rep_len(v, nr * nc))
      wind_field(r)
}
trail_km <- function(d){ # length of each trail, in km
      vapply(split(d, d$group), function(z){
            sum(sqrt((diff(z$x) * 111.32 * cos(z$y[-1] * pi / 180))^2 + (diff(z$y) * 110.57)^2))
      }, numeric(1))
}
# trails stop where they leave the field; keep only trails with the full number of points
full_trails <- function(d){
      n <- table(d$group)
      d[d$group %in% as.integer(names(n)[n == max(n)]), ]
}

test_that("a uniform westerly gives eastward, equal-length trails from block seeds", {
      d <- ggplot2::layer_data(trail_plot(field_from(5, 0), res = 5, fixed_length = TRUE))
      expect_true(all(vapply(split(d, d$group), function(z) all(diff(z$x) > 0) && all(abs(diff(z$y)) < 1e-9),
                             logical(1))))
      expect_equal(length(unique(d$group)), nrow(block_means(ggplot2::fortify(field_from(5, 0)), c("u", "v"), 5)))
      km <- unname(trail_km(full_trails(d)))
      expect_equal(km, rep(km[1], length(km)), tolerance = 1e-6)
      expect_lt(max(trail_km(d)), max(km) + 1e-6) # truncated trails are only ever shorter
      expect_true(all(d$speed == 5))
      expect_true(all(d$progress >= 0 & d$progress <= 1))
})

test_that("length sets trail length relative to seed spacing", {
      f <- field_from(5, 0, nr = 20, nc = 20)
      a <- trail_km(full_trails(ggplot2::layer_data(trail_plot(f, res = 5, length = 1, fixed_length = TRUE))))
      b <- trail_km(full_trails(ggplot2::layer_data(trail_plot(f, res = 5, length = 2, fixed_length = TRUE))))
      expect_equal(median(b) / median(a), 2, tolerance = 0.02)
      expect_error(trail_plot(f, length = 0), "length")
      expect_error(trail_plot(f, res = 0), "res")
      expect_error(trail_plot(f, steps = 0), "steps")
})

test_that("with fixed_length = FALSE, length is proportional to speed", {
      # northern half 10 m/s, southern half 5 m/s, all westerly
      f <- field_from(rep(c(10, 5), each = 200), 0, nr = 20, nc = 20)
      d <- ggplot2::layer_data(trail_plot(f, res = 5, fixed_length = FALSE, direction = "downwind",
                                          length = 0.5))
      # trails seeded well inside the domain (west half), so none are truncated
      d <- d[d$group %in% unique(d$group[d$t == 0 & d$x < -92]), ]
      km <- trail_km(d)
      fast <- tapply(d$speed, d$group, min) == 10
      expect_equal(median(km[fast]) / median(km[!fast]), 2, tolerance = 0.05)
      e <- full_trails(ggplot2::layer_data(trail_plot(f, res = 5, direction = "downwind", fixed_length = TRUE)))
      ke <- unname(trail_km(e))
      expect_equal(ke, rep(ke[1], length(ke)), tolerance = 0.01)
})

test_that("direction controls which way trails are traced", {
      f <- windscape_example("wind_field")
      expect_true(all(ggplot2::layer_data(trail_plot(f, direction = "downwind"))$t >= 0))
      expect_true(all(ggplot2::layer_data(trail_plot(f, direction = "upwind"))$t <= 0))
      expect_true(any(ggplot2::layer_data(trail_plot(f))$t < 0))
})

test_that("facets share seeds and trails", {
      d <- ggplot2::fortify(windscape_example("wind_field"))
      dd <- rbind(cbind(d, version = "a"), cbind(d, version = "b"))
      p <- ggplot2::ggplot(dd, ggplot2::aes(x, y)) + stat_wind_trail(res = 10) + ggplot2::facet_wrap(~version)
      ld <- ggplot2::layer_data(p)
      expect_equal(ld[ld$PANEL == 1, c("x", "y")], ld[ld$PANEL == 2, c("x", "y")], ignore_attr = TRUE)
})

test_that("geom_wind_path draws wind_trails output, ordering points by step", {
      f <- windscape_example("wind_field")
      tr <- wind_trails(f, cbind(c(-92, -86), c(24, 30)), hours = 6, steps = 20)
      shuffled <- tr[sample(nrow(tr)), ]
      p <- ggplot2::ggplot(shuffled, ggplot2::aes(x, y)) + geom_wind_path()
      ld <- ggplot2::layer_data(p)
      expect_equal(length(unique(ld$group)), 2)
      expect_true(all(vapply(split(ld$t, ld$group), function(t) !is.unsorted(t), logical(1))))
      expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
      expect_s3_class(ggplot2::ggplotGrob(p + ggplot2::aes(colour = speed)), "gtable")
})

test_that("stat_wind_trail draws, with and without arrows", {
      for(p in list(trail_plot(), trail_plot(arrow = NULL),
                    trail_plot(mapping = ggplot2::aes(colour = ggplot2::after_stat(speed))))){
            expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
      }
})

test_that("trails default to speed-proportional length; hours sets transport time", {
      f <- field_from(5, 0, nr = 20, nc = 20)
      d <- full_trails(ggplot2::layer_data(trail_plot(f, res = 5, hours = 4, direction = "downwind")))
      expect_true("hours" %in% names(d))
      expect_equal(max(d$hours), 4)
      expect_equal(median(trail_km(d)), 5 * 3.6 * 4, tolerance = 0.01) # 72 km
      expect_true("km" %in% names(ggplot2::layer_data(trail_plot(f, fixed_length = TRUE))))
      expect_error(trail_plot(f, hours = -1), "hours")
})

test_that("seeds can be supplied; wrap is passed through", {
      f <- windscape_example("wind_field")
      sites <- cbind(c(-92, -86), c(24, 30))
      d <- ggplot2::layer_data(trail_plot(f, seeds = sites, hours = 6))
      expect_equal(length(unique(d$group)), 2)
      expect_equal(d[d$t == 0, c("x", "y")], data.frame(x = sites[, 1], y = sites[, 2]), ignore_attr = TRUE)
      expect_error(trail_plot(f, seeds = cbind(1, 2, 3)), "two-column")
      w <- ggplot2::layer_data(trail_plot(field_from(5, 0), seeds = cbind(-90.5, 35), hours = 20,
                                          direction = "downwind", wrap = "horizontal"))
      expect_gt(length(unique(w$group)), 1)
})

test_that("geom_wind_trail computes trails from a wind field, like stat_wind_trail", {
      f <- windscape_example("wind_field")
      base <- ggplot2::ggplot(f, ggplot2::aes(x, y))
      g <- ggplot2::layer_data(base + geom_wind_trail(res = 8))
      s <- ggplot2::layer_data(base + stat_wind_trail(res = 8))
      expect_gt(length(unique(g$group)), 10)
      expect_equal(g[c("x", "y", "group", "t", "speed")], s[c("x", "y", "group", "t", "speed")])
      expect_error(ggplot2::layer_data(base + geom_wind_trail(res = 0)), "res")
})

test_that("field layers reject multi-step series but allow one step per panel", {
      series <- windscape_example("wind_series")
      s2 <- subset_series(series, steps = 1:2)
      base <- ggplot2::ggplot(s2, ggplot2::aes(x, y))
      expect_error(ggplot2::layer_data(base + geom_wind_trail(res = 6)), "more than one wind vector")
      expect_error(ggplot2::layer_data(base + geom_wind_arrow(res = 6)), "more than one wind vector")

      # faceting by time gives one step per panel
      faceted <- base + ggplot2::facet_wrap(~time)
      expect_gt(nrow(ggplot2::layer_data(faceted + geom_wind_trail(res = 6))), 0)
      expect_gt(nrow(ggplot2::layer_data(faceted + geom_wind_arrow(res = 6))), 0)

      # so do groups, e.g. two fields overlaid by color
      d <- ggplot2::fortify(s2)
      grouped <- ggplot2::ggplot(d, ggplot2::aes(x, y, color = factor(step)))
      expect_gt(nrow(ggplot2::layer_data(grouped + geom_wind_arrow(res = 6))), 0)

      # and a one-step series plots like a wind field
      expect_gt(nrow(ggplot2::layer_data(ggplot2::ggplot(subset_series(series, steps = 1), ggplot2::aes(x, y)) +
                                               geom_wind_trail(res = 6))), 0)
})
