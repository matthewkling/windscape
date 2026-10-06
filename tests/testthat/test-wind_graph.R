# wind_graph(), pairwise_least_cost(), least_cost() ------------------------

# cell-center coordinates of row `row`, columns `cols`, of raster r
centers <- function(r, row, cols) terra::xyFromCell(r, terra::cellFromRowCol(r, row, cols))

test_that("wind_graph returns a downwind or upwind wind_graph", {
      r <- uv_rose()
      g <- wind_graph(r)
      expect_s4_class(g, "wind_graph")
      expect_s4_class(g, "TransitionLayer")
      expect_equal(g@direction, "downwind")
      expect_equal(wind_graph(r, direction = "upwind")@direction, "upwind")
})

test_that("downwind travel time along a uniform westerly is exact", {
      r <- uv_rose(nr = 3, nc = 6, u = 5, v = 0)
      xy <- centers(r, 2, c(1, 6))
      d <- pairwise_least_cost(r, xy, snap = TRUE)
      dE <- geosphere::distGeo(c(0, xy[1, 2]), c(1, xy[1, 2])) # m per hop
      expect_equal(d[1, 2], 5 * dE / (5 * 3600), tolerance = 1e-8) # hours, 5 hops at 5 m/s
      expect_equal(d[2, 1], Inf)                                    # no westward wind
      expect_equal(unname(diag(d)), c(0, 0))
})

test_that("upwind least_cost gives travel time to the site", {
      r <- noisy_rose()
      xy <- terra::xyFromCell(r, c(5, 40, 77, 160))
      d <- pairwise_least_cost(r, xy, snap = TRUE)
      up <- least_cost(r, xy[1, , drop = FALSE], direction = "upwind")
      expect_equal(terra::extract(up, xy)[, 1], d[, 1], tolerance = 1e-8)
})

test_that("least-cost functions accept a prebuilt graph with a matching direction", {
      r <- noisy_rose()
      site <- terra::xyFromCell(r, 77)
      up <- wind_graph(r, direction = "upwind")
      expect_equal(terra::values(least_cost(up, site, direction = "upwind")),
                   terra::values(least_cost(r, site, direction = "upwind")))
      expect_error(least_cost(up, site), "upwind wind_graph")
      expect_error(least_cost_paths(up, site, n = 10), "upwind wind_graph")
      expect_error(least_cost(up, site, direction = "upwind", wrap = TRUE), "already a wind_graph")
      expect_error(least_cost(methods::as(r, "SpatRaster"), site), "wind_rose")
      expect_error(least_cost(r, site, direction = "sideways"), "should be one of")
})

test_that("rate = TRUE returns inverse cost distances", {
      r <- noisy_rose()
      xy <- terra::xyFromCell(r, c(5, 40, 77))
      g <- r
      expect_equal(pairwise_least_cost(g, xy, snap = TRUE, rate = TRUE),
                   1 / pairwise_least_cost(g, xy, snap = TRUE))
})

test_that("least_cost agrees with pairwise_least_cost", {
      r <- noisy_rose()
      xy <- terra::xyFromCell(r, c(5, 40, 77))
      g <- r
      s <- least_cost(g, xy[1, , drop = FALSE])
      expect_s4_class(s, "SpatRaster")
      expect_equal(names(s), "hours")
      d <- pairwise_least_cost(g, xy, snap = TRUE)
      expect_equal(terra::extract(s, xy)[, 1], d[1, ], tolerance = 1e-8)
      sr <- least_cost(g, xy[1, , drop = FALSE], rate = TRUE)
      expect_equal(names(sr), "rate")
      expect_equal(terra::values(sr), 1 / terra::values(s), ignore_attr = TRUE)
})

test_that("wrap connects the eastern and western edges", {
      r <- uv_rose(nr = 3, nc = 4, u = 5, v = 0)
      xy <- centers(r, 2, c(4, 1))
      open <- pairwise_least_cost(r, xy, snap = TRUE)
      wrapped <- pairwise_least_cost(r, xy, snap = TRUE, wrap = TRUE)
      dE <- geosphere::distGeo(c(0, xy[1, 2]), c(1, xy[1, 2]))
      expect_equal(open[1, 2], Inf)
      expect_equal(wrapped[1, 2], dE / (5 * 3600), tolerance = 1e-8) # one hop across the suture
})

# least_cost_paths() ------------------------

test_that("path costs equal least-cost distances, with vertices from origin to destination", {
      r <- noisy_rose()
      g <- wind_graph(r)
      from <- terra::xyFromCell(r, c(5, 77))
      to <- terra::xyFromCell(r, c(40, 120, 160))
      p <- least_cost_paths(g, from, to)
      expect_named(p, c("trail", "site", "to", "step", "hours", "x", "y"))
      expect_equal(length(unique(p$trail)), 6)                    # all 2 x 3 pairs
      cost <- matrix(gdistance::costDistance(g, from, to), 2, 3)
      fin <- p[!duplicated(p$trail, fromLast = TRUE), ]
      expect_equal(fin$hours, cost[cbind(fin$site, fin$to)], tolerance = 1e-8)
      for(d in split(p, p$trail)){
            expect_equal(terra::cellFromXY(r, as.matrix(d[1, c("x", "y")])), terra::cellFromXY(r, from[d$site[1], , drop = FALSE]))
            expect_equal(terra::cellFromXY(r, as.matrix(d[nrow(d), c("x", "y")])), terra::cellFromXY(r, to[d$to[1], , drop = FALSE]))
            expect_false(is.unsorted(d$hours))
            expect_equal(d$step, seq_len(nrow(d)) - 1)
      }
})

test_that("a uniform westerly gives straight eastward paths; westward pairs are omitted", {
      r <- uv_rose(nr = 3, nc = 6, u = 5, v = 0)
      g <- wind_graph(r)
      west <- terra::xyFromCell(r, terra::cellFromRowCol(r, 2, 1))
      east <- terra::xyFromCell(r, terra::cellFromRowCol(r, 2, 6))
      p <- least_cost_paths(r, west, east)
      expect_equal(nrow(p), 6)
      expect_true(all(diff(p$x) > 0))
      expect_equal(unique(p$y), unname(west[, 2]))
      expect_warning(q <- least_cost_paths(r, east, west), "no path")
      expect_equal(nrow(q), 0)
})

test_that("paths are directional", {
      r <- noisy_rose()
      g <- wind_graph(r)
      a <- terra::xyFromCell(r, 5)
      b <- terra::xyFromCell(r, 160)
      ab <- least_cost_paths(g, a, b)
      ba <- least_cost_paths(g, b, a)
      expect_false(isTRUE(all.equal(max(ab$hours), max(ba$hours))))
})

test_that("upwind paths run toward the sites, in the direction of travel", {
      r <- noisy_rose()
      a <- terra::xyFromCell(r, 5)
      b <- terra::xyFromCell(r, 160)
      down <- least_cost_paths(r, a, b)                        # a to b, downwind
      up <- least_cost_paths(r, b, a, direction = "upwind")    # also a to b, traced from b
      expect_equal(up[c("step", "hours", "x", "y")], down[c("step", "hours", "x", "y")], tolerance = 1e-8)
      expect_equal(up$site[1], 1)
      expect_equal(up$to[1], 1)
})

test_that("nearest and matched pair options", {
      r <- noisy_rose()
      g <- wind_graph(r)
      from <- terra::xyFromCell(r, c(5, 77, 100))
      to <- terra::xyFromCell(r, c(40, 120, 160))
      cost <- matrix(gdistance::costDistance(g, from, to), 3, 3)
      n <- least_cost_paths(g, from, to, pairs = "nearest")
      chosen <- unique(n[, c("site", "to")])
      expect_equal(nrow(chosen), 3)
      expect_equal(chosen$to, unname(apply(cost, 1, which.min)))
      m <- least_cost_paths(g, from, to, pairs = "matched")
      expect_equal(unique(m[, c("site", "to")]), data.frame(site = 1:3, to = 1:3), ignore_attr = TRUE)
      expect_error(least_cost_paths(g, from, to[1:2, ], pairs = "matched"), "same number of rows")
})

test_that("input handling: data frames, SpatVectors, same-cell pairs, validation", {
      r <- noisy_rose()
      g <- wind_graph(r)
      from <- terra::xyFromCell(r, 5)
      to <- terra::xyFromCell(r, c(5, 160))     # first destination is the origin's cell
      p <- least_cost_paths(g, as.data.frame(from), terra::vect(to, crs = "EPSG:4326"))
      expect_equal(unique(p$to), 2)
      expect_error(least_cost_paths(methods::as(r, "SpatRaster"), from, to), "wind_rose")
      expect_error(least_cost_paths(g, cbind(1, 2, 3), to), "two-column")
      expect_error(least_cost_paths(g, from, to, pairs = "some"))
})

test_that("paths draw with geom_wind_path", {
      r <- noisy_rose()
      g <- wind_graph(r)
      p <- least_cost_paths(g, terra::xyFromCell(r, 5), terra::xyFromCell(r, c(40, 160)))
      plt <- ggplot2::ggplot(p, ggplot2::aes(x, y)) + geom_wind_path(ggplot2::aes(color = hours))
      ld <- ggplot2::layer_data(plt)
      expect_equal(length(unique(ld$group)), 2)
      expect_s3_class(ggplot2::ggplotGrob(plt), "gtable")
})

test_that("least_cost_paths defaults to a grid of destinations", {
      rose <- windscape_example("wind_rose")
      site <- cbind(-105, 40)
      site_x <- terra::xFromCell(rose, terra::cellFromXY(rose, site))
      expect_silent(p <- least_cost_paths(rose, site, n = 60))
      k <- length(unique(p$trail))
      expect_gt(k, 40)
      expect_lte(k, 75)
      expect_true(all(p$x[p$step == 0] == site_x))
      expect_equal(least_cost_paths(wind_graph(rose), site, n = 60), p)
      up <- least_cost_paths(rose, site, n = 60, direction = "upwind")
      expect_true(all(up$x[!duplicated(up$trail, fromLast = TRUE)] == site_x)) # paths end at the site
      expect_error(least_cost_paths(rose, site, pairs = "matched"), "requires `to`")
})

test_that("wind_graph edges average the two cells' conductance in the edge direction", {
      rose <- windscape_example("wind_rose")
      v <- terra::values(rose)
      tm <- gdistance::transitionMatrix(wind_graph(rose))
      nc <- terra::ncol(rose)
      a <- 10 * nc + 20                       # an interior cell
      e <- a + 1                              # its eastern neighbor
      s <- a + nc                             # its southern neighbor
      expect_equal(tm[a, e], mean(v[c(a, e), "E"]))
      expect_equal(tm[a, s], mean(v[c(a, s), "S"]))
      expect_equal(tm[a, s + 1], mean(v[c(a, s + 1), "SE"]))
})

test_that("an upwind graph is the transpose of the downwind graph", {
      rose <- windscape_example("wind_rose")
      dn <- gdistance::transitionMatrix(wind_graph(rose))
      up <- gdistance::transitionMatrix(wind_graph(rose, "upwind"))
      expect_equal(up, Matrix::t(dn))
})

test_that("wrap adds edges across the left and right edges, skipping NA cells", {
      rose <- windscape_example("wind_rose")
      nc <- terra::ncol(rose)
      first <- 5 * nc + 1                     # first column of row 6
      last <- 6 * nc                          # last column of row 6
      expect_equal(gdistance::transitionMatrix(wind_graph(rose))[last, first], 0)
      expect_gt(gdistance::transitionMatrix(wind_graph(rose, wrap = TRUE))[last, first], 0)

      holes <- rose
      v <- terra::values(holes)
      v[first, ] <- NA
      terra::values(holes) <- v
      tm <- gdistance::transitionMatrix(wind_graph(holes, wrap = TRUE))
      expect_false(anyNA(tm@x))
      expect_equal(sum(tm[first, ]) + sum(tm[, first]), 0)
})

test_that("wind_graph rejects an unknown direction", {
      expect_error(wind_graph(windscape_example("wind_rose"), "sideways"), "should be one of")
})

test_that("sites and points outside the grid are handled, not misaligned", {
      r <- noisy_rose()
      xy <- terra::xyFromCell(r, c(5, 40, 77))
      out <- cbind(0, 0)
      d <- pairwise_least_cost(r, xy, snap = TRUE)
      expect_warning(d2 <- pairwise_least_cost(r, rbind(xy[1, ], out, xy[2:3, ]), snap = TRUE), "outside")
      expect_equal(d2[-2, -2], d, ignore_attr = TRUE)
      expect_true(all(is.na(d2[2, ])) && all(is.na(d2[, 2])))

      s <- least_cost(r, xy[1, , drop = FALSE])
      expect_warning(s2 <- least_cost(r, rbind(out, xy[1, ])), "outside")
      expect_equal(terra::values(s2), terra::values(s))
      expect_error(suppressWarnings(least_cost(r, out)), "all sites")

      p <- least_cost_paths(r, xy[1, , drop = FALSE], xy[2:3, ])
      expect_warning(p2 <- least_cost_paths(r, xy[1, , drop = FALSE], rbind(out, xy[2:3, ])), "outside")
      expect_equal(p2$to, p$to + 1)
      expect_equal(p2[c("step", "hours", "x", "y")], p[c("step", "hours", "x", "y")])
})
