# downscale(), weight_conductance() ------------------------

test_that("nearest-neighbor downscaling replicates cells and scales conductance by fact", {
      r <- noisy_rose(nr = 4, nc = 5)
      d <- downscale(r, 3, method = "near")
      expect_equal(dim(d), c(12, 15, 8))
      expect_equal(terra::values(d)[1, ], terra::values(r)[1, ] * 3)
      # a fine cell in the interior of coarse cell (2, 3)
      expect_equal(terra::values(d)[terra::cellFromRowCol(d, 5, 8), ],
                   terra::values(r)[terra::cellFromRowCol(r, 2, 3), ] * 3)
})

test_that("downscaling preserves travel time across the domain", {
      r <- uv_rose(nr = 3, nc = 6, u = 5, v = 0)
      d <- downscale(r, 3, method = "near")
      xy <- terra::xyFromCell(r, terra::cellFromRowCol(r, 2, c(1, 6)))
      t_coarse <- pairwise_least_cost(wind_graph(r), xy, snap = TRUE)[1, 2]
      xy_fine <- terra::xyFromCell(d, terra::cellFromXY(d, xy))
      t_fine <- pairwise_least_cost(wind_graph(d), xy_fine, snap = TRUE)[1, 2]
      # same physical distance between the two fine cells as between coarse cell centers
      expect_equal(t_fine, t_coarse, tolerance = 1e-3)
})

test_that("downscale validates input", {
      r <- noisy_rose()
      expect_error(downscale(methods::as(r, "SpatRaster"), 2), "wind_rose")
      expect_error(downscale(r, c(2, 3)), "single")
})

test_that("weight_conductance multiplies conductance by a weight layer", {
      r <- noisy_rose(nr = 4, nc = 5)
      w <- terra::rast(r, nlyrs = 1, vals = rep(c(1, 0.1), length.out = 20))
      x <- weight_conductance(r, w)
      expect_equal(terra::values(x), terra::values(r) * rep(c(1, 0.1), length.out = 20),
                   ignore_attr = TRUE)
})
