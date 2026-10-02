# vrcd() ------------------------

test_that("vrcd refines only nearby site pairs", {
      r <- noisy_rose(nr = 8, nc = 10)
      ll <- cbind(c(-98.5, -98.3, -93.5), c(33.5, 33.6, 35.5)) # sites 1-2 ~22 km apart
      v <- quietly(vrcd(r, ll, threshold_km = 50, max_nodes = 400))
      expect_named(v, c("wind_dist", "wind_dist_coarse", "point_dist", "cell_dist", "cell_dist_coarse"))
      expect_equal(v$point_dist, point_distance(ll))
      expect_equal(v$wind_dist_coarse, least_cost_distance(wind_graph(r), ll, adjust = FALSE),
                   ignore_attr = TRUE)
      far <- cbind(c(1, 2, 3, 3), c(3, 3, 1, 2))
      expect_equal(v$wind_dist[far], v$wind_dist_coarse[far])
      expect_true(all(is.finite(v$wind_dist[1:2, 1:2])))
      expect_true(v$wind_dist[1, 2] != v$wind_dist_coarse[1, 2])
})
