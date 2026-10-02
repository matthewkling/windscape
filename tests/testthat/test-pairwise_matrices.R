# pairwise_ratios(), pairwise_means() ------------------------

test_that("pairwise_ratios of a vector gives x[i] / x[j]", {
      x <- c(1, 2, 4)
      r <- pairwise_ratios(x, log = FALSE)
      expect_equal(r[2, 3], 0.5)
      expect_equal(r[3, 1], 4)
      expect_equal(diag(r), rep(1, 3))
      expect_equal(pairwise_ratios(x), log(r))
})

test_that("pairwise_ratios of a matrix gives x[i, j] / x[j, i]", {
      m <- matrix(c(1, 2, 3, 4, 1, 6, 8, 9, 1), 3)
      r <- pairwise_ratios(m, log = FALSE)
      expect_equal(r[1, 2], m[1, 2] / m[2, 1])
      lr <- pairwise_ratios(m)
      expect_equal(lr, -t(lr))   # log ratios are antisymmetric
})

test_that("pairwise_means symmetrizes a matrix", {
      m <- matrix(1:9, 3)
      p <- pairwise_means(m)
      expect_equal(p, t(p))
      expect_equal(p[1, 2], (m[1, 2] + m[2, 1]) / 2)
})
