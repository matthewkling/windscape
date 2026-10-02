# mantel_test() ------------------------

sym_matrix <- function(x){
      n <- which(cumsum(1:100) == length(x)) + 1
      m <- matrix(0, n, n)
      m[lower.tri(m)] <- x
      m[upper.tri(m)] <- t(m)[upper.tri(t(m))]
      m
}

test_that("statistic is the correlation of off-diagonal elements", {
      set.seed(1)
      x <- matrix(runif(64), 8)
      y <- matrix(runif(64), 8)
      off <- row(x) != col(x)
      m <- mantel_test(x, y, nperm = 9)
      expect_equal(m$stat, cor(x[off], y[off]))
      expect_equal(mantel_test(x, y, nperm = 9, method = "spearman")$stat,
                   cor(x[off], y[off], method = "spearman"))
      expect_length(m$perm, 9)
})

test_that("partial statistic is the correlation of residuals on controls", {
      set.seed(2)
      x <- sym_matrix(runif(45))
      y <- sym_matrix(runif(45))
      z1 <- sym_matrix(runif(45))
      z2 <- sym_matrix(runif(45))
      off <- row(x) != col(x)
      Z <- cbind(z1[off], z2[off])
      expected <- cor(residuals(lm(x[off] ~ Z)), residuals(lm(y[off] ~ Z)))
      expect_equal(mantel_test(x, y, list(z1, z2), nperm = 9)$stat, expected)
})

test_that("strong association is significant and p-values follow the alternative", {
      set.seed(3)
      x <- sym_matrix(runif(45))
      y <- x + sym_matrix(rnorm(45, 0, 0.05))
      m <- mantel_test(x, y, nperm = 199)
      expect_lt(m$p.value, 0.02)
      expect_equal(mantel_test(x, y, nperm = 199, alternative = "greater")$p.value, 1 - m$quantile,
                   tolerance = 0.05)
      expect_gt(mantel_test(x, y, nperm = 199, alternative = "less")$p.value, 0.95)
})

test_that("p-values are roughly uniform under the null", {
      set.seed(4)
      p <- replicate(40, mantel_test(sym_matrix(runif(45)), sym_matrix(runif(45)), nperm = 99)$p.value)
      expect_gt(mean(p < 0.5), 0.25)
      expect_lt(mean(p < 0.5), 0.75)
})
