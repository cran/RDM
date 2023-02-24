library("copula")

test_that("Check bandwidth equality", {
  nVector <- c(3, 5);
  for(n in nVector) {
    X <- expand.grid(seq(1/n^2, 1, length.out = (n^2-1)), seq(1/n^2, 1, length.out = (n^2-1)))
    N <- computeBandwidth(X, sL = 0.5, sU = 0.5, method="cvasym")
    M <- computeBandwidth(X, sL = 0.5, sU = 0.5, method="cvsym")
    expect_equal(N[1], N[2])
    expect_equal(N[1], M)
  }
})


test_that("Test pobs replacement", {

  n <- 100
  X <- rCopula(n, normalCopula(.5))
  U <- pobs(X[, 1])
  expect_equal(U, rank(X[, 1], ties.method = "max")/(n+1))
  u <- runif(n)
  X <- cbind(u, u)
  U <- pobs(X[, 1])
  expect_equal(U, rank(X[, 1], ties.method = "max")/(n+1))
})
