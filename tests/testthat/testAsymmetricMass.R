library("qad")
library("copula")

createDSMatrix <- function(family, parameter = 0, m = 5, n = 5) {
  cop <- switch(family,
                Cp=fhCopula("upper", dim = 2),
                Cm=fhCopula("lower", dim = 2),
                NC=normalCopula(parameter, dim = 2),
                FC=frankCopula(param = parameter, dim = 2),
                MO=moCopula(param=c(parameter, 0.4))
  )

  A <- matrix(0, nrow = m, ncol = n)
  for(k in 1:(m-1)) {
    for (j in 1:(n-1)) {
      A[k, j] <- pCopula(c(k/m, j/n), cop) -pCopula(c((k-1)/m, j/n), cop) - pCopula(c(k/m, (j-1)/n), cop) + pCopula(c((k-1)/m, (j-1)/n), cop)
    }
    A[k, n] <- 1/m - pCopula(c(k/m, (n-1)/n), cop) + pCopula(c((k-1)/m, (n-1)/n), cop)
  }
  for (j in 1:(n-1)) {
    A[m, j] <- 1/n -pCopula(c((m-1)/m, j/n), cop) + pCopula(c((m-1)/m, (j-1)/n), cop)
  }
  A[m, n] <- 1 - (m-1)/m - (n-1)/n + pCopula(c((m-1)/m, (n-1)/n), cop)
  return(A)
}

test_that("Checkerboard mass (Symmetric case)", {
  n <- 100
  X <- expand.grid(seq(0, 1, length.out = n), seq(0, 1, length.out = n))

  A <- ECBC(X[, 1], X[, 2])
  B <- checkerboardDensity(X[, 1], X[, 2], sqrt(n), sqrt(n))
  expect_equal(A, B)
})

test_that("Estimate checkerboard mass for independence (Asymmetric case)", {
  m <- 7^2
  n <- 13^2
  X <- expand.grid(seq(0, 1, length.out = m), seq(0, 1, length.out = n))

  A <- checkerboardDensity(X[, 1], X[, 2], resolution1=sqrt(m), resolution2=sqrt(n))
  expect_equal(A*sqrt(m*n), matrix(1, nrow=7, ncol=13))
})


test_that("Estimate checkerboard mass for complete dependence (Asymmetric case)", {
  m <- 7^2
  n <- 13^2
  u <- seq(0, 1, length.out=n)
  X <- cbind(u, u)

  A <- checkerboardDensity(X[, 1], X[, 2], resolution1=sqrt(m), resolution2=sqrt(n))
  expect_equal(A, createDSMatrix("Cp", 0, 7, 13))
})


test_that("Checkerboard mass index (Asymmetric)", {
  m <- 7^2
  n <- 13^2
  u <- seq(0, 1, length.out=n)
  X <- cbind(u, u)

  A <- checkerboardDensity(X[, 1], X[, 2], resolution1=sqrt(m), resolution2=sqrt(n))
  U <- rank(X[, 1], ties.method = "max")
  V <- rank(X[, 2], ties.method = "max")
  for(k in 1:sqrt(m)) {
	  for(l in 1:sqrt(n)) {
		  expect_equal(A[k, l], checkerboardDensityIndex(U, V, k, l, resolution1=sqrt(m), resolution2=sqrt(n)))
	  }
  }

})
