library("copula")
suppressMessages(library("Rfast"))

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

naiveSortDSMatrix <- function(A) {
  N1 <- dim(A)[1]
  N2 <- dim(A)[2]
  #In case either dimension is less or equal to 1, there is no rearrangement necessary
  if(min(N1, N2) <= 1) {
    return(A)
  }

  cumSumA <- t(colCumSums(t(A)))
  #Sort each column of cumulative sums in decreasing order
  sortedA <- apply(cumSumA, 2, sort, decreasing=TRUE)

  #Calculate the corresponding (generalized) doubly stochastic matrix from the sorted cumulative sums
  sortedDS <- matrix(nrow=N1, ncol=N2)
  sortedDS[, 1] <- sortedA[, 1]

  for (l in 2:N2) {
    sortedDS[, l] <- sortedA[, l] - sortedA[, l-1]
  }
  return(sortedDS)
}

test_that("implementation", {
  nVec <- c(10, 50, 100)
  for(N in nVec) {
    U <- rCopula(N, normalCopula(0.5))
    A <- checkerboardDensity(U[, 1], U[, 2], resolution1=N, resolution2=(N/2))
    expect_equal(sortDSMatrix(A), naiveSortDSMatrix(A))
    U <- rCopula(N, fhCopula("upper", dim = 2))
    A <- checkerboardDensity(U[, 1], U[, 2], resolution1=(N/2), resolution2=(N/2))
    expect_equal(sortDSMatrix(A), naiveSortDSMatrix(A))
    U <- rCopula(N, moCopula(param=c(0, 0.4)))
    A <- checkerboardDensity(U[, 1], U[, 2], resolution1=(N), resolution2=(N/2))
    expect_equal(sortDSMatrix(A), naiveSortDSMatrix(A))
  }
})

test_that("independence", {
  dVec <- c(1, 2, 3, 25)
  for (d in dVec) {
    A <- matrix( rep(1, len=d*d), nrow=d)/d
    expect_equal(sortDSMatrix(A), A)
  }
})

test_that("complete positive dependence", {
  dVec <- c(1, 2, 3, 25)
  for (d in dVec) {
    A <- diag(d)
    expect_equal(sortDSMatrix(A), A)
  }
})

test_that("complete negative dependence", {
  dVec <- c(2, 3, 11, 25)
  for (d in dVec) {
    A <- diag(d)[d:1, ]
    expect_equal(sortDSMatrix(A), diag(d))
  }
})

test_that("convex combinations", {
  combVec <- seq(from=0, to=1, length.out=3)
  for (comb in combVec) {
    dVec <- c(1, 2, 3, 25)
    for (d in dVec) {
      A1 <- matrix( rep(1, len=d*d), nrow=d)/d
      A2 <- diag(d)[d:1, ]
      expect_equal(sortDSMatrix((1-comb)*A1 + (comb)*A2), (1-comb)*A1 + (comb)*diag(d))
    }
  }
})

test_that("gaussian dependence", {
  d1Vec <- c(3, 5, 10)
  d2Vec <- c(3, 4, 13)
  for(d1 in d1Vec) {
    for (d2 in d2Vec) {
      A <- createDSMatrix("NC", -0.55, d1, d2)
      B <- createDSMatrix("NC", 0.55, d1, d2)
      expect_equal(sortDSMatrix(A), B)
    }
  }
})

test_that("input checks implemented", {
  A <- data.frame(1, 2, 3)
  expect_error(sortDSMatrix(A, checkInput = TRUE))
})

