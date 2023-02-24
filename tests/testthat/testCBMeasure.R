library("copula")
library("qad")

kendall_naive <- function(A) {
  m <- dim(A)[1]
  n <- dim(A)[2]
  Q <- 0
  for (k in 1:m) {
    for (j in 1:n) {
      for (i in 1:k) {
        for (l in 1:j) {
          p1 <- ifelse(k == i, 0.5, 1);
          p2 <- ifelse(j == l, 0.5, 1);
          Q <- Q + A[k, l]*A[i, j]*p1*p2
        }
      }
    }
  }
  return(1 - 4*Q)
}

# createDSMatrix <- function(family, parameter = 0, m = 5, n = 5) {
#   cop <- switch(family,
#                 Cp=fhCopula("upper", dim = 2),
#                 Cm=fhCopula("lower", dim = 2),
#                 NC=normalCopula(parameter, dim = 2),
#                 FC=frankCopula(param = parameter, dim = 2),
#                 MO=moCopula(param=c(parameter, 0.4))
#   )
#
#   A <- matrix(0, nrow = m, ncol = n)
#   for(k in 1:m) {
#     for (j in 1:n) {
#       A[k, j] <- pCopula(c(k/m, j/n), cop) -pCopula(c((k-1)/m, j/n), cop) - pCopula(c(k/m, (j-1)/n), cop) + pCopula(c((k-1)/m, (j-1)/n), cop)
#     }
#   }
#   return(A)
# }
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


test_that("independence (Spearman)", {
  dVec <- c(1, 2, 3, 25)
  for (d in dVec) {
    A <- matrix( rep(1, len=d*d), nrow=d)/d^2
    expect_equal(computeCBMeasure(A, method="spearman"), 0)
  }
})

test_that("complete negative dependence (Spearman)", {
  dVec <- c(2, 3, 25)
  for (d in dVec) {
    A <- diag(d)[d:1, ]/d
    expect_equal(computeCBMeasure(A, method="spearman"), -(d^2 - 1)/(d^2))
  }
})

test_that("complete positive dependence (Spearman)", {
  dVec <- c(1, 2, 3, 25)
  for (d in dVec) {
    A <- diag(d)/d
    expect_equal(computeCBMeasure(A, method="spearman"), (d^2 - 1)/(d^2))
  }
})

test_that("convex combinations (Spearman)", {
  combVec <- seq(from=0, to=1, length.out=3)
  for (comb in combVec) {
    dVec <- c(1, 2, 3, 25)
    for (d in dVec) {
      A1 <- matrix( rep(1, len=d*d), nrow=d)/d^2
      A2 <- diag(d)/d
      expect_equal(computeCBMeasure((1-comb)*A1 + (comb)*A2, method="spearman"), comb*(d^2 - 1)/(d^2))
    }
  }
})

test_that("independence (Kendall)", {
  dVec <- c(1, 2, 3, 25)
  for (d in dVec) {
    A <- matrix( rep(1, len=d*d), nrow=d)/(d^2)
    expect_equal(computeCBMeasure(A, method="kendall"), 0, info=paste0("d = ", d))
  }
})

test_that("complete negative dependence (Kendall)", {
  dVec <- c(2, 3, 25)
  for (d in dVec) {
    A <- diag(d)[d:1, ]/(d)
    expect_equal(computeCBMeasure(A, method="kendall"), -(d - 1)/d)
  }
})

test_that("complete positive dependence (Kendall)", {
  dVec <- c(1, 2, 3, 25)
  for (d in dVec) {
    A <- diag(d)/d
    expect_equal(computeCBMeasure(A, method="kendall"), (d - 1)/d)
  }
})

test_that("Gaussian dependence (Kendall)", {
  d1Vec <- c(3, 5, 10)
  d2Vec <- c(3, 4, 13)
  for(d1 in d1Vec) {
    for (d2 in d2Vec) {
      A <- createDSMatrix("NC", 0.5, d1, d2)
      expect_equal(computeCBMeasure(A, method="kendall"), kendall_naive(A))
    }
  }
})

test_that("independence (DSS)", {
  dVec <- c(1, 2, 3, 25)
  for (d in dVec) {
    for (d2 in dVec) {
      A <- matrix( rep(1, len=d*d2), nrow=d)/(d*d2)
      expect_equal(computeCBMeasure(A, method="dss"), 0)
    }
  }
})

test_that("complete positive dependence (DSS)", {
  dVec <- c(2, 3, 19)
  for (d in dVec) {
    A <- diag(d)/(d)
    expect_equal(computeCBMeasure(A, method="dss"), (d-1)/d)
  }
})

test_that("DSS and Chatterjee coincide", {
  dVec <- c(2, 3, 25)
  for (d in dVec) {
    A <- createDSMatrix("NC", parameter = 0.6, m = 2*d, n=d)
    expect_equal(computeCBMeasure(A, method="dss"), computeCBMeasure(A, method="chatterjee"))
  }
})

test_that("complete positive dependence (DSS, asymmetric)", {
  dVec <- c(2, 4, 8)
  for (d in dVec) {
    A <- createDSMatrix("Cp", m = 2*d, n=d)
    expect_equal(computeCBMeasure(A, method="dss"), 1 - 1/d)
  }
})

test_that("Gaussian dependence (DSS, asymmetric)", {
  N <- 25
  for(p in c(0.25, 0.5, 0.75)) {
    A <- createDSMatrix("NC", p, N, N)
    expect_lte(abs(computeCBMeasure(A, method="dss") -  3*asin((1+p^2)/2)/pi + 0.5), 0.05)
  }
})

test_that("independence (BKR)", {
  dVec <- c(1, 2, 3, 25)
  for (d in dVec) {
    for (d2 in dVec) {
      A <- matrix( rep(1, len=d*d2), nrow=d)/(d*d2)
      expect_equal(computeCBMeasure(A, method="bkr"), 0)
    }
  }
})

test_that("independence (Zeta_1)", {
  dVec <- c(1, 2, 3, 25)
  for (d in dVec) {
    for (d2 in dVec) {
      A <- matrix( rep(1, len=d*d2), nrow=d)/(d*d2)
      expect_equal(computeCBMeasure(A, method="zeta1"), 0)
    }
  }
})

test_that("negative dependence (Zeta_1)", {
  A <- diag(20)/20
  expect_equal(computeCBMeasure(A, method="zeta1"), 3*D1.ECBC(A, matrix(1/400, nrow=20, ncol=20)))
})




test_that("checks implemented", {
  A <- data.frame(1, 2, 3)
  expect_error(computeCBMeasure(A, checkInput = TRUE))
})
