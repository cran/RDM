test_that("Estimate RDM for independence (Fixed)", {
  nVector <- c(2, 5, 13);
  for(n in nVector) {
    X <- expand.grid(seq(0, 1, length.out = n^2), seq(0, 1, length.out = n^2))
    r <- rdm(X, method = "spearman", bandwidth_method = "fixed", bandwidth_parameter = 0.5)
    expect_equal(r, 0)
  }
})

test_that("Estimate RDM for independence (Symmetric Cross Validation)", {
  nVector <- c(3, 5);
  for(n in nVector) {
    X <- expand.grid(seq(1/n^2, 1, length.out = (n^2-1)), seq(1/n^2, 1, length.out = (n^2-1)))
    r <- rdm(X, method = "spearman", bandwidth_method = "cvsym", bandwidth_parameter = c(0.5, 0.5))
    expect_equal(r, 0)
  }
})


test_that("Estimate RDM for complete positive dependence (Fixed, symmetric)", {
  nVector <- c(2, 5, 13);
  for(n in nVector) {
    u <- seq(0, 1, length.out=n^2)
    X <- cbind(u, u)

    r <- rdm(X, method = "spearman", bandwidth_method = "fixed", bandwidth_parameter = 0.5)
    expect_equal(r, (n^2 - 1)/(n^2))
  }
})

test_that("All methods (RDM)", {
  X <- cbind(runif(50), runif(50))
  result <- rdm(X, method="all", bandwidth_method = "fixed", bandwidth_parameter = 0.4)
  expect_equal(result$spearman, rdm(X, method="spearman", bandwidth_method = "fixed", bandwidth_parameter = 0.4))
  expect_equal(result$kendall, rdm(X, method="kendall", bandwidth_method = "fixed", bandwidth_parameter = 0.4))
})

test_that("Permutation test (RDM)", {
  X <- cbind(1:100, 1:100)
  result <- rdm(X, method="spearman", bandwidth_method = "fixed", bandwidth_parameter = 0.5)
  expect_lt(rdm(X, method="spearman", bandwidth_method = "fixed", bandwidth_parameter = 0.4, permutation = TRUE)$pvalue, 0.05)
})

test_that("Invalid bandwidth (RDM)", {
  n <- 5
  u <- seq(0, 1, length.out=n^2)
  X <- cbind(u, u)
  expect_error(rdm(X, method = "spearman", bandwidth_method = "fixed", bandwidth_parameter = 'a', checkInput=TRUE))
  expect_error(rdm(X, method = "spearman", bandwidth_method = "fixed", bandwidth_parameter = -0.5, checkInput=TRUE))
  expect_error(rdm(X, method = "spearman", bandwidth_method = "fixed", bandwidth_parameter = 0.51, checkInput=TRUE))
  expect_error(rdm(X, method = "spearman", bandwidth_method = "cv", bandwidth_parameter = c(0.5, 0.4), checkInput=TRUE))
  expect_error(rdm(X, method = "spearman", bandwidth_method = "cv", bandwidth_parameter = 0.5, checkInput=TRUE))
  expect_error(rdm(X, method = "spearman", bandwidth_method = "cv", bandwidth_parameter = c(0.4, 0.6), checkInput=TRUE))
  expect_error(rdm(X, method = "spearman", bandwidth_method = "cv", bandwidth_parameter = c(-0.1, 0.4), checkInput=TRUE))
  expect_error(rdm(X, method = "spearman", bandwidth_method = "cvsym", bandwidth_parameter = c(0.5, 0.4), checkInput=TRUE))
  expect_error(rdm(X, method = "spearman", bandwidth_method = "cvsym", bandwidth_parameter = c(0.4, 0.6), checkInput=TRUE))
  expect_error(rdm(X, method = "spearman", bandwidth_method = "cvsym", bandwidth_parameter = c(-0.1, 0.4), checkInput=TRUE))
  expect_error(rdm(X, method = "all", bandwidth_method = "fixed", bandwidth_parameter = 0.5, permutation = TRUE, checkInput=TRUE))
})

