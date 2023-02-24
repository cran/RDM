#' @title Compute bandwidth via cross-validation
#' @description An implementation of the cross-validation principle for the bandwidth selection as presented in Strothmann, Dette and Siburg (2022) <arXiv:2201.03329>.
#' @details This function computes the optimal bandwidth given the bivariate observations \eqn{X} of length \eqn{N}.
#' Currently, there are two different algorithms implemented:
#' \itemize{
#'   \item "cvsym" - Computes the optimal bandwidth choice for a square checkerboard mass density according to the cross-validation principle. The bandwidth is a natural number between \eqn{N^{sL}, ..., N^{sU}}
#'   \item "cvasym" - Computes the optimal bandwidth choice \eqn{(N_1, N_2)} for a non-square checkerboard mass density according to the cross-validation principle.
#'   The bandwidths \eqn{N_1, N_2} are natural numbers between \eqn{N^{sL}, ..., N^{sU}} and may possibly attain different values.
#' }
#' @param X A bivariate data.frame containing the observations. Each row contains one observation.
#' @param sL Lower bound \eqn{N^{sL}} for the possible bandwidth parameters (where \eqn{N} is the number of observations).
#' @param sU Upper bound \eqn{N^{sU}} for the possible bandwidth parameters (where \eqn{N} is the number of observations).
#' @param method "cvsym" uses either a symmetric cross-validation principle (N_1 = N_2) and "cvasym" uses an asymmetric cross-validation principle (i.e. \eqn{N_1} and \eqn{N_2} may attain different values).
#' @param reduce In case reduce is set to TRUE, the parameter is chosen from N, N+2, ... instead of N, N+1, N+2, ...
#' @return The chosen bandwidth depending on the data.frame X.
#' @export computeBandwidth
#' @examples
#' n <- 20
#' X <- cbind(runif(n), runif(n))
#' computeBandwidth(X, sL = 0.25, sU = 0.5, method="cvsym", reduce=TRUE)
computeBandwidth <- function(X, sL, sU, method = c("cvsym", "cvasym"), reduce=TRUE) {
  type <- match.arg(method, c("cvsym", "cvasym"))

  if(type == "cvsym") {
    return(.computeSymmetricCrossValidation(X, sL, sU))
  }
  return(.computeAsymmetricCrossValidation(X, sL, sU, reduce=reduce))
}


#Private function to compute the asymmetric version, i.e. calculation of N1 and N2, of the cross validation principle.
.computeAsymmetricCrossValidation <- function(X, sL, sU, reduce = TRUE) {
  n <- dim(X)[1];

  lb <- floor(n^(sL))
  if(lb <= 2) {
    lb <- 2
  }
  ub <- floor(n^(sU))
  if(ub <= lb) {
    ub <- lb
  }
  nVector <- lb:ub;
  lN <- length(nVector)

  if(lN > 10 && reduce) {
    nVector <- seq(lb, ub, 2)
    lN <- length(nVector)
  }

  results <- matrix(0, nrow = lN, ncol = lN)
  #Removed after testing, thus removing @importFrom copula "pobs"
  #U <- pobs(X[, 1])
  #V <- pobs(X[, 2])

  X[, 1] <- rank(X[, 1], ties.method = "max")
  X[, 2] <- rank(X[, 2], ties.method = "max")

  U <- X[, 1]/(n+1)
  V <- X[, 2]/(n+1)

  tempX <- matrix(0, nrow=(n-1), ncol=n)
  tempY <- matrix(0, nrow=(n-1), ncol=n)
  #Calculate the updated rank in case the kth element is removed. This is necessary due to the direct call to asymmetric_checkerboard_mass
  #Then, any element with a higher rank moves down one element, any element with a lower rank remains the same
  for (k in 1:n) {
    idx <- (X[, 1] >= X[k, 1])
    tempX[, k] <- (X[, 1] - idx)[-k]
    idx <- (X[, 2] >= X[k, 2])
    tempY[, k] <- (X[, 2] - idx)[-k]
  }


  for (k in 1:lN) {
    for (l in 1:lN) {
      N1 <- nVector[k]
      N2 <- nVector[l]
      tmp <- 0
      kVector <- ceiling(N1*U)
      lVector <- ceiling(N2*V)
      for (i in 1:n) {
		    tmp <- tmp + N1*N2*(round(asymmetric_checkerboard_index(tempX[, i], tempY[, i], kVector[i], lVector[i], N1, N2),15))
      }
      results[k, l] <- (N1*N2)*sum((round(asymmetric_checkerboard_mass(X[, 1], X[, 2], N1, N2),15))^2) - 2*tmp/n
    }
  }
  idx <- which(results == min(results), arr.ind = TRUE)
  return (c(nVector[idx[1]], nVector[idx[2]]))
}


#Private function to compute the symmetric version, i.e. calculation of N = N1 = N2, of the cross validation principle.
.computeSymmetricCrossValidation <- function(X, sL, sU) {
  n <- dim(X)[1];

  lb <- floor(n^(sL))
  ub <- floor(n^(sU))
  nVector <- lb:ub;
  lN <- length(nVector)

  results <- rep.int(0, lN)

#  U <- pobs(X[, 1])
#  V <- pobs(X[, 2])

  X[, 1] <- rank(X[, 1], ties.method = "max")
  X[, 2] <- rank(X[, 2], ties.method = "max")

  U <- X[, 1]/(n+1)
  V <- X[, 2]/(n+1)

  tempX <- matrix(0, nrow=(n-1), ncol=n)
  tempY <- matrix(0, nrow=(n-1), ncol=n)
  #Calculate the updated rank in case the kth element is removed. This is necessary due to the direct call to asymmetric_checkerboard_mass
  #Then, any element with a higher rank moves down one element, any element with a lower rank remains the same
  for (k in 1:n) {
    idx <- (X[, 1] >= X[k, 1])
    tempX[, k] <- (X[, 1] - idx)[-k]
    idx <- (X[, 2] >= X[k, 2])
    tempY[, k] <- (X[, 2] - idx)[-k]
  }

  for (k in 1:lN) {
    N <- nVector[k]
    tmp <- 0
    kVector <- ceiling(N*U)
    lVector <- ceiling(N*V)
    for (i in 1:n) {
		  tmp <- tmp + N*N*(round(asymmetric_checkerboard_index(tempX[, i], tempY[, i], kVector[i], lVector[i], N, N),15))
    }
    results[k] <- (N*N)*sum((round(asymmetric_checkerboard_mass(X[, 1], X[, 2], N, N),15))^2) - 2*tmp/n;
  }
  idx <- which(results == min(results), arr.ind = TRUE)

  return (nVector[idx])
}


