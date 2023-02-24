#' @title Rearranged dependence measure
#' @description This function estimates the asymmetric dependence between \eqn{X} and \eqn{Y} using the rearranged dependence measure \eqn{R_\mu(X, Y)} for different possible underlying measures \eqn{\mu}.
#' A value of 0 characterizes independence of \eqn{X} and \eqn{Y}, while a value of 1 characterizes a functional relationship between \eqn{X} and \eqn{Y}, i.e. \eqn{Y = f(X)}.
#' @details This function estimates \eqn{R_\mu(X, Y)} using the empirical checkerboard mass density \eqn{A}.
#' To arrive at \eqn{R_\mu(X, Y)}, \eqn{A} is appropriately sorted and then evaluated for the underlying measure.
#' The estimated \eqn{R_\mu} always takes values between 0 and 1 with
#' \itemize{
#'  \item \eqn{R_\mu(X, Y) = 0} if and only if \eqn{X} and \eqn{Y} are independent.
#'  \item \eqn{R_\mu(X, Y) = 1} if and only if \eqn{Y = f(X)} for some measurable function \eqn{f}.
#' }
#' Currently, the following underlying measures are implemented:
#' \itemize{
#'  \item{"spearman"}{ Implements the concordance measure Spearman's \eqn{\rho} (which is identical to the \eqn{L_1}-Schweizer-Wolff-measure),}
#'  \item{"kendall"}{ Implements the concordance measure Kendall's \eqn{\tau},}
#'  \item{"bkr"}{ Implements the Blum–Kiefer–Rosenblatt \eqn{R}, also known as the \eqn{L^2}-Schweizer-Wolff-measure <doi:10.1214/aos/1176345528>,}
#'  \item{"dss"}{ Implements the Dette-Siburg-Stoimenov measure of complete dependence <doi:10.1111/j.1467-9469.2011.00767.x>, also known as Chatterjee's \eqn{\xi} <doi:10.1080/01621459.2020.1758115>,}
#'  \item{"zeta1"}{ Implements the \eqn{\zeta_1}-measure of complete dependence established by W. Trutschnig <doi:10.1016/j.jmaa.2011.06.013>.}
#' }
#' The estimation of the checkerboard mass density \eqn{A} depends on the choice of the bandwidth for the checkerboard copula.
#' For a detailed discussion of "cv" and "cvsym", see \code{\link{computeBandwidth}}.
#' @param X A bivariate data.frame containing the observations. Each row contains one bivariate observation.
#' @param method Options include "spearman", "kendall", "bkr", "dss", "chatterjee" and "zeta1".The option "all" returns the value for all aforementioned methods.
#' @param bandwidth_method A character string indicating the use of either a cross-validation principle (square or non-square) or a fixed bandwidth (oftentimes called resolution).
#' @param bandwidth_parameter A numerical vector which contains the necessary optional parameters for the exponent of the chosen bandwidth method.
#' In case of N observations, the bandwidth_parameter \eqn{(s_1, s_2)} determines a lower bound \eqn{N^{s_1}} and upper bound \eqn{N^{s_2}} for the cross-validation methods
#' or a single number s for the fixed bandwidth method resulting in \eqn{N^s}.
#' The parameters have to lie in \eqn{(0, 1/2)} and fulfil \eqn{s_1 < s_2}.
#' @param permutation Whether or not to perform a permutation test
#' @param npermutation Number of repetitions of the permutation test
#' @param checkInput Whether or not to perform validity checks of the input
#' @return The estimated value of the rearranged dependence measure
#' @examples
#' n <- 50
#' X <- cbind(runif(n), runif(n))
#' rdm(X, method="spearman", bandwidth_method="fixed", bandwidth_parameter=.3)
#' n <- 20
#' U <- runif(n)
#' rdm(cbind(U, U), method="spearman", bandwidth_method="cv", bandwidth_parameter=c(0.25, 0.5))
#' @export rdm
rdm <- function (X, method = c("spearman", "kendall", "dss","zeta1", "bkr", "all"),
                 bandwidth_method = c("fixed", "cv", "cvsym"), bandwidth_parameter = 0.5,
                 permutation = FALSE, npermutation = 1000,
                 checkInput = FALSE) {

  n <- dim(X)[1]
  bm <- match.arg(bandwidth_method, c("cv", "cvsym", "fixed"))

  if(checkInput) {
    if(!is.numeric(bandwidth_parameter)) {
      stop("Non-numeric bandwidth")
    }
    if(bm == "cv") {
      if(length(bandwidth_parameter) < 2) {
        stop("Missing bandwidth parameter")
      } else if(bandwidth_parameter[1] <= 0 || bandwidth_parameter[2] > 0.5 || bandwidth_parameter[2] < bandwidth_parameter[1]) {
        stop("Bandwidth must lie between (0, 0.5) and fulfil b_1 <= b_2")
      }
    } else if(bm == "cvsym"){
      if(bandwidth_parameter[1] <= 0 || bandwidth_parameter[2] > 0.5 || bandwidth_parameter[2] < bandwidth_parameter[1]) {
        stop("Bandwidth must lie between (0, 0.5) and fulfil b_1 <= b_2")
      }
    } else {
      if(bandwidth_parameter <= 0 || bandwidth_parameter > 0.5) {
        stop("Bandwidth must lie between (0, 0.5)")
      }
    }
    if(permutation) {
      compm <- match.arg(method, c("spearman", "kendall", "dss","zeta1", "bkr", "all"))
      if(compm == "all") {
        stop("Permutation test is only possible for one measure at a time")
      }
    }
  }

  #Calculate the bandwidth given the bandwidth method and parameter
  N <- switch(bm,
    cv=computeBandwidth(X, bandwidth_parameter[1], bandwidth_parameter[2], method="cvasym"),
    cvsym = computeBandwidth(X, bandwidth_parameter[1], bandwidth_parameter[2], method="cvsym"),
    fixed = floor(n^bandwidth_parameter))

  A <- NULL
  if(bandwidth_method == "cv") {
    A <- checkerboardDensity(X[, 1], X[, 2], resolution1=N[1], resolution2 = N[2])
  } else {
    A <- checkerboardDensity(X[, 1], X[, 2], resolution1=N, resolution2 = N)
  }
  AS <- sortDSMatrix(as.matrix(A))

  if(method=="all") {
    result <- NULL
    result$spearman <- computeCBMeasure(AS, method="spearman")
    result$kendall <- computeCBMeasure(AS, method="kendall")
    result$dss <- computeCBMeasure(AS, method="dss")
    result$bkr <- computeCBMeasure(AS, method="bkr")
    result$zeta1 <- computeCBMeasure(AS, method="zeta1")
    return(result)
  }

  if(permutation) {
    result <- NULL
    result$rdm <- computeCBMeasure(AS, method)
    randomPermutation <- rep.int(0, npermutation)
    for(k in 1:npermutation) {
      sample <- cbind(X[, 1], X[sample.int(n), 2])
      randomPermutation[k] <- rdm(sample, method=method, bandwidth_method=bandwidth_method, bandwidth_parameter = bandwidth_parameter, checkInput = FALSE)
    }
    result$pvalue <- sum(randomPermutation >= result$rdm)/npermutation
    return(result)
  }

  return(computeCBMeasure(AS, method))
}
