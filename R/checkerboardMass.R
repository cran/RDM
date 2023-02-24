#' @title Estimate the checkerboard mass density
#' @description Estimate a non-square checkerboard mass density
#' @details This implementation modifies the code of build_checkerboard_weights() published in 'qad', version 1.0.4, available at \href{https://CRAN.R-project.org/package=qad}{https://CRAN.R-project.org/package=qad},
#' to allow for non-square checkerboard mass densities.
#' For more details on the implementation see \code{\link[qad]{ECBC}} and for more information on the implemented changes, see the file 'src/code.cpp'.
#' @param X First coordinate of the observations.
#' @param Y Second coordinate of the observations.
#' @param resolution1 A natural number specifying the resolution of the first component.
#' @param resolution2 A natural number specifying the resolution of the second component.
#' @return The estimated checkerboard mass density.
#' @export checkerboardDensity
#' @useDynLib RDM, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @examples
#' checkerboardDensity(runif(20), runif(20), 3, 3)
checkerboardDensity <- function(X, Y, resolution1, resolution2) {
  X <- rank(X, ties.method = "max")
  Y <- rank(Y, ties.method = "max")
  return(round(asymmetric_checkerboard_mass(X, Y, resolution1, resolution2),15))
}


#' @title Estimate a single entry of the checkerboard mass density
#' @description Estimate the value \eqn{A_{kl}} of the non-square checkerboard mass density.
#' @details This implementation modifies the code of build_checkerboard_weights() published in 'qad', version 1.0.4, available at \href{https://CRAN.R-project.org/package=qad}{https://CRAN.R-project.org/package=qad},
#' to allow for the evaluation of a single index of the non-square checkerboard mass densities.
#' For more details on the implementation see \code{\link[qad]{ECBC}} and for more information on the implemented changes, see the file 'src/code.cpp'.
#' @param X First coordinate of the observations.
#' @param Y Second coordinate of the observations.
#' @param k Index of the first component.
#' @param l Index of the second component.
#' @param resolution1 A natural number specifying the resolution of the first component.
#' @param resolution2 A natural number specifying the resolution of the second component.
#' @return The estimated checkerboard mass density \eqn{A_{kl}}.
#' @export checkerboardDensityIndex
#' @examples
#' U <- runif(20)
#' V <- runif(20)
#' checkerboardDensity(U, V, 3, 3)
#' checkerboardDensityIndex(U, V, 1, 2, 3, 3)
checkerboardDensityIndex <- function(X, Y, k, l, resolution1, resolution2) {
  X <- rank(X, ties.method = "max")
  Y <- rank(Y, ties.method = "max")
  return(round(asymmetric_checkerboard_index(X, Y, k, l, resolution1, resolution2), 15))
}

