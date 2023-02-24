#' @title Sort a (possibly non-square) doubly stochastic matrix
#' @description Sorts an arbitrary doubly stochastic \eqn{N_1 \times N_2} matrix A into the matrix \eqn{A^\uparrow} such that the induced checkerboard copula \eqn{C(A^\uparrow)} is stochastically increasing.
#' @details The algorithm to sort a doubly stochastic matrix \eqn{A} is given in Strothmann, Dette and Siburg (2022) <arXiv:2201.03329>.
#' Since this implementation does not depend on the appropriate scaling of the matrix \eqn{A}, both doubly stochastic matrices and checkerboard mass densities are admissible inputs.
#' @param A A (possibly non-square) doubly stochastic matrix or (possibly non-square) checkerboard mass density.
#' @return The sorted version \eqn{A^\uparrow} of the matrix \eqn{A}.
#' @importFrom Rfast "colCumSums" "colSort" "coldiffs"
#' @export sortDSMatrix
#' @examples
#' n <- 4
#' A <- diag(n)[n:1, ]
#' print(A)
#' sortDSMatrix(A)
sortDSMatrix <- function(A) {

  N1 <- dim(A)[1]
  N2 <- dim(A)[2]
  #In case either dimension is less or equal to 1, there is no rearrangement necessary
  if(min(N1, N2) <= 1) {
    return(A)
  }

  cumSumA <- t(colCumSums(t(A)))
  #Sort each column of cumulative sums in decreasing order
  sortedA <- colSort(cumSumA, descending = TRUE)
  sortedA[, 2:N2] <- coldiffs(sortedA)
  return(sortedA)
}
