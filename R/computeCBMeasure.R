#' @title Dependence measures for the checkerboard copula
#' @description Computes \eqn{\mu(C^{\#}(A))} for some underlying measure for the checkerboard copula \eqn{C^{\#}(A)}.
#' This measure depends only on the input matrix A.
#' @details This function computes \eqn{\mu(C^{\#}(A))} for one of several underlying measures for a given checkerboard copula \eqn{C^{\#}(A)}.
#' Most importantly, the value only depends on the (possibly non-square) matrix \eqn{A} and implicitly assumes the form of \eqn{C^{\#}(A)} given in Strothmann, Dette and Siburg (2022) <arXiv:2201.03329>.
#' Currently, the following underlying measures are implemented:
#' \itemize{
#'  \item{"spearman"}{ Implements the concordance measure Spearman's \eqn{\rho},}
#'  \item{"kendall"}{ Implements the concordance measure Kendall's \eqn{\tau},}
#'  \item{"bkr"}{ Implements the Blum–Kiefer–Rosenblatt \eqn{R}, also known as the \eqn{L^2}-Schweizer-Wolff-measure <doi:10.1214/aos/1176345528>,}
#'  \item{"dss"}{ Implements the Dette-Siburg-Stoimenov measure of complete dependence <doi:10.1111/j.1467-9469.2011.00767.x>, also known as Chatterjee's \eqn{\xi} <doi:10.1080/01621459.2020.1758115>,}
#'  \item{"zeta1"}{ Implements the \eqn{\zeta_1}-measure of complete dependence established by W. Trutschnig <doi:10.1016/j.jmaa.2011.06.013>.}
#' }
#' @param A A (possibly non-square) checkerboard mass density.
#' @param method Determines the underlying dependence measure. Options include "spearman", "kendall", "bkr", "dss", "chatterjee" and "zeta1".
#' @return The value of \eqn{\mu(C^{\#}(A))}. For a sorted A, this corresponds to the rearranged dependence measure \eqn{R_{\mu}(C^{\#}(A))}.
#' @export computeCBMeasure
#' @importFrom Rfast "colCumSums"
#' @examples
#' n <- 10
#' A <- diag(n)/n
#' computeCBMeasure(A, method="spearman")

computeCBMeasure <- function(A, method = c("spearman", "kendall", "bkr", "dss", "zeta1")) {
  N1 <- dim(A)[1]
  N2 <- dim(A)[2]
  #In case N1 or N2 are less than 2, the value is necessarily zero for all admissable measures mu
  if(min(N1, N2) <= 1) {
    return(0)
  }

  method <- match.arg(method, c("spearman", "kendall", "dss", "chatterjee", "zeta1", "bkr"))

  if(method == "kendall") {
    T1 <- rbind(rep.int(0, N2), colCumSums(A)[1:(N1-1), ])
    T2 <- t(rbind(rep.int(0, N1), colCumSums(t(A))[1:(N2-1), ]))
    Q <- sum(T1*T2 + 0.5*A*(T1+T2 + 0.5*A))
    return(1 - 4*Q)
  } else if(method == "dss" || method == "chatterjee") {
    J = matrix(rep(1:N2, N2), nrow=N2, ncol=N2)
    X = ((N2+0.5)*matrix(1, nrow=N2, ncol=N2) - pmax(t(J), J) - diag(N2)/6)
    W = t(A) %*% A;
    return(6*N1*sum(W*X)/N2-2)
  } else if(method == "bkr") {
    B <- (A - 1/(N1*N2))
    J = matrix(rep(1:N1, N1), nrow=N1, ncol=N1)
    P1 = ((N1+0.5)*matrix(1, nrow=N1, ncol=N1) - pmax(t(J), J) - diag(N1)/6)
    J = matrix(rep(1:N2, N2), nrow=N2, ncol=N2)
    P2 = ((N2+0.5)*matrix(1, nrow=N2, ncol=N2) - pmax(t(J), J) - diag(N2)/6)
    return(sqrt(90*sum(B*(P1 %*% B %*% P2))/(N1*N2)))
  } else if(method == "zeta1") {
    B <- (N1*N2*A-1)/N1
    K <- cbind(0, t(colCumSums(t(B)))[, 1:(N2-1)])/N2
    B[abs(B) < 1e-16] <- 0
    H <- K/B
    H[!is.finite(H)] <- 0
    result <- sum(abs(K[B == 0])/N2)
    idx <- which(H > 0)
    result <- result + sum(abs(B[idx]/(2*N2^2) + K[idx]/N2))
    idx <- which(H <= -1/N2)
    result <- result + sum(abs(B[idx]/(2*N2^2) + K[idx]/N2))
    idx <- which(((B > 0 & K <= 0) | (B < 0 & K >= 0)) & H > -1/N2)
    result <- result + sum(abs(K[idx]^2/B[idx] + B[idx]/(2*N2^2) + K[idx]/N2))
    return(3*(result))
  }
  #Spearmans rho coincide
  T1 <- (2*N1 - 2*matrix(rep(1:N1, N2), nrow=N1, ncol=N2)+1)/(2*N1)
  T2 <- (2*N2 - 2*matrix(rep(1:N2, N1), nrow=N2, ncol=N1)+1)/(2*N2)
  W <- T1*t(T2)

  return(12*sum(A*W)-3)
}


