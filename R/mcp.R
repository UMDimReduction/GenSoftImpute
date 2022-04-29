#' Multiple Correspondence Analysis
#'
#' The function takes in a data matrix and performs Multiple Correspondence
#' Analysis (MCP) using Generalized Singular Value Decomposition.
#'
#' @param X a data matrix to decompose
#' @param k total number of components to return. If \code{NULL}, the full set
#'   of nonzero components are returned.
#' @param tol default is \code{.Machine$double.eps^0.5}. A tolerance level for
#'   eliminating effectively zero (small variance) singular value components.
#'
#' @return A list with three elements:
#'   \item{d}{A vector containing the singular values of X above the
#'   tolerance threshold.}
#'   \item{u}{Left (rows) generalized singular vectors. Dimensions are
#'   \code{nrow(X)} by k.}
#'   \item{v}{Right (columns) generalized singular vectors. Dimensions are
#'   \code{ncol(X)} by k.}
#' @export
MCP <- function(X, k = NULL, tol = .Machine$double.eps^0.5) {
  # Prepare data
  observed_matrix <- X/sum(X)
  row_frequencies <- rowSums(observed_matrix)
  col_frequencies <- colSums(observed_matrix)
  expected_matrix <- row_frequencies %o% col_frequencies
  deviations_matrix <- observed_matrix - expected_matrix
  # Generalized SVD
  Amat <- deviations_matrix
  Mmat <- diag(1/row_frequencies)
  Wmat <- diag(1/col_frequencies)

  Mmat_sqrt <- diag(1/sqrt(row_frequencies))
  Wmat_sqrt <- diag(1/sqrt(col_frequencies))

  Atilde <- Mmat_sqrt %*% Amat %*% Wmat_sqrt

  Atilde_svd <- svd(Atilde)
  keep <- if (is.null(k)) which(Atilde_svd$d > tol) else k

  Umat <- solve(Mmat_sqrt) %*% Atilde_svd$u[,keep]
  Vmat <- solve(Wmat_sqrt) %*% Atilde_svd$v[,keep]

  return(list(d = Atilde_svd$d[keep],
              u = Umat,
              v = Vmat))
}

#----
#' @export
MCP2 <- function(X, k = NULL, tol = .Machine$double.eps^0.5) {
  # Prepare data
  # observed_matrix <- X/sum(X)
  # row_frequencies <- rowSums(observed_matrix)
  # col_frequencies <- colSums(observed_matrix)
  # expected_matrix <- row_frequencies %o% col_frequencies
  # deviations_matrix <- observed_matrix - expected_matrix
  # # Generalized SVD
  # Amat <- deviations_matrix
  # Mmat <- diag(1/row_frequencies)
  # Wmat <- diag(1/col_frequencies)
  #
  # Mmat_sqrt <- diag(1/sqrt(row_frequencies))
  # Wmat_sqrt <- diag(1/sqrt(col_frequencies))
  #
  # Atilde <- Mmat_sqrt %*% Amat %*% Wmat_sqrt
  #
  # Atilde_svd <- svd(Atilde)
  # keep <- if (is.null(k)) which(Atilde_svd$d > tol) else k
  #
  # Umat <- solve(Mmat_sqrt) %*% Atilde_svd$u[,keep]
  # Vmat <- solve(Wmat_sqrt) %*% Atilde_svd$v[,keep]
  #
  # return(list(d = Atilde_svd$d[keep],
  #             u = Umat,
  #             v = Vmat))
  gen_svd(X, diag(nrow(X)), diag(ncol(X)))
}
