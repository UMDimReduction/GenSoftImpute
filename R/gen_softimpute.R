#' Soft Imputation
#'
#' @param X matrix
#' @param ncp Number of components
#' @param lambdas Tuning parameters
#' @param maxiter Maximum number of iterations
#' @param tol Tolerance level
#' @export
soft_impute <- function(X, ncp = 2L, lambdas = c(10, 1, 0.1),
                        maxiter = 100, tol = 0.001) {
  ncp <- min(c(ncp, dim(X)))
  mask <- (is.na(X)) * 1L
  result <- softimpute_cpp(X, mask, ncp, lambdas, tol, maxiter)
  return(drop(result))
}

#----
#' @export
#' @param M Row constraint matrix
#' @param W Column constraint matrix
#' @rdname soft_impute
generalized_soft_impute <- function(X, M = diag(nrow(X)),
                                    W = diag(ncol(X)),
                                    ncp = 2L,
                                    lambdas = c(10, 1, 0.1),
                                    maxiter = 100, tol = 0.001) {
  ncp <- min(c(ncp, dim(X)))
  mask <- (is.na(X)) * 1L
  result <- gensoftimpute_cpp(X, M, W, mask, ncp,
                              lambdas, tol, maxiter)
  return(drop(result))
}
