#' Soft Imputation
#'
#' @param X matrix
#' @param lambdas Tuning parameters
#' @param maxiter Maximum number of iterations
#' @param tol Tolerance level
#' @export
soft_impute <- function(X, lambdas = c(10, 1, 0.1),
                        maxiter = 100, tol = 0.001) {
  n <- nrow(X)
  p <- ncol(X)
  IDmat <- (is.na(X)) * 1L

  # Minit <- X
  # Minit[is.na(Minit)] <- mean(Minit, na.rm = TRUE)
  result <- gensoftimpute_cpp(X, IDmat, lambdas, tol, maxiter)
  return(result)
}
