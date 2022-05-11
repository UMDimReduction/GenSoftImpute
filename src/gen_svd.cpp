# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]

// [[Rcpp::export]]
Rcpp::List gen_svd_for_R(const Rcpp::NumericMatrix &X,
                   const Rcpp::NumericMatrix &M,
                   const Rcpp::NumericMatrix &W) {
  // Convert to arma::mat
  const arma::mat Xmat = Rcpp::as<arma::mat>(X);
  const arma::mat Mmat = Rcpp::as<arma::mat>(M);
  const arma::mat Wmat = Rcpp::as<arma::mat>(W);

  // Compute Cholesky decompositions
  const arma::mat Mmat_sqrt = arma::chol(Mmat);
  const arma::mat Wmat_sqrt = arma::chol(Wmat);

  // Generalized SVD
  const arma::mat Xtilde = Mmat_sqrt * Xmat * Wmat_sqrt;

  arma::mat U;
  arma::vec s;
  arma::mat V;

  arma::svd(U, s, V, Xtilde);

  return Rcpp::List::create(Rcpp::Named("svalues") = s.as_col(),
                            Rcpp::Named("left_svectors") = inv(Mmat_sqrt) * U,
                            Rcpp::Named("right_svectors")  = inv(Wmat_sqrt) * V);
}

void gen_svd(arma::mat &U, arma::vec &s, arma::mat &V,
             const arma::mat &X,
             const arma::mat &M,
             const arma::mat &W) {
  // Compute Cholesky decompositions
  const arma::mat Mmat_sqrt = arma::chol(M);
  const arma::mat Wmat_sqrt = arma::chol(W);

  // Generalized SVD
  const arma::mat Xtilde = Mmat_sqrt * X * Wmat_sqrt;

  arma::svd(U, s, V, Xtilde);

  s = s.as_col();
  U = inv(Mmat_sqrt) * U;
  V = inv(Wmat_sqrt) * V;
}

void gen_svd(arma::mat &U, arma::vec &s, arma::mat &V,
             const arma::mat &X,
             const arma::vec &M,
             const arma::vec &W) {
  // Compute Cholesky decompositions
  const arma::mat Mmat_sqrt = arma::diagmat(arma::sqrt(M));
  const arma::mat Wmat_sqrt = arma::diagmat(arma::sqrt(W));

  // Generalized SVD
  const arma::mat Xtilde = Mmat_sqrt * X * Wmat_sqrt;

  arma::svd(U, s, V, Xtilde);

  s = s.as_col();
  U = inv(Mmat_sqrt) * U;
  V = inv(Wmat_sqrt) * V;
}
