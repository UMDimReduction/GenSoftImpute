# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]

// [[Rcpp::export]]
Rcpp::List gen_svd(const Rcpp::NumericMatrix &X,
                   const Rcpp::NumericMatrix &M,
                   const Rcpp::NumericMatrix &W) {
  // Convert to arma::mat
  arma::mat Xmat = Rcpp::as<arma::mat>(X);
  arma::mat Mmat = Rcpp::as<arma::mat>(M);
  arma::mat Wmat = Rcpp::as<arma::mat>(W);

  // Compute Cholesky decompositions
  arma::mat Mmat_sqrt = arma::chol(Mmat);
  arma::mat Wmat_sqrt = arma::chol(Wmat);

  // Generalized SVD
  arma::mat Xtilde = Mmat_sqrt * Xmat * Wmat_sqrt;

  arma::mat U;
  arma::vec s;
  arma::mat V;

  arma::svd(U, s, V, Xtilde);

  arma::mat Umat = inv(Mmat_sqrt) * U;
  arma::mat Vmat = inv(Wmat_sqrt) * V;

  arma::vec d = s.as_col();
  return Rcpp::List::create(Rcpp::Named("svalues") = d,
                            Rcpp::Named("left_svectors") = Umat,
                            Rcpp::Named("right_svectors")  = Vmat);
}
