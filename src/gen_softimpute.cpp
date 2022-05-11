# include <RcppArmadillo.h>
# include "gen_svd.h"
// [[ Rcpp :: depends ( RcppArmadillo )]]

void soft_thresholding(arma::vec &s, const double lambda){
  for (int i = 0; i < s.n_elem; i++) {
    if (s(i) > lambda){
      s(i) -= lambda;
    } else {
      s(i) = 0.0;
    }
  }
}

// [[Rcpp::export]]
arma::cube softimpute_cpp(arma::mat &X, const arma::Mat<unsigned short int> &mask,
                          const arma::vec &lambdas, const double tol,
                          const int maxiter) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  const int K = lambdas.n_elem;

  arma::cube output(n, p, K, arma::fill::zeros);
  arma::mat Z_current(n, p, arma::fill::zeros);
  arma::mat Z_next(n, p, arma::fill::zeros);

  arma::uvec id_mis = find(mask == 1);

  arma::mat U;
  arma::vec s;
  arma::mat V;

  for (int k = 0; k < K; k++) {
    double lambdak = lambdas(k);
    double diff = 1.0;
    int iter = 0;

    while (diff > tol && iter < maxiter){
      iter++;
      X(id_mis) = Z_current(id_mis);

      svd(U, s, V, X);
      soft_thresholding(s, lambdak);
      Z_next = U * arma::diagmat(s) * V.t();

      diff = std::pow(arma::norm(Z_current - Z_next,
                                 "fro")/arma::norm(Z_next, "fro"), 2);

      Z_current = Z_next;
    }

    output.slice(k) = Z_current;
  }

  return(output);
}

// [[Rcpp::export]]
arma::cube gensoftimpute_cpp(arma::mat &X, const arma::mat &M, const arma::mat &W,
                             const arma::Mat<unsigned short int> &mask,
                             const arma::vec &lambdas, const double tol,
                             const int maxiter) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  const int K = lambdas.n_elem;

  arma::cube output(n, p, K, arma::fill::zeros);
  arma::mat Z_current(n, p, arma::fill::zeros);
  arma::mat Z_next(n, p, arma::fill::zeros);

  arma::uvec id_mis = find(mask == 1);

  arma::mat U;
  arma::vec s;
  arma::mat V;

  for (int k = 0; k < K; k++) {
    double lambdak = lambdas(k);
    double diff = 1.0;
    int iter = 0;

    while (diff > tol && iter < maxiter){
      iter++;
      X(id_mis) = Z_current(id_mis);

      gen_svd(U, s, V, X, M, W);
      soft_thresholding(s, lambdak);
      Z_next = U * arma::diagmat(s) * V.t();

      diff = std::pow(arma::norm(Z_current - Z_next,
                                 "fro")/arma::norm(Z_next, "fro"), 2);

      Z_current = Z_next;
    }

    output.slice(k) = Z_current;
  }

  return(output);
}
