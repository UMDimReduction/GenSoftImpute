#ifndef GENSVD_H
#define GENSVD_H

Rcpp::List gen_svd_for_R(const Rcpp::NumericMatrix &X,
                         const Rcpp::NumericMatrix &M,
                         const Rcpp::NumericMatrix &W);

void gen_svd(arma::mat &U, arma::vec &s, arma::mat &V,
             const arma::mat &X,
             const arma::mat &M,
             const arma::mat &W);

#endif
