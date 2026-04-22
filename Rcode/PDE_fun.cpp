#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat PDE_update_Rcpp(arma::mat h1,
                          arma::mat h2,
                          arma::mat cmatrix) {
  return h1 * cmatrix + h2 % (cmatrix % cmatrix);
}

