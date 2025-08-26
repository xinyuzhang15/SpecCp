#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// A Rcpp function for normalizing vector.

// [[Rcpp::export]]
NumericVector mysd(NumericVector x) {
  arma::vec vec(x.begin(), x.size(), false);  // Use Armadillo's vec (no copy)
  double norm_val = norm(vec, 2);  // Compute L2 norm

  // If the norm is zero, return the original vector
  if (norm_val == 0) {
    return x;
  } else {
    arma::vec standardized_vec = vec / norm_val;

    // Convert to NumericVector and remove names or dim attributes
    NumericVector result = as<NumericVector>(wrap(standardized_vec));
    result.attr("dim") = R_NilValue;  // Remove any dimension attributes
    return result;
  }
}
