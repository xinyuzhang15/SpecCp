#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// A Rcpp function for matrix power update.

// [[Rcpp::export]]
List powerIteration0(const arma::mat& X, int k, Nullable<NumericVector> v1_init = R_NilValue, int max_iter = 10, double lambda_diff_threshold = 1e-06, bool trace = false) {
  int p = X.n_cols;
  arma::vec v1;
  if (v1_init.isNull()) {
    v1 = arma::vec(p, arma::fill::ones) / std::sqrt(p);
  } else {
    v1 = as<arma::vec>(v1_init);
  }
  if (k < 1 || k > p) {
    stop("If specified, k must be between 1 and %d", p);
  }
  double lambda1 = NA_REAL;
  for (int i = 0; i < max_iter; ++i) {
    if (trace) {
      Rcout << "Iteration " << i + 1 << std::endl;
    }
    v1 = X * v1;
    v1 /= arma::norm(v1, 2);
    if (k != p) {
      arma::uvec indices = arma::sort_index(arma::abs(v1), "ascend");
      v1.elem(indices.head(p - k)).zeros();
      v1 /= arma::norm(v1, 2);
    }
    double lambda1_est = arma::as_scalar(v1.t() * X * v1);
    if (i > 0 && std::abs(lambda1_est - lambda1) <= lambda_diff_threshold * std::abs(lambda1_est)) {
      if (trace) {
        Rcout << "Lambda1 diff less than threshold at iteration " << i + 1 << ", exiting." << std::endl;
      }
      break;
    }
    lambda1 = lambda1_est;
  }
  return List::create(Named("v1") = wrap(v1), Named("lambda1") = lambda1, Named("num_iter") = max_iter);
}

