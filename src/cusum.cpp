#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Rcpp functions for computing CUSUM for general vectors and symmetric tensors (the latter is specific for our spectral matrix case)

// [[Rcpp::export]]
arma::cube cusum_symm_tensor(const arma::cube& z) {
  int n = z.n_rows;
  int len = z.n_slices;
  cube cs(n, n, len - 1, fill::zeros);
  double sqrtlen = std::sqrt(len);
  double factor = std::sqrt(len - 1) / sqrtlen;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) { // Only compute for lower triangle and diagonal
      double iplus = z(i, j, 0);
      double iminus = accu(z.tube(i, j)) - z(i, j, 0); // Use accu to sum the elements

      cs(i, j, 0) = factor * (iplus - iminus / (len - 1));

      if (len - 2 > 0) {
        for (int k = 1; k < len - 1; ++k) {
          factor = std::sqrt(k + 1) * std::sqrt(len - k - 1) / sqrtlen;
          iplus += z(i, j, k);
          iminus -= z(i, j, k);
          cs(i, j, k) = factor * (iplus / (k + 1) - iminus / (len - k - 1));
        }
      }

      // Fill the upper triangle using symmetry for all slices
      if (i != j) {
        for (int k = 0; k < len - 1; ++k) {
          cs(j, i, k) = cs(i, j, k);
        }
      }
    }
  }

  return cs;
}

// [[Rcpp::export]]
NumericVector cusum_vec(const arma::vec& z) {
  int len = z.n_elem;
  arma::vec cs(len - 1, fill::zeros);
  double sqrtlen = std::sqrt(len);
  double iplus = z(0);
  double iminus = accu(z) - z(0);
  double factor = std::sqrt(len - 1) / sqrtlen;

  cs(0) = factor * (iplus - iminus / (len - 1));

  if (len - 2 > 0) {
    for (int k = 1; k < len - 1; ++k) {
      factor = std::sqrt(k + 1) * std::sqrt(len - k - 1) / sqrtlen;
      iplus += z(k);
      iminus -= z(k);
      cs(k) = factor * (iplus / (k + 1) - iminus / (len - k - 1));
    }
  }

  return NumericVector(cs.begin(), cs.end());
}

