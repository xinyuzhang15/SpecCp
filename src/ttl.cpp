#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// A Rcpp version of the rTensor::ttl function, i.e., contracted (m-Mode) product between a Tensor of arbitrary number of modes and a list of matrices. The result is folded back into Tensor.

// Helper function to unfold a 3D tnsr along mode m (1, 2, or 3)
// [[Rcpp::export]]
arma::mat unfold_tnsr(const arma::cube& tnsr, int mode) {
  int I = tnsr.n_rows, J = tnsr.n_cols, K = tnsr.n_slices;

  if (mode == 1) {
    // Mode-1 unfolding: result is I x (J * K)
    arma::mat unfolded(I, J * K);
    for (int k = 0; k < K; ++k) {
      for (int j = 0; j < J; ++j) {
        for (int i = 0; i < I; ++i) {
          unfolded(i, j + k * J) = tnsr(i, j, k);
        }
      }
    }
    return unfolded;
  } else if (mode == 2) {
    // Mode-2 unfolding: result is J x (I * K)
    arma::mat unfolded(J, I * K);
    for (int k = 0; k < K; ++k) {
      for (int j = 0; j < J; ++j) {
        for (int i = 0; i < I; ++i) {
          unfolded(j, i + k * I) = tnsr(i, j, k);
        }
      }
    }
    return unfolded;
  } else if (mode == 3) {
    // Mode-3 unfolding: result is K x (I * J)
    arma::mat unfolded(K, I * J);
    for (int j = 0; j < J; ++j) {
      for (int i = 0; i < I; ++i) {
        for (int k = 0; k < K; ++k) {
          unfolded(k, i + j * I) = tnsr(i, j, k);
        }
      }
    }
    return unfolded;
  } else {
    stop("Invalid mode. Mode must be 1, 2, or 3.");
  }
}

// Helper function to fold a matrix back into a tnsr with updated dimensions
arma::cube fold_tnsr(const arma::mat& unfolded, arma::uvec dims, int mode) {
  int I = dims[0], J = dims[1], K = dims[2];
  arma::cube folded(I, J, K);

  if (mode == 1) {
    // Mode-1 folding
    for (int k = 0; k < K; ++k) {
      for (int j = 0; j < J; ++j) {
        for (int i = 0; i < I; ++i) {
          folded(i, j, k) = unfolded(i, j + k * J);
        }
      }
    }
  } else if (mode == 2) {
    // Mode-2 folding
    for (int k = 0; k < K; ++k) {
      for (int j = 0; j < J; ++j) {
        for (int i = 0; i < I; ++i) {
          folded(i, j, k) = unfolded(j, i + k * I);
        }
      }
    }
  } else if (mode == 3) {
    // Mode-3 folding
    for (int j = 0; j < J; ++j) {
      for (int i = 0; i < I; ++i) {
        for (int k = 0; k < K; ++k) {
          folded(i, j, k) = unfolded(k, i + j * I);
        }
      }
    }
  } else {
    stop("Invalid mode. Mode must be 1, 2, or 3.");
  }

  return folded;
}

// [[Rcpp::export]]
arma::cube ttl_rcpp(arma::cube tnsr, List list_mat, IntegerVector ms) {
  arma::uvec tnsr_dims = {tnsr.n_rows, tnsr.n_cols, tnsr.n_slices};
  int num_mats = list_mat.size();

  if (ms.size() != num_mats) {
    stop("The length of 'ms' and 'list_mat' must match.");
  }

  for (int i = 0; i < num_mats; ++i) {
    int m = ms[i] - 1; // Adjust for 0-based indexing in C++

    // Convert R matrix to Armadillo matrix
    NumericMatrix mat_r = list_mat[i];
    arma::mat mat(mat_r.begin(), mat_r.nrow(), mat_r.ncol(), false, true);

    // Check compatibility of dimensions
    if (tnsr_dims[m] != mat.n_cols) {
      stop("Matrix columns do not match tnsr mode dimension.");
    }

    // Prepare unfolded tnsr
    arma::mat unfolded = unfold_tnsr(tnsr, m + 1);

    // Matrix multiplication
    arma::mat result = mat * unfolded;

    // Update dimensions after multiplication
    tnsr_dims[m] = mat.n_rows;

    // Fold back to tnsr
    tnsr = fold_tnsr(result, tnsr_dims, m + 1);
  }

  return tnsr;
}
