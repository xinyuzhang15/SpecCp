#include <RcppArmadillo.h>
#include <complex>

using namespace Rcpp;
using namespace arma;

// A Rcpp function to estimate the spectral density matrix at a consecutive period, for the refinement step.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cx_cube bs_consecutive(const arma::mat& ys, const arma::vec& w, const arma::vec& freq,
                         const arma::vec& band, const arma::uvec& band_ind, const arma::uvec& at,
                         const arma::mat& projections, int R, int L, bool coherency) {
  int m = ys.n_cols;
  int nfreq = freq.n_elem;
  int nband = band.n_elem;
  int nat = at.n_elem;
  arma::cube x2(m, m, R + 1, arma::fill::zeros);
  arma::cx_cube bs(1, 1, nband * nat, arma::fill::zeros);

  std::complex<double> i(0.0, 1.0); // Define the imaginary unit once

  for (int ttt = 0; ttt < nat; ++ttt) {
    arma::cx_cube eval(m, m, nfreq, arma::fill::zeros);
    arma::cx_cube eval_band(m, m, nband, arma::fill::zeros);
    int tt = at[ttt] - 1;

    for (int h = 0; h <= R; ++h) {
      arma::uvec ind_l = arma::regspace<arma::uvec>(tt - L + 1 + h, tt);
      if (ttt == 0) { // the start point
        x2.slice(h) = ys.rows(ind_l - h).t() * ys.rows(ind_l);
      } else { // iteration, minus the first one, and add the new one
        x2.slice(h) = x2.slice(h) - ys.row(ind_l[0] - 1 - h).t() * ys.row(ind_l[0] - 1) +
          ys.row(ind_l[L - h - 1] - h).t() * ys.row(ind_l[L - h - 1]);
      }

      for (int f = 0; f < nfreq; ++f) {
        std::complex<double> exp_term = std::exp(-i * static_cast<double>(h) * freq[f]);
        eval.slice(f) += x2.slice(h) * (w[h + R] * exp_term);
        if (h > 0) {
          eval.slice(f) += x2.slice(h).t() * (w[h + R] * std::conj(exp_term));
        }
      }
    }
    eval /= (2 * M_PI * L);

    if (coherency) {
      for (int ii = 0; ii < nfreq; ++ii) {
        arma::cx_vec diag_vals = eval.slice(ii).diag();
        arma::vec Is = arma::sqrt(1 / arma::real(diag_vals));
        arma::mat Is_diag = arma::diagmat(Is);
        eval.slice(ii) = Is_diag * eval.slice(ii) * Is_diag;
        eval.slice(ii).diag().ones();
      }
    }

    for (int i = 0; i < nband; ++i) {
      arma::uvec indices = arma::find(band_ind == (i + 1));
      eval_band.slice(i) = (band[i] == 1) ? eval.slice(indices(0)) : arma::mean(eval.slices(indices), 2);
      bs.slice(ttt * nband + i) = projections.row(i) * eval_band.slice(i) * projections.row(i).t();
    }
  }

  return bs;
}

