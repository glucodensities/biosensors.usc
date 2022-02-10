// ConfidenceBand.h: biosensors.usc glue
//
// Copyright (C) 2019 - 2021  Juan C. Vidal and Marcos Matabuena
//
// This file is part of biosensors.usc.
//
// biosensors.usc is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// biosensors.usc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with biosensors.usc.  If not, see <http://www.gnu.org/licenses/>.

#ifndef _CONFIDENCE_BAND_H // include guard
#define _CONFIDENCE_BAND_H

#include <stdlib.h>
#include <RcppArmadillo.h>
#include "WassersteinRegression.h"


namespace bio {

struct confidence_struct {
  arma::mat Qpred;
  arma::mat Q_lx;
  arma::mat Q_ux;
  arma::mat fpred;
};

inline arma::uword sumNumbers(arma::mat x) {
  arma::uword value = 0;
  arma::uword n = x.n_rows;
  arma::uword m = x.n_cols;
  for (arma::uword i=0; i < n; i++) {
    for (arma::uword j=0; j < m; j++) {
      if (!std::isnan(x(i,j))) {
        value++;
      }
    }
  }
  return value;
}

inline double quantile(arma::vec x, double p) {
  arma::mat y = arma::sort(x);
  arma::uword m = sumNumbers(y); // check for NaN
  double pp = p * m + 0.5;
  arma::uword pi = std::max(std::min((arma::uword)std::floor(pp), m-1), (arma::uword)1);
  double pr = std::max(std::min (pp - pi, 1.0), 0.0);
  return (1-pr) * y(pi-1) + pr * y(pi);
}


/**
 * This function computes intrinsic confidence bands for Wasserstein regression.
 * Inputs:
 *     xfit  - nxp matrix of predictor values for fitting (do not include a column for the intercept)
 *     xpred - kxp vector of input values for regressors for prediction.
 *     Q_obs - nxm matrix of quantile functions. Q_obs(i, :) is a 1xm vector of quantile function values on grid t_vec.
 *     q_obs - nxm matrix of quantile density functions. q_obs(i, :) is a 1xm vector of quantile density function values on grid t_vec.
 *     t_vec - 1xm vector - common grid for all quantile density functions in Q_obs, q_obs, q_prime_obs.
 *     alpha - 100*(1 - alpha) is the significant level
 *     delta - in (0, 1/2), the boundary control value
 * Outputs:
 *   A structure with the following fields:
 *     Q_lx  - lower bound of confidence bands in terms of density functions
 *     Q_ux  - upper bound of confidence bands in terms of density functions
 *     Qpred - fitted density function at xpred.
 */
inline confidence_struct confidence_band(const arma::mat xfit, const arma::mat xpred, const arma::mat Q_obs, const arma::mat q_obs, const arma::vec t_vec, const double alpha) {
  arma::uword n = xfit.n_rows;
  arma::uword k = xpred.n_rows;
  arma::uword m = t_vec.n_elem;


  // ===============  1) compute fitted values  ================== //
  regression_struct res = wasserstein_regression(xfit, q_obs, Q_obs.col(0), xpred, t_vec, 1e-6);
  arma::mat fpred = res.fpred;
  arma::mat Qfit = res.Qfit;
  arma::mat Qpred = res.Qpred;

  // ===============  2) compute  D(s, t)  ======================= //
  arma::mat Xmat = arma::join_horiz(arma::ones(n,1), xfit);
  arma::mat Sigma = Xmat.t() * Xmat / n;
  arma::mat Q_res = Q_obs - Qfit;
  arma::field<arma::mat> D_cell(m,m);
  for (arma::uword i=0; i < m; i++) {
    // to save time, compute lower triangle C_cell
    for (arma::uword j=0; j <= i; j++) {
      D_cell(i,j) = Xmat.t() * arma::diagmat(Q_res.col(i) % Q_res.col(j)) * Xmat / n;
      D_cell(j,i) = D_cell(i,j);
    }
  }

  int R = 1000;
  // arma::vec m_alpha(k); m_alpha.zeros();
  arma::vec m_alpha = arma::zeros(k);
  arma::mat se = arma::zeros(k, m);
  arma::mat C_x = arma::zeros(m, m);

  for (arma::uword l=0; l < k; l++) {
    arma::vec x_star = solve(Sigma, arma::join_horiz(arma::ones(1,1), xpred.row(l)).t());
    for (arma::uword i=0; i < m; i++) {
      for (arma::uword j=0; j <= i; j++) {
        C_x(i,j) = arma::conv_to<double>::from(x_star.t() * D_cell(i,j) * x_star);
        C_x(j,i) = C_x(i,j);
      }
    }
    // Compute m_alpha
    arma::mat aux = arma::diagmat(C_x);
    arma::mat C_x_diag = sqrt(aux.diag());

    // Compute eigenfunctions of R_x
    arma::vec eigValues;
    arma::mat eigFuns;
    eig_sym(eigValues, eigFuns, C_x);

    // Note: discard the negative eigenvalues and corresponding eigenvectors
    arma::uvec logi_index = find(eigValues > 0);

    eigValues = eigValues(logi_index);
    eigFuns = eigFuns.cols(logi_index);

    // Compute m_alpha
    arma::uvec index_robust = find(eigValues > 0.001 * sum(eigValues));
    eigValues = eigValues(index_robust);
    eigFuns = eigFuns.cols(index_robust);

    // Generate independent normal variable / FPC scores
    arma::uword dim_gau = eigValues.n_elem;
    arma::mat FPCs = arma::randn(R, dim_gau) * arma::diagmat(sqrt(eigValues));

    // Nsimu number of Gaussian Processes
    arma::mat GaussinProcess = FPCs * eigFuns.t();     // R x m matrix
    //
    // Get maximum of Gaussian Processes
    arma::mat sequence_max = max((abs(GaussinProcess) * arma::diagmat(pow(C_x_diag,-1))).t()); // 1 x Nsimu vector contains maximum value of each GP

    // get 1-alpha percentile in the maximum sequence
    m_alpha(l) = quantile(arma::conv_to<arma::vec>::from(sequence_max), 1-alpha);

    se.row(l) = arma::conv_to<arma::rowvec>::from(C_x_diag * 1/sqrt(n));
  }

  // ==================   3) compute Q_lx and Q_ux      ================== %%
  arma::mat Q_lx = arma::zeros(k, m);
  arma::mat Q_ux = arma::zeros(k, m);

  arma::mat H = arma::diagmat(arma::ones(1,m));

  arma::mat A_subtrahend = arma::join_horiz(arma::zeros(m,1), H.cols(0, H.n_cols-2));
  arma::mat A = H - A_subtrahend;
  A = A.rows(0, A.n_cols-2);
  arma::vec m_zeros(m, arma::fill::zeros);

  for (arma::uword i=0; i < k; i++) {
    arma::rowvec Qi_lx = Qpred.row(i) - m_alpha(i) * se.row(i);
    arma::rowvec Qi_ux = Qpred.row(i) + m_alpha(i) * se.row(i);

    arma::uvec aux = find(arma::diff(Qi_lx) < 0);
    if (aux.n_elem > 0) {
      arma::vec b = arma::diff(Qpred.row(i)).t();
      arma::vec lb = -1 * m_alpha(i) * se.row(i).t();
      arma::vec delta_Q = quadprog(H, -lb, -A, b, lb, m_zeros);
      arma::rowvec Qi_lx = Qpred.row(i) + delta_Q.t();
    }

    aux = find(arma::diff(Qi_ux) < 0);
    if (aux.n_elem > 0) {
      arma::vec b = arma::diff(Qpred.row(i)).t();
      arma::vec ub = m_alpha(i) * se.row(i).t();
      arma::vec delta_Q = quadprog(H, -ub, -A, b, m_zeros, ub);
      // arma::vec delta_Q = quadprog(DENSE_AUL, H, -ub, -A, b, m_zeros, ub);
      arma::rowvec Qi_ux = Qpred.row(i) + delta_Q.t();
    }

    Q_lx.row(i) = Qi_lx;
    Q_ux.row(i) = Qi_ux;
  }

  confidence_struct result;
  result.Qpred = Qpred;
  result.Q_lx = Q_lx;
  result.Q_ux = Q_ux;
  result.fpred = fpred;
  return result;
}



}

#endif
