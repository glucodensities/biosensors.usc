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

  // std::cout << "Q_obs: \n" << Q_obs << std::endl;
  // std::cout << "q_obs: \n" << q_obs << std::endl;

  // ===============  1) compute fitted values  ================== //
  regression_struct res = wasserstein_regression(xfit, q_obs, Q_obs.col(0), xpred, t_vec, 1e-6);
  arma::mat fpred = res.fpred;
  arma::mat Qfit = res.Qfit;
  arma::mat Qpred = res.Qpred;
  // std::cout << "fpred: \n" << fpred << std::endl;
  // std::cout << "Qfit: \n" << Qfit << std::endl;
  // std::cout << "Qpred: \n" << Qpred << std::endl;

  // ==== INICIO BORRAR - Matrices para comprobar el resto del método
  // arma::mat fpred = {{0.0000,  0.0000,  0.0000,  0.0000,  6.5844,  0.0000,  0.0005,  7.6998, 9.0053, 0.0000, 0.0000},
  // {0.0000,  0.0000,  0.0000, 10.0000,  9.9995,  9.9959,  9.9937,  9.9999, 0.0000, 0.0000, 0.0000},
  // {0.0000,  0.0000,  0.0000,  9.9783,  9.9625,  0.0000,  0.1643,  9.9275, 0.0001, 0.0000, 0.0000},
  // {0.0000,  0.0000,  0.0000, 10.0000, 10.0000,  0.0000,  9.9945, 10.0000, 9.9997, 0.0000, 0.0000},
  // {0.0000,  0.0000, 10.0000, 10.0000, 10.0000, 10.0000, 10.0000, 10.0000, 0.0001, 0.0000, 0.0000}};
  //
  // arma::mat Qfit = {{-1.5956, -1.4293, -1.3148, -1.2294, -1.1893, -1.0604, -0.9314, -0.9314, -0.9314, -0.8894, -0.6585},
  // {-1.6540, -1.1768, -0.9757, -0.9757, -0.9757, -0.9757, -0.9757, -0.9757, -0.9671, -0.9302, -0.5962},
  // {-0.7262, -0.6922, -0.6245, -0.5569, -0.4872, -0.4252, -0.3731, -0.3199, -0.2688, -0.2261, -0.1722}};
  //
  // arma::mat Qpred = {{-1.1959, -1.0122, -0.8900, -0.8264, -0.8139, -0.7763, -0.7377, -0.7367, -0.7367, -0.6895, -0.4372},
  // {-1.2263, -0.9780, -0.8222, -0.7866, -0.7866, -0.7866, -0.7866, -0.7866, -0.7765, -0.7368, -0.4802},
  // {-1.4302, -1.1768, -0.9900, -0.9276, -0.9276, -0.8982, -0.8688, -0.8688, -0.8643, -0.8187, -0.5468},
  // {-1.8766, -1.5210, -1.3193, -1.2687, -1.2687, -1.2057, -1.1426, -1.1426, -1.1426, -1.1041, -0.7976},
  // {-1.3553, -1.0462, -0.9087, -0.9087, -0.9087, -0.9087, -0.9087, -0.9087, -0.9010, -0.8816, -0.6229}};
  // ==== FIN BORRAR


  // ===============  2) compute  D(s, t)  ======================= //
  arma::mat Xmat = arma::join_horiz(arma::ones(n,1), xfit);
  // std::cout << "Xmat: \n" << Xmat << std::endl;
  arma::mat Sigma = Xmat.t() * Xmat / n;
  // std::cout << "Sigma: \n" << Sigma << std::endl;
  arma::mat Q_res = Q_obs - Qfit;
  // std::cout << "Q_res: \n" << Q_res << std::endl;
  arma::field<arma::mat> D_cell(m,m);
  for (arma::uword i=0; i < m; i++) {
    // to save time, compute lower triangle C_cell
    for (arma::uword j=0; j <= i; j++) {
      D_cell(i,j) = Xmat.t() * arma::diagmat(Q_res.col(i) % Q_res.col(j)) * Xmat / n;
      D_cell(j,i) = D_cell(i,j);
    }
  }
  // std::cout << "D_cell(0,0): \n" << D_cell(0,0) << std::endl;

  // for l = 1:k  3.1) compute m_alpha and se
  int R = 1000;
  // arma::vec m_alpha(k); m_alpha.zeros();
  arma::vec m_alpha = arma::zeros(k);
  arma::mat se = arma::zeros(k, m);
  arma::mat C_x = arma::zeros(m, m);

  // ==== INICIO BORRAR
  // k = 1;
  // ==== FIN BORRAR
  for (arma::uword l=0; l < k; l++) {
    arma::vec x_star = solve(Sigma, arma::join_horiz(arma::ones(1,1), xpred.row(l)).t());
    // ==== INICIO BORRAR
    // arma::vec x_star = {-0.2660, 0.6969, 0.1451, -0.2437, -1.6740};
    // ==== FIN BORRAR

    // std::cout << "x_star: \n" << x_star << std::endl;
    for (arma::uword i=0; i < m; i++) {
      for (arma::uword j=0; j <= i; j++) {
        C_x(i,j) = arma::conv_to<double>::from(x_star.t() * D_cell(i,j) * x_star);
        C_x(j,i) = C_x(i,j);
      }
    }
    // std::cout << "C_x: \n" << C_x << std::endl;
    // Compute m_alpha
    arma::mat aux = arma::diagmat(C_x);
    arma::mat C_x_diag = sqrt(aux.diag());
    // std::cout << "C_x_diag: \n" << C_x_diag << std::endl;

    // Compute eigenfunctions of R_x
    arma::vec eigValues;
    arma::mat eigFuns;
    eig_sym(eigValues, eigFuns, C_x);
    // std::cout << "eigValues: \n" << eigValues << std::endl;
    // std::cout << "eigFuns: \n" << eigFuns << std::endl;

    // Note: discard the negative eigenvalues and corresponding eigenvectors
    arma::uvec logi_index = find(eigValues > 0);
    // std::cout << "logi_index: \n" << logi_index << std::endl;

    eigValues = eigValues(logi_index);
    eigFuns = eigFuns.cols(logi_index);
    // std::cout << "eigValues 2: \n" << eigValues << std::endl;
    // std::cout << "eigFuns 2: \n" << eigFuns << std::endl;

    // Compute m_alpha
    arma::uvec index_robust = find(eigValues > 0.001 * sum(eigValues));
    eigValues = eigValues(index_robust);
    eigFuns = eigFuns.cols(index_robust);
    // std::cout << "eigValues 3: \n" << eigValues << std::endl;
    // std::cout << "eigFuns 3: \n" << eigFuns << std::endl;

    // std::cout << "eigValues size: " << eigValues.n_elem << std::endl;
    // std::cout << "eigFuns size: " << eigFuns.n_rows << ", " << eigFuns.n_cols << std::endl;

    // Generate independent normal variable / FPC scores
    arma::uword dim_gau = eigValues.n_elem;
    arma::mat FPCs = arma::randn(R, dim_gau) * arma::diagmat(sqrt(eigValues));
    // std::cout << "diagmat(sqrt(eigValues)): \n" << arma::diagmat(sqrt(eigValues)) << std::endl;
    // std::cout << "FPCs size: " << FPCs.n_rows << ", " << FPCs.n_cols << std::endl;

    // Nsimu number of Gaussian Processes
    arma::mat GaussinProcess = FPCs * eigFuns.t();     // R x m matrix
    // std::cout << "GaussinProcess: " << GaussinProcess.n_rows << ", " << GaussinProcess.n_cols << std::endl;
    //
    // Get maximum of Gaussian Processes
    arma::mat sequence_max = max((abs(GaussinProcess) * arma::diagmat(pow(C_x_diag,-1))).t()); // 1 x Nsimu vector contains maximum value of each GP
    // std::cout << "sequence_max: " << sequence_max.n_rows << ", " << sequence_max.n_cols << std::endl;

    // get 1-alpha percentile in the maximum sequence
    m_alpha(l) = quantile(arma::conv_to<arma::vec>::from(sequence_max), 1-alpha);
    // std::cout << "m_alpha: " << m_alpha << std::endl;

    se.row(l) = arma::conv_to<arma::rowvec>::from(C_x_diag * 1/sqrt(n));
    // std::cout << "se: " << se.row(l) << std::endl;
  }

  // ==================   3) compute Q_lx and Q_ux      ================== %%
  arma::mat Q_lx = arma::zeros(k, m);
  arma::mat Q_ux = arma::zeros(k, m);

  arma::mat H = arma::diagmat(arma::ones(1,m));
  // std::cout << "H: \n" << H << std::endl;

  arma::mat A_subtrahend = arma::join_horiz(arma::zeros(m,1), H.cols(0, H.n_cols-2));
  // std::cout << "A_subtrahend: \n" << A_subtrahend << std::endl;
  arma::mat A = H - A_subtrahend;
  A = A.rows(0, A.n_cols-2);
  arma::vec m_zeros(m, arma::fill::zeros);
  // std::cout << "A: \n" << A << std::endl;

  for (arma::uword i=0; i < k; i++) {
    arma::rowvec Qi_lx = Qpred.row(i) - m_alpha(i) * se.row(i);
    // std::cout << "Qi_lx: \n" << Qi_lx << std::endl;
    arma::rowvec Qi_ux = Qpred.row(i) + m_alpha(i) * se.row(i);
    // std::cout << "Qi_ux: \n" << Qi_ux << std::endl;

    arma::uvec aux = find(arma::diff(Qi_lx) < 0);
    // std::cout << "diff(Qi_lx):\n " << arma::diff(Qi_lx) << std::endl;
    // std::cout << "aux.n_elem: " << aux.n_elem << std::endl;
    if (aux.n_elem > 0) {
      arma::vec b = arma::diff(Qpred.row(i)).t();
      // std::cout << "b: " << b << std::endl;
      arma::vec lb = -1 * m_alpha(i) * se.row(i).t();
      // std::cout << "lb: " << lb << std::endl;
      // std::cout << "H: " << H << std::endl;
      // std::cout << "lb: " << lb << std::endl;
      // std::cout << "A: " << A << std::endl;
      // std::cout << "b: " << b << std::endl;
      arma::vec delta_Q = quadprog(H, -lb, -A, b, lb, m_zeros);
      // arma::vec delta_Q = quadprog(DENSE_AUL, H, -lb, -A, b, lb, m_zeros);
      // std::cout << "delta_Q: \n" << delta_Q << std::endl;
      arma::rowvec Qi_lx = Qpred.row(i) + delta_Q.t();
      // std::cout << "Qi_lx: \n" << Qi_lx << std::endl;
    }

    aux = find(arma::diff(Qi_ux) < 0);
    //std::cout << "diff(Qi_ux):\n " << arma::diff(Qi_ux) << std::endl;
    // std::cout << "aux.n_elem: " << aux.n_elem << std::endl;
    if (aux.n_elem > 0) {
      arma::vec b = arma::diff(Qpred.row(i)).t();
      // std::cout << "b: " << b << std::endl;
      arma::vec ub = m_alpha(i) * se.row(i).t();
      // std::cout << "ub: " << ub << std::endl;
      arma::vec delta_Q = quadprog(H, -ub, -A, b, m_zeros, ub);
      // arma::vec delta_Q = quadprog(DENSE_AUL, H, -ub, -A, b, m_zeros, ub);
      // std::cout << "delta_Q: \n" << delta_Q << std::endl;
      arma::rowvec Qi_ux = Qpred.row(i) + delta_Q.t();
      // std::cout << "Qi_ux: \n" << Qi_ux << std::endl;
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
