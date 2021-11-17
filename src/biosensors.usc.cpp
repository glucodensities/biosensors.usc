// biosensors.usc.cpp: biosensors.usc glue
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

#include <iostream>
#include <iterator>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "AlglibSolvers.h"
#include "WassersteinRegression.h"
#include "NadarayaRegression.h"
#include "RidgeRegression.h"
#include "ConfidenceBand.h"


//' This function perform Frechet regression with the Wasserstein distance.
//'
//' @param xfit A nxp matrix of predictor values for fitting (do not include a column for the intercept).
//' @param q A nxm matrix of quantile density functions. q(i, :) is a 1xm vector of quantile density function values on an equispaced grid on [0, 1].
//' @param Q0 A 1xn array of quantile function values at 0.
//' @param xpred A kxp matrix of input values for regressors for prediction.
//' @param t A 1xm vector - common grid for all quantile density functions in q.  If missing, defaults to linspace(0, 1, m);  For best results, should use a finer grid than for quantle estimation, especially near the boundaries.
//' @param qdmin A positive lower bound on the estimated quantile densites. Defaults to 1e-6.
//'
//' @return An object containing the components:
//' \code{xpred} See input of same name.
//' \code{qpred} A kxN array.  qpred(l,:) is the regression prediction of q (the quantile density) given X = xpred(l, :)'
//' \code{Qpred} A kxm array. Qpred(l, :) is the regression prediction of Q given X = xpred(l, :)'
//' \code{fpred} A kxm array. fpred(l, :) is the regression prediction of f (the density) given X = xpred(l, :)', evaluated on the grid Qpred(l, :)
//' \code{xfit} See input of same name.
//' \code{qfit} A nxN array. qfit(l, :) is the regression prediction of q given X = xfit(l, :)'
//' \code{Qfit} A nxm array. Qfit(l, :) is the regression prediction of Q given X = xfit(l, :)'
//' \code{fpred} A kxm array. fpred(l, :) is the regression prediction of f (the density) given X = xpred(l, :)', evaluated on the grid Qfit(l, :)
//' \code{QP_used} A flag indicating whether OLS fits all satisfied the constraints (=0) or if the quadratic program was used in fitting (=1).
// [[Rcpp::export]]
Rcpp::List cpp_wasserstein_regression(const arma::mat xfit, const arma::mat q, const arma::mat Q0,
                 const arma::mat xpred, const arma::vec t, const double qdmin) {
  bio::regression_struct result = bio::wasserstein_regression(xfit, q, Q0, xpred, t, qdmin);
  return Rcpp::List::create(
    Rcpp::Named("q")       = q,
    Rcpp::Named("Q0")      = Q0,
    Rcpp::Named("t")       = t,
    Rcpp::Named("qdmin")   = qdmin,
    Rcpp::Named("xfit")    = result.xfit,
    Rcpp::Named("xpred")   = result.xpred,
    Rcpp::Named("Qfit")    = result.Qfit,
    Rcpp::Named("Qpred")   = result.Qpred,
    Rcpp::Named("qfit")    = result.qfit,
    Rcpp::Named("qpred")   = result.qpred,
    Rcpp::Named("ffit")    = result.ffit,
    Rcpp::Named("fpred")   = result.fpred,
    Rcpp::Named("QP_used") = result.QP_used
  );
}



//' This function computes intrinsic confidence bands for Wasserstein regression.
//'
//' @param xfit A nxp matrix of predictor values for fitting (do not include a column for the intercept).
//' @param xpred A kxp matrix of input values for regressors for prediction.
//' @param Q_obs A nxm matrix of quantile functions. Q_obs(i, :) is a 1xm vector of quantile function values on grid t_vec.
//' @param q_obs A nxm matrix of quantile density functions. q_obs(i, :) is a 1xm vector of quantile density function values on grid t_vec.
//' @param t_vec A 1xm vector - common grid for all quantile density functions in Q_obs, q_obs, q_prime_obs.
//' @param alpha The significant level is 100*(1 - alpha).
//'
//' @return An object containing the components:
//' \code{Q_lx} Lower bound of confidence bands in terms of density functions.
//' \code{Q_ux} Upper bound of confidence bands in terms of density functions.
//' \code{Qpred} Fitted density function at xpred.
// [[Rcpp::export]]
Rcpp::List cpp_confidence_band(const arma::mat xfit, const arma::mat xpred, const arma::mat Q_obs,
                     const arma::mat q_obs, const arma::vec t_vec, const double alpha) {
  bio::confidence_struct result = bio::confidence_band(xfit, xpred, Q_obs, q_obs, t_vec, alpha);
  return Rcpp::List::create(
    Rcpp::Named("xfit")   = xfit,
    Rcpp::Named("xpred")  = xpred,
    Rcpp::Named("Q_obs")  = Q_obs,
    Rcpp::Named("t_vec")  = t_vec,
    Rcpp::Named("alpha")  = alpha,
    Rcpp::Named("Qpred")  = result.Qpred,
    Rcpp::Named("Q_lx")   = result.Q_lx,
    Rcpp::Named("Q_ux")   = result.Q_ux,
    Rcpp::Named("fpred")  = result.fpred
  );
}




// [[Rcpp::export]]
Rcpp::List cpp_nadayara_regression(const arma::mat X, const arma::mat t, const arma::mat Y, const arma::mat hs,
                        const arma::umat indices_1, const arma::umat indices_2) {

  bio::nadaraya_struct result = bio::nadayara_regression(X, t, Y, hs, indices_1, indices_2);
  return Rcpp::List::create(
    Rcpp::Named("coefficients") = result.coefficients,
    Rcpp::Named("residuals")    = result.residuals,
    Rcpp::Named("r2")           = result.r2,
    Rcpp::Named("error")        = result.error,
    Rcpp::Named("r2_global")    = result.r2_global
  );
}

// [[Rcpp::export]]
arma::mat cpp_nadayara_prediction(const arma::mat X, const arma::mat t, const arma::mat Y, const arma::mat hs,
                                   const arma::umat indices_1, const arma::umat indices_2) {

  return bio::nadayara_predicion(X, t, Y, hs, indices_1, indices_2);
}


// [[Rcpp::export]]
Rcpp::List cpp_ridge_regression(const arma::mat dist, const arma::vec Y, const arma::mat W,
                            const arma::vec w, const arma::vec lambdas, const arma::vec sigmas) {

  bio::ridge_struct result = bio::ridge_regression(dist, Y, W, w, lambdas, sigmas);

  return Rcpp::List::create(
    Rcpp::Named("best_alphas") = result.best_alphas,
    Rcpp::Named("sigmas") = result.sigmas,
    Rcpp::Named("predictions") = result.predictions,
    Rcpp::Named("r2")= result.r2,
    Rcpp::Named("error") = result.error,
    Rcpp::Named("predictions_cross") = result.predictions_cross,
    Rcpp::Named("best_kernel") = result.best_kernel,
    Rcpp::Named("best_sigma") = result.best_sigma,
    Rcpp::Named("best_lambda") = result.best_lambda
  );
}




