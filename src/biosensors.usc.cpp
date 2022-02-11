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

#include "NadarayaRegression.h"
#include "RidgeRegression.h"


// [[Rcpp::export]]
Rcpp::List cpp_nadayara_regression(const arma::mat X, const arma::mat t, const arma::mat Y, const arma::mat hs,
                        const arma::umat indices_1, const arma::umat indices_2) {

  bio::nadaraya_struct result = bio::nadayara_regression(X, t, Y, hs, indices_1, indices_2);
  return Rcpp::List::create(
    Rcpp::Named("prediction")   = result.prediction,
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




