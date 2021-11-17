// RidgeRegression.h: biosensors.usc glue
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

#ifndef _RIDGE_REGRESSION_H // include guard
#define _RIDGE_REGRESSION_H

#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>


namespace bio {

struct ridge_struct {
  arma::mat predictions;
  arma::mat predictions_cross;
  arma::vec error;
  arma::vec r2;
  arma::mat best_kernel;
  arma::vec best_alphas;
  arma::vec sigmas;
  double best_sigma;
  arma::vec lambdas;
  double best_lambda;
};

ridge_struct ridge_regression(arma::mat distance, arma::vec Y, arma::mat W, arma::vec w,
                                 arma::vec lambdas, arma::vec sigmas){
  arma::uword n = distance.n_rows;
  arma::uword np = lambdas.n_elem; //lambda.size();
  arma::uword np2 = sigmas.n_elem; //sigmas.size();
  arma::uword combinations = np * np2;
  arma::mat kernel(n,n);
  arma::mat best_kernel(n,n);
  arma::vec best_alphas(n);
  double best_sigma = 0;
  double best_lambda = 0;
  double taum2 = 0;
  arma::mat H;
  arma::mat results(combinations, n);
  arma::mat predictions(combinations, n);
  arma::mat residuos(combinations, n);
  arma::vec error(combinations);
  arma::mat predictions_cross(combinations, n);
  arma::vec r2(combinations);
  arma::vec alpha;
  double aux;
  arma::vec residuosiniciales(n);
  arma::vec media(n);
  media.fill(arma::sum((Y%w)/arma::sum(w)));
  residuosiniciales= Y-media;
  double cuenta  = 0;
  double cuenta1 = 0;
  double cuenta2 = 0;
  double best_error = 100000000;
  arma::vec  auxpred;

  int contar = -1;
  for(arma::uword j=0; j < np2; j++){
    taum2 = sigmas(j);
    kernel = arma::exp(-distance*(1/double(taum2)));
    for(arma::uword i=0; i < np; i++){
      contar += 1;
      aux = lambdas(i);
      alpha=(inv(W * kernel + aux * arma::eye(n,n)))*((W*Y));
      results.row(i) = alpha.t();
      auxpred= kernel*alpha;
      predictions.row(contar) = auxpred.t();
      residuos.row(contar) = (Y-auxpred).t();
      cuenta1 = arma::sum(w % (Y-auxpred) % (Y-auxpred));
      cuenta2 = arma::sum(w % (residuosiniciales % residuosiniciales));
      cuenta= 1- cuenta1 / cuenta2;
      r2(contar) = cuenta;
      H= kernel * (arma::inv(W * kernel + aux * arma::eye(n,n))) * W;
      arma::mat auxn = arma::inv(arma::eye(n, n) - arma::diagmat(H));
      arma::mat auxn2 = (Y - auxpred) % auxn.diag();
      predictions_cross.row(contar) = auxn2.t();
      double sal = arma::sum(w % auxn2 % auxn2) / arma::sum(w);
      if(sqrt(sal) < best_error) {
        best_kernel = kernel;
        best_sigma = taum2;
        best_lambda = aux;
        best_error = sqrt(sal);
        best_alphas = alpha;
      }
      error(contar) = sqrt(sal);
    }
  }

  ridge_struct result;
  result.predictions = predictions;
  result.predictions_cross = predictions_cross;
  result.error = error;
  result.r2 = r2;
  result.best_kernel = best_kernel;
  result.best_alphas = best_alphas;
  result.sigmas = sigmas;
  result.best_sigma = best_sigma;
  result.lambdas = lambdas;
  result.best_lambda = best_lambda;
  return result;
}






}

#endif
