// NadarayaRegression.h: biosensors.usc glue
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

#ifndef _NADARAYA_REGRESSION_H // include guard
#define _NADARAYA_REGRESSION_H

#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>


namespace bio {

struct nadaraya_struct {
  arma::mat prediction;
  arma::mat residuals;
  arma::mat r2;
  arma::mat error;
  arma::mat r2_global;

};


inline arma::mat trapecio(arma::mat X, arma::mat Y) {
  if (X.n_rows != Y.n_rows || X.n_cols != Y.n_cols)
    throw std::invalid_argument("Arguments 'x' and 'y' must be matrices of the same dimension");
  // bool transposed = false;
  // arma::uword m = 1;
  // if (X.n_cols == 1) {
  //   transposed = true;
  // } else {
  // bool transposed = false;
  arma::uword m = 1;
  if (X.n_cols != 1) {
  //   transposed = true;
  // } else {
    X = X.t();
    Y = Y.t();
  }

  arma::uword n = X.n_rows;
  arma::mat C = arma::zeros(n,m);
  arma::mat tmp = (X.rows(1,n-1) - X.rows(0,n-2)) % (Y.rows(1,n-1) + Y.rows(0,n-2));
  arma::mat result;
  if (tmp.n_rows == 1) {
    C.row(1) = 0.5 * tmp;
    result = C.row(n-1);
  } else {
    C.rows(1,n-1) = 0.5 * cumsum(tmp);
    result = C.row(n-1);
  }
  return result;
}


inline arma::mat eucdistance1(arma::mat X, arma::mat t){
  arma::uword n = X.n_rows;
  // arma::uword m = X.n_cols;
  arma::mat distancias(n,n);
  // double distanciaaux=0;
  for(arma::uword i=0; i<n; i++) {
    for(arma::uword j=0; j<n; j++) {
      distancias(i,j) = sqrt(trapecio(t,((X.row(i)-X.row(j))%(X.row(i)-X.row(j))).t())(0));
    }
  }
  return distancias;
}


inline arma::mat eucdistance2(arma::mat X, arma::mat t, arma::mat x){
  arma::uword n = X.n_rows;
  arma::uword nx = x.n_rows;

  arma::mat distancias(nx,n);
  for(arma::uword i=0; i<nx; i++) {
    for(arma::uword j=0; j<n; j++) {
      distancias(i,j) = sqrt(trapecio(t,((x.row(i)-X.row(j))%(x.row(i)-X.row(j))).t())(0));
    }
  }
  return distancias;
}



inline arma::vec gaussiankernelinicial(arma::vec x,double h){
  arma::vec salida;
  salida = double(2/sqrt(double(2) * M_PI))*exp(-(x/h)%(x/h)*0.5);
  return salida;
}


inline arma::vec triangular(arma::vec x, double h) {
  int l = x.size();
  arma::vec salida(l);

  for(int i=0; i<l; i++){
    if(double(x(i)/h) < 0 || double(x(i)/h) > 1) {
      salida(i)=0;
    } else {
      salida(i)=double(double(35)/double(16))*pow(double(1-pow(x(i)/h,2)),3);
    }
  }
  return salida;
}

/*
nadaraya_struct nadayara_prediction(const arma::mat distancias, const arma::mat X, const arma::mat t,
                                    const arma::mat Y, const arma::mat hs) {
  arma::uword n = X.n_rows;
  arma::vec media= mean(Y);
  arma::vec mediavectorial(n);
  mediavectorial.fill(media(0));
  arma::vec residuos=Y-mediavectorial;
  arma::uword nh = hs.n_rows;
  arma::mat prediciones(n,nh);
  arma::mat residuosglobal(n,nh);
  arma::mat R2(nh,1);
  arma::vec aux(n);
  arma::vec aux2(n);
  arma::vec aux3(n);
  arma::vec residuosaux(n);

  for(arma::uword j=0; j<nh; j++) {
    for(arma::uword i=0; i<n; i++) {
      aux2 = distancias.row(i).t();
      aux = gaussiankernelinicial(aux2, hs(j));
      aux3(i) = sum((aux%Y)) / sum(aux);
    }
    prediciones.col(j) = aux3;
    residuosaux= Y - aux3;
    residuosaux= 1 - double(sum(residuosaux%residuosaux) / sum(residuos%residuos));
    residuosglobal.col(j) = Y - aux3;
    R2.row(j) = residuosaux(0);
  }

  nadaraya_struct result;
  result.coefficients = prediciones;
  result.residuals = residuosglobal;
  result.r2 = R2;
  return result;
}*/




nadaraya_struct nadayara_regression(const arma::mat X, const arma::mat t, const arma::mat Y, const arma::mat hs,
                                    const arma::umat indices1, const arma::umat indices2){
  arma::uword n = X.n_rows;
  // arma::uword m = X.n_cols;
  arma::mat distancias(n,n);
  arma::vec media= mean(Y);
  arma::vec mediavectorial(n);
  mediavectorial.fill(media(0));
  arma::vec residuos = Y - mediavectorial;
  arma::uword nh = hs.n_rows;

  // salida rendimiento modelos

  arma::mat prediciones(n,nh);
  arma::mat residuosglobal(n,nh);
  arma::mat R2(nh,1);
  arma::mat error(nh,1);

  // float h = 1;
  arma::uword n1x = indices1.n_rows;
  arma::uword n1p = indices1.n_cols;
  arma::uword n2x = indices2.n_rows;
  distancias= eucdistance1(X, t);

  arma::mat Y1(n1x,1);
  arma::mat Y2(n2x,1);
  arma::mat distancias2(n2x,n1x);
  arma::umat indices1aux(n1x,1);
  arma::umat indices2aux(n2x,1);
  arma::mat R2validacion(nh,n1p);
  arma::mat R2validacionaux(nh,1);

  for(arma::uword l=0; l<n1p;l++) {

    indices1aux= indices1.col(l);
    indices2aux= indices2.col(l);

    Y1= Y.rows(indices1aux);
    Y2= Y.rows(indices2aux);

    arma::vec auxx(n1x);
    arma::vec aux2x(n1x);
    arma::vec aux3x(n2x);
    arma::mat prediciones2(n2x,nh);

    distancias2 = distancias.rows(indices2aux);
    distancias2 = distancias2.cols(indices1aux);

    arma::vec residuosauxx(n1x);

    arma::vec mediax = mean(Y2);
    arma::vec mediavectorialx(n2x);
    mediavectorialx.fill(mediax(0));

    arma::vec residuosx = Y2-mediavectorialx;

    for(arma::uword j=0; j<nh; j++){
      for(arma::uword i=0; i<n2x; i++) {
        aux2x = distancias2.row(i).t();
        auxx = gaussiankernelinicial(aux2x,hs(j));
        aux3x(i) =  sum((auxx%Y1)) / sum(auxx);

      }
      prediciones2.col(j) = aux3x;
      residuosauxx = Y2-aux3x;
      residuosauxx = sum(residuosauxx%residuosauxx);
      //residuosauxx= 1-double(sum(residuosauxx%residuosauxx)/sum(residuosx%residuosx));
      //residuosglobal.col(j)= Y2-aux3x;
      R2validacionaux.row(j) = residuosauxx(0);
    }
    R2validacion.col(l) = R2validacionaux;
  }

  arma::vec aux(n);
  arma::vec aux2(n);
  arma::vec aux3(n);
  arma::vec residuosaux(n);

  for(arma::uword j=0; j<nh; j++) {
    for(arma::uword i=0; i<n; i++) {
      aux2 = distancias.row(i).t();
      aux = gaussiankernelinicial(aux2,hs(j));
      aux3(i) =  sum((aux%Y)) / sum(aux);
    }
    prediciones.col(j) = aux3;

    residuosaux = Y - aux3;
    error.row(j) = sum(residuosaux%residuosaux);

    residuosaux = 1-double(sum(residuosaux%residuosaux) / sum(residuos % residuos));
    residuosglobal.col(j) = Y - aux3;
    R2.row(j) =residuosaux(0);
  }

  nadaraya_struct result;
  result.prediction = prediciones;
  result.residuals = residuosglobal;
  result.r2 = R2;
  result.error = error;
  result.r2_global = R2validacion;
  return result;
}





arma::mat nadayara_predicion(const arma::mat X, const arma::mat t, const arma::mat Y, const arma::mat hs,
                             const arma::umat indices1, const arma::umat indices2){
  arma::uword n = X.n_rows;
  // arma::uword m = X.n_cols;
  arma::mat distancias(n,n);
  arma::vec media= mean(Y);
  arma::vec mediavectorial(n);
  arma::uword nh = hs.n_rows;

  // float h = 1;
  arma::uword n1x = indices1.n_rows;
  arma::uword n1p = indices1.n_cols;
  arma::uword n2x = indices2.n_rows;
  distancias= eucdistance1(X, t);

  // salida de la funcion
  arma::mat prediccionesfinales(n2x,nh);
  arma::mat Y1(n1x,1);
  arma::mat distancias2(n2x,n1x);
  arma::umat indices1aux(n1x,1);
  arma::umat indices2aux(n2x,1);

  for(arma::uword l=0; l<n1p; l++) {
    indices1aux= indices1.col(l);
    indices2aux= indices2.col(l);

    Y1= Y.rows(indices1aux);

    arma::vec auxx(n1x);
    arma::vec aux2x(n1x);
    arma::vec aux3x(n2x);
    arma::mat prediciones2(n2x,nh);

    distancias2 = distancias.rows(indices2aux);
    distancias2 = distancias2.cols(indices1aux);

    for(arma::uword j=0; j<nh; j++){
      for(arma::uword i=0; i<n2x; i++) {
        aux2x = distancias2.row(i).t();
        auxx = gaussiankernelinicial(aux2x,hs(j));
        aux3x(i) =  sum((auxx%Y1)) / sum(auxx);

      }
      prediciones2.col(j) = aux3x;

    }
    prediccionesfinales = prediciones2;

  }

  return prediccionesfinales;
}


}

#endif
