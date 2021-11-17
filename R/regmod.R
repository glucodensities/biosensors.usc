## regmod.R: biosensors.usc glue
##
## Copyright (C) 2019 - 2021  Juan C. Vidal and Marcos Matabuena
##
## This file is part of biosensors.usc.
##
## biosensors.usc is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## biosensors.usc is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with biosensors.usc.  If not, see <http://www.gnu.org/licenses/>.

#' @importFrom osqp solve_osqp osqpSettings

#' @title regmod regression
#' @description Performs the Wasserstein regression.
#' @param data A gd object.
#' @param response The name of the vector of observed values (one of the columns of data$variables).
#' @return An object of class wasserstein containing the components:
#' \code{call} The function call.
#' \code{error} The residuals.
#' \code{prediction} The fitted regression.
#' \code{ub} The upper band of the regression.
#' \code{lb} The lower band of the regression.
#' where \code{error}, \code{prediction}, \code{ub}, and \code{lb} are \code{fdata} objects.
#' @export
regmod_regression <- function(data, predictor) {
  if (!is(data, "biosensor"))
    stop("The data must be an object of biosensor class. @seealso biosensors.usc::load_data")

  if (is.null(data$quantiles))
    stop("The data attribute quantiles can not be NULL")

  if (is.null(data$variables))
    stop("The data attribute variables can not be NULL")

  if (!is(data$quantiles, "fdata"))
    stop("The data attribute quantiles must be of type fdata")

  if (!is(data$variables, "data.frame"))
    stop("The data attribute variables must be of type matrix or array")

  if (!(predictor %in% colnames(data$variables)))
    stop("Error: predictor name is not a colname in data$variables.")

  formulax = NULL
  nas <- tryCatch(
    {
      !is.na(data$variables[, predictor])
    },
    error = function(e) {
      message("The An error occured while computing the wassertein regression:\n", e)
    }
  )
  pred <- as.data.frame(data$variables[nas, predictor])
  cuantil <- as.data.frame(data$quantiles$data[nas, ])

  cuadratico = function(prediciones, cotainferior=-10e-5, cotasuperior=800){
    prediciones = as.matrix(prediciones)
    n = dim(prediciones)[1]
    p = dim(prediciones)[2]

    salida = matrix(0, nrow = n, ncol = p)

    for(i in 1:n){

      P = diag(p)
      A = diag(p)*-1
      for(j in 1:(p-1)){
        A[j,j+1]=1
      }
      A[p,p]=0
      u = rep(cotainferior,p-1)
      l = rep(cotasuperior,p-1)
      u = c(u,cotainferior)
      l = c(l,cotasuperior)
      q = -prediciones[i,]
      # u, l change
      settings <- osqpSettings(verbose = FALSE)
      res <- solve_osqp(P, q, A, u, l, settings)
      res = res$x
      salida[i,] = res
    }
    return(salida)
  }

  n = dim(cuantil)[1]
  pcuantil = dim(cuantil)[2]
  px = length(pred)
  nombrecovariables = paste("X", 1:px, sep="")
  nombrecovariablescuantil = paste("Y", 1:pcuantil, sep="")
  Y = cuantil
  colnames(Y)=nombrecovariablescuantil

  X = pred
  colnames(X) = nombrecovariables


  prediciones = matrix(0, n, ncol = pcuantil)
  beta = matrix(0, nrow = px+1, ncol = pcuantil)

  residuos = prediciones
  residuos2 = prediciones
  sd = prediciones
  varbeta = beta

  data = cbind(X,Y)
  formulaaux = paste("X",1:px,sep="")

  if(is.null(formulax)){

    formulax = formulaaux[1]
    if(px==1){


    } else{

      for(i in 2:px){

        formulax = paste(formulax,formulaaux[i],sep="+")
      }

    }
  }

  # print(formulax)

  residuocuadrado= prediciones

  for(i in 1:pcuantil){


    formulaaux2= paste("Y",i,sep="")
    formulaaux2= paste(formulaaux2,"~",sep="")
    formulamedia= paste(formulaaux2, formulax, sep="")
    formulamedia= as.formula(formulamedia)
    # print(formulamedia)
    mmedia= lm(formulamedia, data= data)
    # print(mmedia)
    beta[,i]= as.numeric(mmedia$coefficients)
    prediciones[,i]= as.numeric(mmedia$fitted.values)
    residuos[,i]= mmedia$residuals
    residuocuadrado[,i]= (mmedia$residuals)*(mmedia$residuals)
    data$Yaux= residuocuadrado[,i]

    formulaaux2= paste("Yaux","~",sep="")
    formulavar= paste(formulaaux2, formulax, sep="")
    formulavar= as.formula(formulavar)

    mmvar=lm(formulavar, data= data)
    sd[,i]= sqrt(pmax(0, as.numeric(mmvar$fitted.values), na.rm = TRUE))
    varbeta[,i]= as.numeric(mmvar$coefficients)

  }



  predicciones2= cuadratico(prediciones)

  databeta= data.frame(X, predicciones2[,1])

  colnames(databeta)[px+1]="Y"



  for(i in 1:pcuantil){


    formulaaux2= "Y"
    formulaaux2= paste(formulaaux2,"~",sep="")
    formulamedia= paste(formulaaux2, formulax, sep="")
    formulamedia= as.formula(formulamedia)
    # print(formulamedia)
    databeta[,px+1]= predicciones2[,i]
    mmedia= lm(formulamedia, data= databeta)
    # print(mmedia)

    # print(mmedia)
    beta[,i]= as.numeric(mmedia$coefficients)
    residuocuadrado[,i]= (Y[,i]-predicciones2[,i])^2
    data$Ynew= residuocuadrado[,i]

    formulaaux2= paste("Ynew","~",sep="")
    formulavar= paste(formulaaux2, formulax, sep="")
    formulavar= as.formula(formulavar)
    mmvar=lm(formulavar, data= data)
    # print(mmvar)
    sdaux = sqrt(pmax(0, as.numeric(mmvar$fitted.values), na.rm = TRUE))
    sdaux[is.nan(sdaux)]=0
    sd[,i]= sdaux
    varbeta[,i]= as.numeric(mmvar$coefficients)

  }

  gd.regmod = list(
    "beta"=beta,
    "varbeta"=varbeta,
    "predcrudo"=prediciones,
    "predmedia"=predicciones2,
    "residuosmedia"=residuos,
    "predsd" = sd)

  representar(gd.regmod$predmedia, cuantil, prediciones)

  class(gd.regmod) <- "bregmod"
  return(gd.regmod)
}



representar <- function(aux, aux2, aux3) {
  par(mfrow= c(1,3))
  plot(fdata(aux), main="Contitional mean quantile curve")
  plot(fdata(aux2), main="Observation quantile curves")
  plot(fdata(aux2-aux), main="Residual curves")
  #plot(fdata(aux3-aux), main="Diff means")
}


#' @title regmod prediction
#' @description Performs the Wasserstein regression.
#' @param data A gd object.
#' @param response The name of the vector of observed values (one of the columns of data$variables).
#' @return An object of class wasserstein containing the components:
#' \code{call} The function call.
#' \code{error} The residuals.
#' \code{prediction} The fitted regression.
#' \code{ub} The upper band of the regression.
#' \code{lb} The lower band of the regression.
#' where \code{error}, \code{prediction}, \code{ub}, and \code{lb} are \code{fdata} objects.
#' @export
regmod_prediction <- function(data, Xpred) {

  if (class(data) != "bregmod")
    stop("The data must be an object of bregmod class. ")

  nx= dim(Xpred)[1]

  unos = rep(1,nx)
  matrizdiseño = cbind(unos,Xpred)
  matrizdiseño = as.matrix(matrizdiseño)
  predcrudo = matrizdiseño%*%data$beta

  cuadratico = function(prediciones, cotainferior=-10e-5, cotasuperior=800){
    prediciones = as.matrix(prediciones)
    n = dim(prediciones)[1]
    p = dim(prediciones)[2]

    data = matrix(0, nrow = n, ncol = p)

    for(i in 1:n){

      P = diag(p)
      A = diag(p)*-1
      for(j in 1:(p-1)){
        A[j,j+1]=1
      }
      A[p,p]=0
      u = rep(cotainferior,p-1)
      l = rep(cotasuperior,p-1)
      u = c(u,cotainferior)
      l = c(l,cotasuperior)
      q = -prediciones[i,]
      # u, l change
      settings <- osqpSettings(verbose = FALSE)
      res <- solve_osqp(P, q, A, u, l, settings)
      res = res$x
      data[i,] = res
    }
    return(data)
  }


  predfinal= cuadratico(predcrudo)


  return(predfinal)

}






