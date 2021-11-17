## ridge.R: biosensors.usc glue
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

#' @importFrom parallelDist parDist


#' @title ridge_regression
#' @description Generates a quantile reression model V + V2 * v + tau * V2 * Q0 where Q0 is a truncated random variable, v = 2 * X, tau = 2 * X, V ~ Unif(-1, 1), V2 ~ Unif(-1, -1), V3 ~ Unif(0.8, 1.2), and E(V|X) = tau * Q0;
#' @param X ToDo.
#' @param Y ToDo.
#' @param w ToDo.
#' @param method The distance measure to be used (@seealso parallelDist::parDist)
#' @param type The kernel type ("gaussian" or "lapla")
#' @return An object of class wasserstein containing the components:
#' \code{call} The function call.
#' \code{error} The residuals.
#' \code{prediction} The fitted regression.
#' \code{ub} The upper band of the regression.
#' \code{lb} The lower band of the regression.
#' @export
ridge_regression = function(data, predictor, w=NULL, method="manhattan", type="gaussian") {
# ridge_regression = function(X, Y, w=1, method="manhattan", type="gaussian") {


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

  if (is.null(w))
    w = rep(1, nrow(cuantil))

  X = as.matrix(cuantil)
  Y = as.matrix(pred)
  #w = ifelse(is.na(w),0,w)
  #Y = ifelse(is.na(Y),0,Y)
  X = X[w>0,]
  Y = Y[w>0,]
  w = w[w>0]
  n = dim(X)[1]
  W = matrix(0, nrow=n, ncol=n)
  diag(W) = w
  W2 = W/sum(W)
  distancia = parallelDist::parDist(X, method = method)
  distancia = as.matrix(distancia)
  if(type == "gaussian") {
    distancia = as.matrix(distancia)
    mediana = median(distancia[distancia>0]^2)
    mediana = sqrt(mediana)
    potencias = seq(0.3,3.5, length=35)
    sigmas = mediana^(potencias)
    lambdas = seq(0.3,2,length=20)
    expandido = expand.grid(lambdas, potencias)
    distancia = (distancia)^2
    mediana = mediana
    m2 = cpp_ridge_regression(distancia, as.vector(Y), as.matrix(W), as.vector(w), as.vector(lambdas), as.vector(sigmas))
    return(m2)
  }
  if(type == "lapla") {
    distancia = as.matrix(distancia)
    mediana = median(distancia[distancia>0]^2)
    mediana = sqrt(mediana)
    potencias = seq(0.3,3.5, length=35)
    sigmas = mediana^(potencias)
    lambdas = seq(0.3,2,length=20)
    expandido = expand.grid(lambdas, potencias)
    m2 = cpp_ridge_regression(distancia, as.vector(Y), as.matrix(W), as.vector(w), as.vector(lambdas), as.vector(sigmas))
    return(m2)
  }
}

# n=1000
# p= 2
# X= matrix(runif(n*p,0,1), nrow= n, ncol=p)
# Y= X[,1]+rnorm(n)
# m=ridge_analysis(X, Y ,w= rep(1,n))


#' @title ridge_prediction
#' @description Generates a quantile reression model V + V2 * v + tau * V2 * Q0 where Q0 is a truncated random variable, v = 2 * X, tau = 2 * X, V ~ Unif(-1, 1), V2 ~ Unif(-1, -1), V3 ~ Unif(0.8, 1.2), and E(V|X) = tau * Q0;
#' @param m ToDo.
#' @param X ToDo.
#' @param Xpred ToDo.
#' @param method The distance measure to be used (@seealso parallelDist::parDist)
#' @param type The kernel type ("gaussian" or "lapla")
#' @return An object of class wasserstein containing the components:
#' \code{call} The function call.
#' \code{error} The residuals.
#' \code{prediction} The fitted regression.
#' \code{ub} The upper band of the regression.
#' \code{lb} The lower band of the regression.
#' @export
ridge_prediction = function(m, X, Xpred, method="manhattan", type="gaussian"){

  alpha = m$best_alphas
  Xnew = rbind(X,Xpred)
  n1 = dim(X)[1]
  n2 = n1+dim(Xpred)[1]

  sel1 = (n1+1):n2
  sel2 = 1:n1


  distancia = parallelDist::parDist(Xnew, method = method)
  distancia = as.matrix(distancia)
  distancia = distancia[sel1,sel2]

  if(type == "gaussian"){
    distancia = as.matrix(distancia)
    potencias = seq(0.3,3.5, length=35)
    sigmas = m$best_sigma
    distancia = (distancia)^2
    mediana = sigmas
    kernel = exp(-distancia*(1/mediana))
    pred = kernel%*%as.matrix(alpha)
  }

  if(type == "lapla"){
    distancia = as.matrix(distancia)
    potencias = seq(0.3,3.5, length=35)
    sigmas = m$best_sigma
    distancia = distancia
    mediana = sigmas
    kernel = exp(-distancia*(1/mediana))
    pred = kernel%*%as.matrix(alpha)
  }

  return(pred)


}


# ntest = 10
# Xpred = matrix(runif(ntest*p,0,1), nrow= ntest, ncol=p)
# print(ridge_prediction(m=m,X=X,Xpred))
# print(Xpred%*%as.matrix(c(1,0)))
