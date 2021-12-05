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
#' @importFrom stats median


#' @title ridge_regression
#' @description Performs a Ridge regression.
#' @param data A biosensor object.
#' @param response The name of the scalar response. The response must be a column name in data$variables.
#' @param w A weight function.
#' @param method The distance measure to be used (@seealso parallelDist::parDist). By default manhattan distance.
#' @param type The kernel type ("gaussian" or "lapla"). By default gaussian distance.
#' @return An object containing the components:
#' \code{best_alphas} Best coefficients obtained with leave-one-out cross-validation criteria.
#' \code{best_kernel} The kernel matrix of the best solution.
#' \code{best_sigma} The sigma parameter of the best solution.
#' \code{best_lambda} The lambda parameter of the best solution.
#' \code{sigmas} The sigma parameters used in the fitting according to the median heuristic fitting criteria.
#' \code{predictions} A matrix of predictions.
#' \code{r2} R-square of the different models fitted.
#' \code{error} Mean squared-error of the different models fitted.
#' \code{predictions_cross} A matrix of predictions obtained with leave-one-out cross-validation criteria.
#' @usage
#' ridge_regression(data, response, w=NULL, method="manhattan", type="gaussian")
#' @examples
#' # Data extracted from the paper: Hall, H., Perelman, D., Breschi, A., Limcaoco, P., Kellogg, R.,
#' # McLaughlin, T., Snyder, M., Glucotypes reveal new patterns of glucose dysregulation, PLoS
#' # biology 16(7), 2018.
#' file1 = system.file("extdata", "data_1.csv", package = "biosensors.usc")
#' file2 = system.file("extdata", "variables_1.csv", package = "biosensors.usc")
#' data = load_data(file1, file2)
#' regm = ridge_regression(data, "BMI")
#' @export
ridge_regression = function(data, response, w=NULL, method="manhattan", type="gaussian") {
# ridge_regression = function(X, Y, w=1, method="manhattan", type="gaussian") {

  nas <- tryCatch(
    {
      !is.na(data$variables[, response])
    },
    error = function(e) {
      message("The An error occured while computing the wassertein regression:\n", e)
    }
  )
  pred <- as.data.frame(data$variables[nas, response])
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

