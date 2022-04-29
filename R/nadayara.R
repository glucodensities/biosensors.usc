## regression: biosensors.usc glue
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

# result.predition = prediciones;
# result.residuals = residuosglobal;
# result.r2 = R2;
# result.error = error;
# result.r2_global = R2validacion;

#' @importFrom graphics par plot axis box mtext legend
#' @importFrom stats complete.cases
#' @importFrom methods is

#' @title nadayara_regression
#' @description Functional non-parametric Nadaraya-Watson regression with 2-Wasserstein distance, using as predictor the distributional representation and as response a scalar outcome.
#' @param data A biosensor object.
#' @param response The name of the scalar response. The response must be a column name in data$variables.
#' @return An object of class bnadaraya:
#' \code{prediction} The Nadaraya-Watson prediction for each point of the training data at each h=seq(0.8, 15, length=200).
#' \code{r2} R2 estimation for the training data at each h=seq(0.8, 15, length=200).
#' \code{error} Standard mean-squared error after applying leave-one-out cross-validation for the training data at each h=seq(0.8, 15, length=200).
#' \code{data} A data frame with biosensor raw data.
#' \code{response} The name of the scalar response.
#' @usage
#' nadayara_regression(data, response)
#' @examples
#' # Data extracted from the paper: Hall, H., Perelman, D., Breschi, A., Limcaoco, P., Kellogg, R.,
#' # McLaughlin, T., Snyder, M., Glucotypes reveal new patterns of glucose dysregulation, PLoS
#' # biology 16(7), 2018.
#' file1 = system.file("extdata", "data_1.csv", package = "biosensors.usc")
#' file2 = system.file("extdata", "variables_1.csv", package = "biosensors.usc")
#' data = load_data(file1, file2)
#' nada = nadayara_regression(data, "BMI")
#' @export
nadayara_regression <- function(data, response) {
  if (!is(data, "biosensor"))
    stop("Error: data must be an object of biosensor class. @seealso biosensors.usc::load_data")

  if (is.null(data$quantiles))
    stop("The data attribute quantiles can not be NULL")

  if (is.null(data$variables))
    stop("The data attribute variables can not be NULL")

  if (!is(data$quantiles, "fdata"))
    stop("The data attribute quantiles must be of type fdata")

  if (!is(data$variables, "data.frame"))
    stop("The data attribute variables must be of type matrix or array")

  if (!(response %in% colnames(data$variables)))
    stop("Error: response name is not a colname in data$variables.")

  nas <- tryCatch(
    {
      !is.na(data$variables[, response])
    },
    error = function(e) {
      message("An error occured while computing the wassertein regression:\n", e)
    }
  )

  hs <- seq(0.8, 15, length = 200)
  X <- data$quantiles$data[nas, ]
  n <- dim(X)[1]
  p <- dim(X)[2]
  t <- data$quantiles$argvals
  Y <- data$variables[nas, response]

  conjunto <- data.frame(X, Y)
  aux <- stats::complete.cases(conjunto)
  conjunto <- conjunto[stats::complete.cases(conjunto), ]
  X <- conjunto[, 1:p]
  Y <- conjunto[, p + 1]
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  t <- as.matrix(t)
  hs <- as.matrix(hs)

  n <- dim(X)[1]
  p <- dim(X)[2]

  training_data <- n - 1
  test_data <- 1
  cross <- n

  indices <- 0:(n - 1)
  indices1 <- matrix(0, nrow = training_data, ncol = cross)
  indices2 <- matrix(0, nrow = test_data, ncol = cross)
  indicesaux <- 1:n

  for (j in 1:cross) {
    generar <- indicesaux[-c(j)]
    indices1[, j] <- indices[generar]
    indices2[, j] <- indices[c(j)]
  }

  indices1 <- as.matrix(indices1)
  indices2 <- as.matrix(indices2)

  res <- cpp_nadayara_regression(X, t, Y, hs, indices1, indices2)

  predictivo <- apply(res$r2_global, 1, function(x) {
    sqrt(mean(x))
  })
  ventanas <- hs
  cuantos <- sum(is.nan(predictivo))
  ventanas2 <- ventanas[!is.nan(predictivo)]
  predictivo2 <- predictivo[!is.nan(predictivo)]
  R2 <- res$r2[!is.nan(predictivo)]

  rango1 <- c(min(predictivo2), max(predictivo2))
  rango2 <- c(min(R2), max(R2))

  time <- ventanas2
  cell.density <- predictivo2
  betagal.abs <- R2
  betagal.abs
  ## add extra space to right margin of plot within frame

  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))                 # code line i + 1

  graphics::par(mfrow = c(1, 1))


  graphics::plot(time, betagal.abs,
                 col = "#0073C2FF",
                 xlab = "Smoothing-parameter", ylab = NA, type = "p",
                 main = "Performance model vs. smoothing-parameter\n(leave-one-out cross-validation)"
  )

  graphics::mtext(side = 2, line = 3, "R-square")

  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))                 # code line i + 1

  new_par <- old_par <- par("mar")
  new_par[4] <- old_par[2]
  graphics::par(mar = new_par)

  graphics::par(new = T)
  graphics::par(mar = old_par)

  graphics::plot(time, cell.density,
                 col = "#FC4E07", xlab = NA, ylab = NA, axes = FALSE, type = "p"
  )

  graphics::axis(side = 4)
  graphics::mtext(side = 4, line = 3, 'Standard deviation')

  graphics::legend("topright",
                   legend = c("R-square", "Error"),
                   pch = c(1, 1), col = c("#0073C2FF", "#FC4E07")
  )

  gd.regression <- list(prediction = res$prediction, r2 = R2, error = cell.density, data = data, response = response)
  class(gd.regression) <- "bnadaraya"
  return(gd.regression)
}




#' @title nadayara_prediction
#' @description Functional non-parametric Nadaraya-Watson prediction with 2-Wasserstein distance.
#' @param nadaraya A Nadaraya regression object.
#' @param Qpred Quantile curves that will be used in the predictions
#' @param hs Smoothing parameters for the predictions, by default hs = seq(0.8, 15, length = 200)
#' @return An object of class bnadarayapred:
#' \code{prediction} The Nadaraya-Watson prediction for the test data at each value of hs.
#' \code{hs} Hs values used for the prediction.
#' @usage
#' nadayara_prediction(nadaraya, Qpred, hs=NULL)
#' @examples
#' # Data extracted from the paper: Hall, H., Perelman, D., Breschi, A., Limcaoco, P., Kellogg, R.,
#' # McLaughlin, T., Snyder, M., Glucotypes reveal new patterns of glucose dysregulation, PLoS
#' # biology 16(7), 2018.
#' file1 = system.file("extdata", "data_1.csv", package = "biosensors.usc")
#' file2 = system.file("extdata", "variables_1.csv", package = "biosensors.usc")
#' data = load_data(file1, file2)
#' nada = nadayara_regression(data, "BMI")
#' # Example of prediction with the column mean of quantiles
#' npre = nadayara_prediction(nada, t(colMeans(data$quantiles$data)))
#' @export
nadayara_prediction <- function(nadaraya, Qpred, hs=NULL){
  if (!is(nadaraya, "bnadaraya"))
    stop("The data must be an object of bnadaraya class. ")

  if (is.null(Qpred))
    stop("The data attribute Qpred can not be NULL")

  if (!(is(Qpred, "matrix") || is(Qpred, "array")))
    stop("The data must be an object of bnadaraya class. ")

  # falta meter excepciones Qpred

  if(is.null(hs)==TRUE){
    hs <- seq(0.8, 15, length = 200)
  }

  X <- as.matrix(nadaraya$data$quantiles$data)
  n <- dim(X)[1]
  p <- dim(X)[2]
  t <- nadaraya$data$quantiles$argvals
  Y <- nadaraya$data$variables[, nadaraya$response]

  Xtest= as.matrix(Qpred)
  nXtest= dim(Xtest)[1]
  indicesXtest=(n+1):(n+nXtest)

  conjunto <- data.frame(X, Y)
  conjunto <- as.matrix(conjunto)
  aux <- stats::complete.cases(conjunto)
  conjunto <- conjunto[stats::complete.cases(conjunto), ]
  X <- conjunto[, 1:p]
  Y <- conjunto[, p + 1]
  X <- rbind(X,Xtest)
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  t <- as.matrix(t)
  hs <- as.matrix(hs)

  training_data <- n
  test_data <- nXtest
  cross <- 1

  indices <- 0:(n - 1)
  indices1 <- matrix(0, nrow = training_data, ncol = cross)
  indices2 <- matrix(0, nrow = test_data, ncol = cross)
  indicesaux <- 1:n

  for (j in 1:cross) {
    indices1[, j] <- indices
    indices2[, j] <- indicesXtest-1
  }

  indices1 <- as.matrix(indices1)
  indices2 <- as.matrix(indices2)

  res <- cpp_nadayara_prediction(X, t, Y, hs, indices1, indices2)

  gd.pred <- list(prediction = res, hs = hs)
  class(gd.pred) <- "bnadarayapred"
  return(gd.pred)
}
