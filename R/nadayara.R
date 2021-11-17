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

#' @importFrom graphics par plot axis box mtext legend
#' @importFrom stats complete.cases

#' @title regression_analysis
#' @description Performs an analysis of the performance of the regression vs a smoothing parameter.
#' @param data A biosensor object.
#' @param predictor The name of the vector of observed values (one of the columns of data$variables).
#' @return An object of class performance containing the components:
#' \code{call} The function call
#' \code{error} The residuals.
#' \code{prediction} The fitted regression.
#' \code{ub} The upper band of the regression.
#' \code{lb} The lower band of the regression.
#' where \code{error}, \code{prediction}, \code{ub}, and \code{lb} are \code{fdata} objects.
#' @export
nadayara_regression <- function(data, predictor) {
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

  if (!(predictor %in% colnames(data$variables)))
    stop("Error: predictor name is not a colname in data$variables.")

  nas <- tryCatch(
    {
      !is.na(data$variables[, predictor])
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
  Y <- data$variables[nas, predictor]

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

  graphics::par(mfrow = c(1, 1))


  graphics::plot(time, betagal.abs,
                 col = "#0073C2FF",
                 xlab = "Smoothing-parameter", ylab = NA, type = "p",
                 main = "Performance model vs. smoothing-parameter\n(leave-one-out cross-validation)"
  )

  graphics::mtext(side = 2, line = 3, "R-square")

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

  gd.regression <- list(nadayara = res, r2 = R2, error = cell.density, data = data, predictor = predictor)
  class(gd.regression) <- "bnadaraya"
  return(gd.regression)
}




#' @title nadayara_prediction
#' @description Performs an analysis of the performance of the regression vs a smoothing parameter.
#' @param data A biosensor object.
#' @param predictor The name of the vector of observed values (one of the columns of data$variables).
#' @return An object of class performance containing the components:
#' \code{call} The function call
#' \code{error} The residuals.
#' \code{prediction} The fitted regression.
#' \code{ub} The upper band of the regression.
#' \code{lb} The lower band of the regression.
#' where \code{error}, \code{prediction}, \code{ub}, and \code{lb} are \code{fdata} objects.
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
  Y <- nadaraya$data$variables[, nadaraya$predictor]

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

  gd.pred <- list(prediction = res, windows = hs)
  class(gd.pred) <- "bpred"
  return(gd.pred)
}
