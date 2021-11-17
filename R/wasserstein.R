## wass_regression.R: biosensors.usc glue
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

#' @importFrom fda.usc fdata
#' @importFrom graphics par plot lines
#' @importFrom stats complete.cases
#'
#' @title wasserstein_analysis
#' @description Performs the Wasserstein regression analysis.
#' @param data A gd object.
#' @param predictor The name of the vector of observed values (one of the columns of data$variables).
#' @return An object of class wasserstein containing the components:
#' \code{error} The residuals.
#' \code{prediction} The fitted regression.
#' \code{ub} The upper band of the regression.
#' \code{lb} The lower band of the regression.
#' where \code{error}, \code{prediction}, \code{ub}, and \code{lb} are \code{fdata} objects.
#' @export
wasserstein_regression <- function(data, predictor) {
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

  wass <- wasserstein(data, predictor)
  band <- confidence_band(data, predictor)

  Qp <- fda.usc::fdata(band$Qpred, argvals = band$t)
  Ql <- fda.usc::fdata(band$Q_lx, argvals = band$t)
  Qu <- fda.usc::fdata(band$Q_ux, argvals = band$t)

  graphics::par(mfrow = c(1, 1))
  graphics::plot(wass$error, main = "Residual functions vs fitted values", xlab = "t", ylab = "X(t)")
  graphics::par(mfrow = c(1, 1))
  graphics::plot(Qp, main = "Confidence band of mean values", xlab = "t", ylab = "X(t)")
  graphics::lines(Qu, lty = 2, col = "red")
  graphics::lines(Ql, lty = 2, col = "red")

  gd.wasserstein <- list(prediction = Qp, ub = Qu, lb = Ql, regression = wass, data = data, predictor = predictor)
  class(gd.wasserstein) <- "bwasserstein"
  return(gd.wasserstein)
}


#' @title wasserstein_prediction
#' @description Performs the Wasserstein regression prediction.
#' @param data A gd object.
#' @param predictor The name of the vector of observed values (one of the columns of data$variables).
#' @return An object of class wasserstein containing the components:
#' \code{error} The residuals.
#' \code{prediction} The fitted regression.
#' \code{ub} The upper band of the regression.
#' \code{lb} The lower band of the regression.
#' where \code{error}, \code{prediction}, \code{ub}, and \code{lb} are \code{fdata} objects.
#' @export
wasserstein_prediction <- function(reg, xpred) {
  if (class(reg) != "bwasserstein")
    stop("Error: data must be an object of bwasserstein class. ")

  object <- cpp_wasserstein_regression(reg$regression$xfit, reg$regression$q, reg$regression$Q0, xpred,
                                       reg$regression$t, reg$regression$qdmin)

  plot(fdata(object$Qpred))
  return(object$Qpred)
}




wasserstein <- function(data, predictor) {
  nas <- tryCatch(
    {
      !is.na(data$variables[, predictor])
    },
    error = function(e) {
      message("An error occured while computing the wassertein regression:\n", e)
    }
  )
  y <- data$variables[nas, predictor]

  real <- data$quantiles
  real$data <- real$data[nas, ]

  t <- seq(0, 1, length = ncol(real$data))
  h <- t[2]

  xfit <- as.matrix(y)
  q <- derivative(real$data, h)
  Q0 <- as.matrix(real$data[, 1])
  xpred <- t(as.matrix(c(mean(xfit))))
  qdmin <- 1e-6

  object <- cpp_wasserstein_regression(xfit, q, Q0, xpred, t, qdmin)

  predicho <- fda.usc::fdata(object$Qfit, argvals = t)
  error <- real - predicho

  salida = list(
    "q"       = q,
    "Q0"      = Q0,
    "t"       = t,
    "qdmin"   = qdmin,
    "xfit"    = object$xfit,
    "xpred"   = object$xpred,
    "Qfit"    = object$Qfit,
    "Qpred"   = object$Qpred,
    "qfit"    = object$qfit,
    "qpred"   = object$qpred,
    "ffit"    = object$ffit,
    "fpred"   = object$fpred,
    "error"   = error
  )

  return(salida)
}


confidence_band <- function(data, predictor) {
  nas <- tryCatch(
    {
      !is.na(data$variables[, predictor])
    },
    error = function(e) {
      message("An error occured while computing the confidence bands of the wassertein regression:\n", e)
    }
  )
  y <- data$variables[nas, predictor]

  real <- data$quantiles
  real$data <- real$data[nas, ]

  t <- seq(0, 1, length = ncol(real$data))
  h <- t[2]

  xfit <- as.matrix(y)
  q <- derivative(real$data, h)
  xpred <- t(as.matrix(c(mean(xfit))))
  q0_obs <- as.matrix(q)
  Q0_obs <- as.matrix(real$data)

  return(cpp_confidence_band(xfit, xpred, Q0_obs, q0_obs, t, 0.05))
}


derivative_vector <- function(datos, h) {
  if (!is.numeric(datos)) {
    stop("Parameter datos must be numeric")
  }
  resultado <- c()
  resultado[1] <- (-3 * datos[1] + 4 * datos[2] - datos[3]) / (2 * h)
  for (i in 2:(length(datos) - 1)) {
    resultado[i] <- (datos[i + 1] - datos[i - 1]) / (2 * h)
  }
  final <- length(datos)
  resultado[final] <- (datos[final - 2] - 4 * datos[final - 1] + 3 * datos[final]) / (2 * h)
  return(resultado)
}

derivative <- function(ds, h) {
  resultado <- c()
  for (i in 1:nrow(ds)) {
    resultado <- c(resultado, derivative_vector(as.numeric(ds[i, ]), h))
  }
  return(matrix(nrow = nrow(ds), ncol = ncol(ds), data = t(resultado), byrow = TRUE))
}
