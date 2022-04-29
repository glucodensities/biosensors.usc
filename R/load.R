## load.R: biosensors.usc glue
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
#' @importFrom utils read.csv
#' @importFrom stats runif approxfun
#' @importFrom truncnorm qtruncnorm

#' @title load_data
#' @description R function to read biosensors data from a csv files.
#' @param filename_fdata A csv file with the functional data. The csv file must have long format with, at least, the following three columns: id, time, and value, where the id identifies the individual, the time indicates the moment in which the data was captured, and the value is a monitor measure.
#' @param filename_variables A csv file with the clinical variables. The csv file contains a row per individual and must have a column id identifying this individual.
#' @return A biosensor object:
#' \code{data} A data frame with biosensor raw data.
#' \code{densities} A functional data object (fdata) with a non-parametric density estimation.
#' \code{quantiles} A functional data object (fdata) with the empirical quantile estimation.
#' \code{variables} A data frame with the covariates.
#' @examples
#' # Data extracted from the paper: Hall, H., Perelman, D., Breschi, A., Limcaoco, P., Kellogg, R.,
#' # McLaughlin, T., Snyder, M., Glucotypes reveal new patterns of glucose dysregulation, PLoS
#' # biology 16(7), 2018.
#' file1 = system.file("extdata", "data_1.csv", package = "biosensors.usc")
#' file2 = system.file("extdata", "variables_1.csv", package = "biosensors.usc")
#' data = load_data(file1, file2)
#' names(data)
#' head(data$quantiles)
#' head(data$variables)
#' plot(data$quantiles, main="Quantile curves")
#' @export
load_data <- function(filename_fdata, filename_variables=NULL) {

  df <- process_data(filename_fdata)
  id_quantiles <- unique(df$id)

  if (!("value" %in% colnames(df)))
    stop("The csv file filename_fdata must have a column named 'value'.")

  if (!("id" %in% colnames(df)))
    stop("The csv file filename_fdata must have a column named 'id'.")

  min_val <- min(df$value, na.rm=TRUE)
  max_val <- max(df$value, na.rm=TRUE)
  t1 <- seq(min_val, max_val, length = 300)
  r1 <- load_density_data(df, t1)
  t2 <- seq(0, 1, length = 300)
  r2 <- load_quantile_data(df, t2)
  r3 <- NULL
  if (!is.null(filename_variables)) {
    r3 <- utils::read.csv(filename_variables)
    ######## ORDENAR VARIABLES IGUAL ORDEN Q DATOS
    id_dataframe = c()
    for (i in id_quantiles) {
      id_dataframe = c(id_dataframe, which(r3$id == i))
    }
    r3 = r3[id_dataframe,]
  }
  data <- list(data = df, densities = r1, quantiles = r2, variables = r3)
  class(data) <- "biosensor"
  return(data)
}

process_data <- function(filename) {
  df <- utils::read.csv(filename)
  digits <- grepl("^[[:digit:]]+", df$value)
  df <- df[digits, ]
  return (df)
}


load_quantile_data <- function(df, t) {
  different <- unique(df$id)
  quantiles_matrix <- matrix(0, ncol = length(t), nrow = length(different))

  counter <- 0
  for (i in different) {
    counter <- counter + 1
    aux <- df[df$id == i, ]
    value <- as.numeric(as.character(aux$value))
    value <- value[!is.na(value)]
    quantiles_matrix[counter, ] <- stats::quantile(value, probs = t)
  }
  return(fda.usc::fdata(quantiles_matrix, argvals = t))
}


load_density_data <- function(df, t) {
  different <- unique(df$id)
  density_matrix <- matrix(0, ncol = length(t), nrow = length(different))
  counter <- 0
  for (i in different) {
    counter <- counter + 1
    aux <- df[df$id == i, ]
    value <- as.numeric(as.character(aux$value))
    value <- value[!is.na(value)]
    density_aux <- stats::density(value, from = min(t), to = max(t))

    approx <- approxfun(density_aux)
    density_aux2 <- approx(t)
    density_aux2[is.na(density_aux2)] <- 0
    density_matrix[counter,] <- density_aux2
  }
  return(fda.usc::fdata(density_matrix, argvals = t))
}


#' @title generate_data
#' @description Generates a quantile regression model V + V2 * v + tau * V3 * Q0 where Q0 is a truncated random variable, v = 2 * X, tau = 2 * X, V ~ Unif(-1, 1), V2 ~ Unif(-1, -1), V3 ~ Unif(0.8, 1.2), and E(V|X) = tau * Q0;
#' @param n Sample size.
#' @param Qp Dimension of the quantile.
#' @param Xp Dimension of covariates where X_i~Unif(0,1).
#' @return A biosensor object:
#' \code{data} NULL.
#' \code{densities} NULL.
#' \code{quantiles} A functional data object (fdata) with the empirical quantile estimation.
#' \code{variables} A data frame with Xp covariates.
#' @examples
#' data = generate_data(n=100, Qp=100, Xp=5)
#' names(data)
#' head(data$quantiles)
#' head(data$variables)
#' plot(data$quantiles, main="Quantile curves")
#' @export
generate_data <- function(n=100, Qp=100, Xp=5) {
  if (n <= 0)
    stop("Error: n must be positive")
  if (Qp <= 0)
    stop("Error: Qp must be positive")
  if (Xp <= 0)
    stop("Error: Xp must be positive")

  Qs <- matrix(0, nrow=n, ncol=Qp)
  teorico <- matrix(0, nrow=n, ncol=Qp)
  X <- matrix(0, nrow=n, ncol=Xp)
  t <- seq(0, 1, length=Qp)
  Q0 <- truncnorm::qtruncnorm(t, -5, 5)

  for(i in 1:n) {
    X[i,] <- runif(Xp)
    v <- sum(2*X[i,])
    tau <- sum(2*X[i,])
    V <- runif(1,-1,1)
    V2 <- runif(1,-1,1)
    V3 <- runif(1,0.8,1.2)
    Qs[i,] <- V+V2*v+tau*V3*Q0
    teorico[i,] <- tau*Q0
  }
  densities = fda.usc::fdata(teorico, argvals = t)
  quantiles = fda.usc::fdata(Qs, argvals = t)

  df = as.data.frame(X)

  data <- list(data = NULL, densities = NULL, quantiles = quantiles, variables = df)
  class(data) <- "biosensor"
  return(data)
}



