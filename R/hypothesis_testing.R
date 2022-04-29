## comparison.R: biosensors.usc glue
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

#' @importFrom fda.usc fdata func.mean func.var semimetric.basis metric.lp
#' @importFrom graphics par plot lines legend
#' @importFrom stats pchisq
#' @importFrom methods is

#' @title hypothesis_testing
#' @description Hypothesis testing between two random samples of distributional representations to detect differences in scale and localization (ANOVA test) or distributional differences (Energy distance).
#' @param data1 A biosensor object. First population.
#' @param data2 A biosensor object. Second population.
#' @param permutations Number of permutations used in the energy distance calibration test.
#' @return An object of class biotest:
#' \code{p1_mean} Quantile mean of the first population.
#' \code{p1_variance} Quantile variance of the first population.
#' \code{p2_mean} Quantile mean of the second population.
#' \code{p2_variance} Quantile variance of the second population.
#' \code{energy_pvalue} P-value of the energy distance test.
#' \code{anova_pvalue} P-value of the ANOvA-Fréchet test.
#' @usage
#' hypothesis_testing(data1, data2, permutations=100)
#' @examples
#' # Data extracted from the paper: Hall, H., Perelman, D., Breschi, A., Limcaoco, P., Kellogg, R.,
#' # McLaughlin, T., Snyder, M., Glucotypes reveal new patterns of glucose dysregulation, PLoS
#' # biology 16(7), 2018.
#' file1 = system.file("extdata", "data_1.csv", package = "biosensors.usc")
#' file2 = system.file("extdata", "variables_1.csv", package = "biosensors.usc")
#' data1 = load_data(file1, file2)
#' file3 = system.file("extdata", "data_2.csv", package = "biosensors.usc")
#' file4 = system.file("extdata", "variables_2.csv", package = "biosensors.usc")
#' data2 = load_data(file3, file4)
#' htest = hypothesis_testing(data1, data2)
#' @export
hypothesis_testing <- function(data1, data2, permutations = 100) {
  if (!is(data1, "biosensor"))
    stop("Error: data1 must be an object of biosensor class. @seealso biosensors.usc::load_data")

  if (!is(data2, "biosensor"))
    stop("Error: data2 must be an object of biosensor class. @seealso biosensors.usc::load_data")

  if (!is(data1$quantiles, "fdata"))
    stop("The data attribute quantiles must be of type fdata in data1")

  if (!is(data1$variables, "data.frame"))
    stop("The data attribute quantiles must be of type fdata in data1")

  if (!is(data2$quantiles, "fdata"))
    stop("The data attribute quantiles must be of type fdata in data2")

  if (!is(data1$variables, "data.frame"))
    stop("The data attribute quantiles must be of type fdata in data2")

  m1 <- fda.usc::func.mean(data1$quantiles)
  m2 <- fda.usc::func.mean(data2$quantiles)
  r1 <- range(m1)
  r2 <- range(m2)
  minys <- min(c(r1[1], r2[1]))
  maxys <- max(c(r1[2], r2[2]))

  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))                 # code line i + 1

  graphics::par(mfrow = c(1, 1))

  graphics::plot(m1, main = "Quantile mean", xlab = "t", ylab = "X(t)", lwd = 1, col = "blue", ylim = c(minys, maxys))
  graphics::lines(m2, lwd = 1, col = "red")
  graphics::legend("topleft", legend = c("Population 1", "Population 2"), col = c("blue", "red"), lty = 1:2, cex = 0.8)

  v1 <- fda.usc::func.var(data1$quantiles)
  v2 <- fda.usc::func.var(data2$quantiles)
  r1 <- range(v1)
  r2 <- range(v2)
  minys <- min(c(r1[1], r2[1]))
  maxys <- max(c(r1[2], r2[2]))

  graphics::par(mfrow = c(1, 1))
  graphics::plot(v1,
    main = "Quantile variance", xlab = "t", ylab = "X(t)",
    lwd = 1, col = "blue", ylim = c(minys, maxys)
  )
  graphics::lines(v2, lwd = 1, col = "red")
  graphics::legend("topleft", legend = c("Population 1", "Population 2"), col = c("blue", "red"), lty = 1:2, cex = 0.8)

  # Energy
  epv <- 0
  tryCatch(
    {
      epv <- calcular_pvalor_energia(data1$quantiles, data2$quantiles, as.numeric(permutations))
    },
    error = function(e) {
      message("An error occured while computing the p-value between both populations:\n", e)
    }
  )

  # Anova
  apv <- 0
  tryCatch(
    {
      apv <- calcular_pvalor_anova(data1$quantiles, data2$quantiles)
    },
    error = function(e) {
      message("An error occured while computing the p-value between both populations:\n", e)
    }
  )



  test <- list(
    p1_mean = m1, p1_variance = v1, p2_mean = m2, p2_variance = v2,
    energy_pvalue = epv, anova_pvalue = apv
  )
  class(test) <- "biotest"
  return(test)
}


calcular_pvalor <- function(D, E, n, m, w1, w2, nrep) {
  accu <- 0
  e <- 1:(n + m)
  for (i in 1:nrep) {
    p1 <- sample(e, n)
    p2 <- e [!e %in% p1]
    # Hay q rehacer las matrices A, B, y C a partir de los valores del sample
    A1 <- matrix(0, n, m)
    B1 <- matrix(0, n, n)
    C1 <- matrix(0, m, m)
    for (j in 1:n) {
      i1 <- p1[j]
      for (k in 1:m) {
        i2 <- p2[k]
        A1[j, k] <- D[i1, i2]
      }
      for (k in 1:n) {
        i2 <- p1[k]
        B1[j, k] <- D[i1, i2]
      }
    }
    for (j in 1:m) {
      i1 <- p2[j]
      for (k in 1:m) {
        i2 <- p2[k]
        C1[j, k] <- D[i1, i2]
      }
    }
    E1 <- 2 * t(w1) %*% A1 %*% w2 - t(w1) %*% B1 %*% w1 - t(w2) %*% C1 %*% w2
    if (E1 > E) {
      accu <- accu + 1
    }
  }
  # return (1/(accu + 1))
  return(accu / (nrep + 1))
}


calcular_pvalor_energia <- function(fun1, fun2, nrep) {
  A <- fda.usc::semimetric.basis(fun1, fun2,
    # type.basis1="fourier", nbasis1=11, type.basis2="fourier", nbasis2=11)
    nbasis1 = 4, nbasis2 = 4
  )
  B <- fda.usc::semimetric.basis(fun1, fun1,
    # type.basis1="fourier", nbasis1=11, type.basis2="fourier", nbasis2=11)
    nbasis1 = 4, nbasis2 = 4
  )
  C <- fda.usc::semimetric.basis(fun2, fun2,
    # type.basis1="fourier", nbasis1=11, type.basis2="fourier", nbasis2=11)
    nbasis1 = 4, nbasis2 = 4
  )

  n <- nrow(fun1)
  m <- nrow(fun2)
  w1 <- rep(1 / n, n)
  w2 <- rep(1 / m, m)

  E <- 2 * t(w1) %*% A %*% w2 - t(w1) %*% B %*% w1 - t(w2) %*% C %*% w2
  D <- rbind(cbind(B, A), cbind(t(A), C))

  return(calcular_pvalor(D, E, n, m, w1, w2, nrep))
}







calcular_pvalor_anova <- function(q1, q2) {
  n1 <- dim(q1)[1]
  n2 <- dim(q2)[1]
  n <- n1 + n2

  m1 <- dim(q1)[2]
  m2 <- dim(q2)[2]
  if (m1 != m2) {
    stop("The second dimension of q1 and q2 must be the same")
  }

  gri <- seq(0, 1, length = m1)

  lambda_1 <- n1 / n
  lambda_2 <- n2 / n

  data <- fda.usc::fdata(rbind(q1$data, q2$data), argvals = gri)

  mean_1 <- fda.usc::func.mean(q1)
  mean_2 <- fda.usc::func.mean(q2)
  mean_global <- fda.usc::func.mean(data)

  distance_1 <- fda.usc::metric.lp(data$data, mean_global$data, p = 2)
  distance_2 <- fda.usc::metric.lp(q1$data, mean_1$data, p = 2)
  distance_3 <- fda.usc::metric.lp(q2$data, mean_2$data, p = 2)


  dispersion_global <- (1 / n) * sum(distance_1^2)
  dispersion_population_1 <- (1 / n1) * sum(distance_2^2)
  dispersion_population_2 <- (1 / n2) * sum(distance_3^2)
  dispersion_total <- dispersion_global - lambda_1 * dispersion_population_1 - lambda_2 * dispersion_population_2

  variance_1 <- (1 / n1) * sum(distance_2^4) - ((1 / n1) * sum(distance_2^2))^2
  variance_2 <- (1 / n2) * sum(distance_3^4) - ((1 / n2) * sum(distance_3^2))^2
  variance_global <- ((lambda_1 * lambda_2) / (variance_1 * variance_2)) * (dispersion_population_1 - dispersion_population_2)^2

  statistical_global <- (n / (lambda_1 / variance_1 + lambda_2 / variance_2)) * variance_global + (n * dispersion_total^2) / (lambda_1 * variance_1 + lambda_2 * variance_2)

  # print(statistical_global)
  # print(statistical_global)
  pvalue <- 1 - stats::pchisq(statistical_global, 1)

  return(pvalue)
}
