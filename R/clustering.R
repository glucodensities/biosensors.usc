## clustering.R: biosensors.usc glue
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

#' @importFrom energy kgroups
#' @importFrom graphics par plot
#' @importFrom stats dist
#' @importFrom methods is

#' @title clustering
#' @description Performs energy clustering with Wasserstein distance using quantile distributional representations as covariates.
#' @param data A biosensor object.
#' @param clusters Number of clusters.
#' @param iter_max Maximum number of iterations.
#' @param restarts Number of restarts.
#' @return An object of class bclustering:
#' \code{data} A data frame with biosensor raw data.
#' \code{result} A kgroups object (see energy library).
#' @usage
#' clustering(data, clusters=3, iter_max=10, restarts=1)
#' @examples
#' # Data extracted from the paper: Hall, H., Perelman, D., Breschi, A., Limcaoco, P., Kellogg, R.,
#' # McLaughlin, T., Snyder, M., Glucotypes reveal new patterns of glucose dysregulation, PLoS
#' # biology 16(7), 2018.
#' file1 = system.file("extdata", "data_1.csv", package = "biosensors.usc")
#' file2 = system.file("extdata", "variables_1.csv", package = "biosensors.usc")
#' data = load_data(file1, file2)
#' clus = clustering(data, clusters=3)
#' @export
clustering <- function(data, clusters = 3, iter_max = 10, restarts = 1) {
  if (!is(data, "biosensor"))
    stop("Error: data must be an object of biosensor class. @seealso biosensors.usc::load_data")

  if (is.null(data$quantiles))
    stop("The data attribute quantiles can not be NULL")

  if (!is(data$quantiles, "fdata"))
    stop("The data attribute quantiles must be of type fdata")

  if (!is.null(data$densities) && class(data$densities) != "fdata")
    stop("The data attribute densities must be of type fdata")


  result <- tryCatch(
    {
      energy::kgroups(data$quantiles$data, clusters, iter.max = iter_max, nstart = restarts)
    },
    warning = function(w) {
      message("A warning occured while clustering the data:\n", w)
    },
    error = function(e) {
      message("An error occured while clustering the data:\n", e)
    }
  )

  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1

  if (is.null(data$densities))
    graphics::par(mfrow = c(result$k, 1))
  else
    graphics::par(mfrow = c(result$k, 2))
  for (i in 1:result$k) {
    graphics::plot(data$quantiles[result$cluster == i],
      main = paste("Cluster ", i, " (quantiles)"),
      xlab = "t", ylab = "X(t)", lwd = 1
    )
    if (!is.null(data$densities))
      graphics::plot(data$densities[result$cluster == i],
        main = paste("Cluster ", i, " (densities)"),
        xlab = "t", ylab = "X(t)", lwd = 1
      )
  }

  bclustering = list(
    "data"   = data,
    "result" = result
  )

  class(bclustering) <- "bclustering"
  return(bclustering)
}



#' @title clustering_prediction
#' @description Predicts the cluster of each element of the objects list
#' @param clustering A gl.clustering object.
#' @param objects Matrix of objects to cluster.
#' @return The clusters to which these objects are assigned.
#' @usage
#' clustering_prediction(clustering, objects)
#' @export
clustering_prediction <- function(clustering, objects) {
  if (!is(clustering, "bclustering"))
    stop("Error: data must be an object of bclustering class. ")

  X = as.matrix(clustering$data$quantiles$data)
  np= length(unique(clustering$result$cluster))
  n= dim(objects)[1]
  nx= dim(X)[1]
  dispersiones= matrix(0, nrow= n, ncol=np)

  for(i in 1:n) {
    for(j in 1:np) {
      clusteraux = X[clustering$result$cluster==j,]
      dataux = t(data.frame(objects[i,]))
      dataux = as.matrix(dataux)
      juntar = rbind(dataux, clusteraux)
      juntar = as.matrix(juntar)
      juntar = dist(juntar)
      juntar = as.matrix(juntar)
      juntaraux= juntar[1,-c(1)]
      juntar2= juntar[-c(1),-c(1)]
      suma= 2*mean(juntaraux)
      nxc= dim(clusteraux)[1]

      suma2= (1/(nxc))^2*sum(juntar2)
      dispersiones[i,j]= suma-suma2
    }
  }

  asignacion= apply(dispersiones, 1, which.min)

  return(asignacion)
}


