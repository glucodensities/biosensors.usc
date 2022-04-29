 <!-- badges: start -->
  [![R-CMD-check](https://github.com/glucodensities/biosensors.usc/workflows/R-CMD-check/badge.svg)](https://github.com/glucodensities/biosensors.usc/actions)
  <!-- badges: end -->

# biosensors.usc

New distributional representations of biosensors data in different statistical modeling tasks. 
Please use reference [1] to cite this package.


## References

1. Matabuena, M., Petersen, A., Vidal, J. C., & Gude, F. (2021). Glucodensities: A new representation of glucose profiles using distributional data analysis. Statistical methods in medical research, 0962280221998064.

2. Matabuena, M., & Petersen, A. (2021). Distributional data analysis with accelerometer data in a NHANES database with nonparametric survey regression models. arXiv preprint arXiv:2104.01165.


## Abstract

The biosensor.usc aims to provide a unified and user-friendly framework for using new distributional representations of biosensors data in different statistical modeling tasks: regression models, hypothesis testing, cluster analysis, visualization, and descriptive analysis.
Distributional representations are a functional extension of compositional time-range metrics and we have used them successfully so far in modeling glucose profiles and accelerometer data. However, these functional representations can be used to represent any biosensor data such as ECG or medical imaging such as fMRI.

## Installation Instructions

### Required software and packages
    
1. R (https://www.r-project.org/)

2. R packages: 
[Rcpp](https://CRAN.R-project.org/package=Rcpp), 
[RcppArmadillo](https://CRAN.R-project.org/package=RcppArmadillo),  
[energy](https://CRAN.R-project.org/package=energy), 
[fda.usc](https://CRAN.R-project.org/package=fda.usc), 
[osqp](https://CRAN.R-project.org/package=osqp), 
[truncnorm](https://CRAN.R-project.org/package=truncnorm), 
[parallelDist](https://CRAN.R-project.org/package=parallelDist),
graphics, stats, methods, utils (required in R >= 2.14).

Please install the required R packages before you install the biosensor.usc package. After the installation of the dependencies, please install the **biosensor.usc** as following steps.

### Install biosensor.usc from source code

Install from source code using devtools library:

```
library("devtools")
install_github("glucodensities/biosensors.usc@main")
```

## Usage Instructions

biosensor.usc is an R package which provides:

1) Loading biosensors data from a csv files. 

2) Generating a quantile regression model V + V2 * v + tau * V3 * Q0 where Q0 is a truncated random variable, v = 2 * X, tau = 2 * X, V ~ Unif(-1, 1), V2 ~ Unif(-1, -1), V3 ~ Unif(0.8, 1.2), and E(V|X) = tau * Q0.

3) Performing a Wasserstein regression using a quantile density function.

4) Performing a prediction from a Wasserstein regression.

5) Performing a Ridge regression using a quantile density function.

6) Performing a functional non-parametric Nadaraya-Watson regression with 2-Wasserstein distance, using as predictor the distributional representation and as response a scalar outcome.

7) Performing a prediction from a functional non-parametric Nadaraya-Watson regression with 2-Wasserstein distance.

8) Performing a hypothesis testing between two random samples of distributional representations to detect differences in scale and localization (ANOVA test) or distributional differences (Energy distance).

9) Performing a energy clustering with Wasserstein distance using quantile distributional representations as covariates.

10) Obtaining the clusters to which a set of object belong using a previously trained energy clustering with Wasserstein distance using quantile distributional representations as covariates.


The following codes show how to call above steps in R.

We also attach a data set example through csv files in the package, extracted from the paper: Hall, H., Perelman, D., Breschi, A., Limcaoco, P., Kellogg, R., McLaughlin, T., Snyder, M., “Glucotypes reveal new patterns of glucose dysregulation”, PLoS biology 16(7), 2018.

This data set has two different types of files. The first one contains the functional data, which csv files must have long format with, at least, the following three columns: id, time, and value, where the id identifies the individual, the time indicates the moment in which the data was captured, and the value is a monitor measure:

```
library(biosensors.usc)
file1 = system.file("extdata", "data_1.csv", package = "biosensors.usc")
```

The second type contains the clinical variables. This csv file must contain a row per individual and must have a column id identifying this individual.

```
file2 = system.file("extdata", "variables_1.csv", package = "biosensors.usc")
```

From these files, biosensor data can be loaded as follow: 

```
data1 = load_data(file1, file2)
```

We also provide a way to generate biosensor data from the aforementioned quantile regression model:

```
data1 = generate_data(n=100, Qp=100, Xp=5)
```

Call the Wasserstein regression :

```
wass = wasserstein_regression(data1, "BMI")
```

Use the previously computed Wasserstein regression to obtain the regression prediction given a kxp matrix of input values, where k is the number of points we do the prediction and p is the dimension of the input variables:

```
xpred = as.matrix(25)
pred = wasserstein_prediction(wass, xpred)
```

Alternatively we can also compute the Wasserstein regression using the following function: 

```
wass = regmod_regression(data1, "BMI")
```

Call the Ridge regression:

```
regm = ridge_regression(data1, "BMI")
```

Call the Nadaraya-Watson regression with 2-Wasserstein distance:

```
 nada = nadayara_regression(data1, "BMI")
```

Use the previously computed Nadaraya-Watson regression to obtain the regression prediction given the quantile curves:

```
npre = nadayara_prediction(nada, t(colMeans(data1$quantiles$data)))
```


Call the Hypothesis testing between two random samples of distributional representations:

```
file3 = system.file("extdata", "data_2.csv", package = "biosensors.usc")
file4 = system.file("extdata", "variables_2.csv", package = "biosensors.usc")
data2 = load_data(file3, file4)
htest = hypothesis_testing(data1, data2)
```

Call the energy clustering with Wasserstein distance using quantile distributional representations as covariates:

```
clus = clustering(data1, clusters=3)
```


Use the previously computed clustering to obtain the clusters of the given objects: 

```
assignments = clustering_prediction(clus, data1$quantiles$data)
```


