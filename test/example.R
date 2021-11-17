
library(biosensors.usc)

setwd("./test")

df <- utils::read.csv("data_1.csv")

g1 = biosensors.usc::load_data("data_1.csv", "variables_1.csv")

g1w = biosensors.usc::wasserstein_regression(g1, "BMI")

g1r = biosensors.usc::nadayara_regression(g1, "BMI")

g2 = biosensors.usc::load_data("data_2.csv", "variables_2.csv")

g1c = biosensors.usc::clustering(g1, clusters = 3, iter_max = 10, restarts = 1)

g1g2 = biosensors.usc::hypothesis_testing(g1, g2)


g1 = generate_data(n =100, Qp = 100, Xp = 1)
g2 = generate_data(n =100, Qp = 100, Xp = 1)

g1rm = regmod_regression(g1, "V1")
g1rmp = regmod_prediction(g1rm, g1$quantiles$data)

g1w = biosensors.usc::wasserstein_regression(g1, "V1")

xpred = t(as.matrix(c(mean(g1$quantiles$data[1,]))))
g1wp = biosensors.usc::wasserstein_prediction(g1w, xpred)
xpred = t(as.matrix(c(mean(g1$quantiles$data[2,]))))
g1wp = biosensors.usc::wasserstein_prediction(g1w, xpred)
xpred = t(as.matrix(c(mean(g1$quantiles$data[3,]))))
g1wp = biosensors.usc::wasserstein_prediction(g1w, xpred)
xpred = t(as.matrix(c(mean(g1$quantiles$data[4,]))))
g1wp = biosensors.usc::wasserstein_prediction(g1w, xpred)
xpred = t(as.matrix(c(mean(g1$quantiles$data[5,]))))
g1wp = biosensors.usc::wasserstein_prediction(g1w, xpred)

g1n = biosensors.usc::nadayara_regression(g1, "V1")
g1np = biosensors.usc::nadayara_prediction(g1n, g1$quantiles$data)
plot(g1n$nadayara$coefficients[1,], g1np$prediction[1,])
g1np = biosensors.usc::nadayara_prediction(g1n, t(colMeans(g1$quantiles$data)))
plot(g1n$nadayara$coefficients[1,], g1np$prediction[1,])



g1c = biosensors.usc::clustering(g1, clusters = 3, iter_max = 10, restarts = 1)
g1cp = biosensors.usc::clustering_prediction(g1c, g1$quantiles$data)


g1g2 = biosensors.usc::hypothesis_testing(g1, g2)



###
rr = biosensors.usc::ridge_regression(g1, "V1")

n=1000
p= 2
X= matrix(runif(n*p,0,1), nrow= n, ncol=p)
Y= X[,1]+rnorm(n)
m=ridgesurvey(X, Y,w= rep(1,n))


