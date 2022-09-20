# jewel
**jewel** package with the implementation ofthe updated version of the [_jewel_ method](https://www.mdpi.com/2227-7390/9/17/2105) for joint estimation of Gaussian graphical models (i.e. networks of conditional dependencies, direct connections between variables) given several datasets. We assume that datasets contain measurements of the same variables collected under different conditions (different equipment, locations, even sub-types of disease).


## Installation
To install from github, run
```
#install.packages("devtools")

library("devtools")
install_github("annaplaksienko/jewel")
```

## Example

This is a minimal working example of all **jewel** package functions. See jewel_manual.pdf for documentation of the package.

```
library(jewel)

#First, generate the data: 3 datasets with 100 variables and 50 samples.
#Their underlying graphs will have a quarter of edges not in common (tune with perc parameter).
data <- generateData_rewire(K = 3, p = 100, n = 50, perc = 0.08)
G_list_true <- data$Graphs
G_common_true <- data$CommonGraph
X <- data$Data

#estimate regularization parameter with Bayesian information criterion
#note that this part can be time consuming
#if you want to skip this step – use any value of lambda as "optimal" in the next steps
lambda_BIC <- estimateLambdaBIC(X)$lambda_opt

#estimate the graphs with lambda_BIC
res_BIC <- jewel(X, lambda1 = lambda_BIC)
G_list_est <- res_BIC$G_list
G_common_est <- res_BIC$CommonG

#now evaluate results
evaluatePerformance(G = G_common_true, G_hat = G_common_est)
mapply(evaluatePerformance, G_list_true, G_list_est)
```

