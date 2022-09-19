# jewel
**jewel** package with the implementation of [_jewel_ method](https://www.mdpi.com/2227-7390/9/17/2105) for joint estimation of Gaussian grahical model from multiple datasets.


## Installation
To install from github, run
```
#install.packages("devtools")

library("devtools")
install_github("annaplaksienko/jewel")
```

## Example

This is a very minimal working example of all **jewel** package functions. See jewel_manual.pdf for documentation of the package.

```
library(jewel)

K <- 3
n <- 50
p <- 100

#generate the data
data <- generateData(K, n, p)
G_true <- data$trueGraph
X <- data$data

#generate the grid of lambda
lambda <- generateLambdaGrid(X, nlambda = 10)

#estimate regularization parameter with two procedures
lambda_BIC <- estimateLambdaBIC(X, lambda)$lambda_opt
lambda_CV <- estimateLambdaCV(X, lambda)$lambda_opt

#estimate the graph with these two regularization parameters
G_BIC <- jewel(X, lambda = lambda_BIC)$EstAdjMat
G_CV <- jewel(X, lambda = lambda_CV)$EstAdjMat

#evaluate results
edges <- sum(G_true)
noedges <- p * (p - 1) / 2 - sum(G_true)
message("True graph has ", edges, " edges and ", noedges, " absence of edges.")
evaluatePerformance(G = G_true, G_hat = G_BIC)
evaluatePerformance(G = G_true, G_hat = G_CV)
```
