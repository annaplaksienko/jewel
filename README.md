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

This is a minimal working example of all **jewel** package functions. See jewel_manual.pdf for documentation of the package. See also real data example for better understanding.

```
library(jewel)

#First, generate the data: 3 datasets with 100 variables and 50 samples.
#Their underlying graphs will have a quarter of edges not in common (tune with perc parameter).
K <- 3
p <- 100
n <- 50
data <- generateData_rewire(K = K, p = p, n = n)
G_list_true <- data$Graphs
G_common_true <- data$CommonGraph
X <- data$Data

#let's assume we have prior information on whether some degrees are hubs
#to simulate that, we'll choose 3% of vertices with the highest degrees and put their degree to 10 ("hub")
#the rest will be put to 1
#we use only one graph because in simulation the degree distribution is the same by construction for all k
true_degrees <- rowSums(G_list_true[[1]])
cut <- sort(true_degrees, decreasing = TRUE)[ceiling(p * 0.03)]
apriori_hubs <- ifelse(true_degrees >= cut, 10, 1)
#now we use that to construct weights for penalization problem
W <- constructWeights(apriori_hubs, K = K)

#estimate the graphs with user chosen lambda1 and weights W and with stability selection procedure
res <- jewel(X, lambda1 = 0.1, W = W, stability = TRUE)
G_list_est <- res$G_list
G_common_est <- res$CommonG

#now evaluate results
evaluatePerformance(G = G_common_true, G_hat = G_common_est)
mapply(evaluatePerformance, G_list_true, G_list_est)
```

