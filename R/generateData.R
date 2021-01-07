#' Generate a scale-free graph and corresponding datasets using it as their Gaussian graphical model
#'
#' Function generates a scale-free graph with \code{p} vertices, \code{K} corresponding precision and covariance matrices, all of dize \code{p} by \code{p}
#' and then for each \code{l}-th element of vector \code{n} generates \code{K} data matrices of size \code{n_l} by \code{p}. For the same underlying graph we can generate several datasets with different sample sizes.
#'
#' @param K Number of data matrices
#' @param n Vector of the sample sizes for each desired set of data matrices. Can be a vector of one element.
#' @param p Number of nodes in the true graph.
#' @param power Power of preferential attachment for Barabasi-Albert algorithm for generation of scale-free graph.
#' @param m Number of edges to add in each time step of Barabasi-Albert algorithm for generation of scale-free graph.
#' @param c Entries in precision matrices are generated from the uniform distribution on the interval \code{[-d, -c] + [c, d]}. The default value is \code{c = 0.2}.
#' @param d Entries in precision matrices are generated from the uniform distribution on the interval \code{[-d, -c] + [c, d]}. The default value is \code{d = 0.8}.
#' @param makePlot If makePlot = FALSE, plotting of the generated true graph is disabled. The default value is TRUE.
#' @param verbose If verbose = FALSE, tracing information printing is disabled. The default value is TRUE.
#'
#' @importFrom igraph barabasi.game as_adjacency_matrix graph_from_adjacency_matrix plot.igraph
#' @importFrom MASS mvrnorm
#' @importFrom stats runif
#' @importFrom Matrix forceSymmetric
#' @importFrom parallel detectCores mclapply
#' @importFrom purrr transpose
#' @importFrom matrixcalc is.positive.definite
#' @importFrom rlist list.flatten

#'
#' @return The following list is returned
#'  \itemize{
#'    \item trueGraph - sparse adjacency matrix of the true graph
#'    \item data - list of lists, \code{K} data matrices of the size \code{n_l} by \code{p} for each sample size element \code{l} of the input vector \code{n}
#'    \item Sigma - list of \code{K} covariance matrices of the size \code{p} by \code{p}
#' }
#' @export

generateData <- function (K, n, p,
                          power = 1, m = 1, c = 0.2, d = 0.8,
                          makePlot = TRUE, verbose = TRUE) {

  generation <- function(k) {

    #change all non-zero entries to samples from uniform distribution
    #on the interval [a, b] united [c, d]

    samp1 <- runif(nn / 2, a, b)
    samp2 <- runif(nn / 2, c, d)
    new_entries <- sample(c(samp1, samp2))

    #precision matrices
    Omega <- G
    Omega[Omega == 1] <- new_entries
    Omega <- forceSymmetric(Omega)

    #ensure positive definiteness
    #diagonal elements = |lambda_min(Omega_k)| + 0.1
    diag(Omega) <- abs(min(eigen(Omega)$values)) + 0.1;

    #check if matrices are positive definite
    #package matrixcalc
    SymPosDefCheck <- is.positive.definite(as.matrix(Omega))
    if (SymPosDefCheck == 0) {
      print("Some matrix is not symmmetric positive definite")
      return("Some matrix is not symmmetric positive definite");
    }

    #covariance matrices
    OmegaInv <- chol2inv(as.matrix(Omega))
    Sigma <- OmegaInv;

    for (i in 1:p) {
      for (j in 1:i) {
        Sigma[i, j] <- OmegaInv[i, j] / sqrt(OmegaInv[i, i] * OmegaInv[j, j]);
        Sigma[j, i] <- OmegaInv[j, i] / sqrt(OmegaInv[i, i] * OmegaInv[j, j]);
      }
    }

    #generate simulation data
    #package MASS
    X <- vector(mode = "list", length = length(n))
    for (i in 1:length(n)) {
      X[[i]]  <- scale(mvrnorm(n = n[i], mu = rep(0, p),
                               Sigma = Sigma, empirical = FALSE))
    }

    names(X) <- paste("X", n, sep = "_")

    return(list(X = X, Sigma = Sigma, Omega = Omega))

  }

  if (verbose) message("1/3 Constructing true graph...")

  G <- as_adjacency_matrix(barabasi.game(n = p,
                                         power = power, m = m,
                                           directed = FALSE))

  if (makePlot == TRUE) plot.igraph(graph_from_adjacency_matrix(G, mode = "undirected"),
                                      vertex.label = NA)

  nn <- sum(G)

  if (verbose) message("2/3 Completed. Constructing precision matrices...")

  #boundaries for sampling from uniform distribution
  # a < b < c< d
  a <- -d;
  b <- -c;


  if (.Platform$OS.type == "windows") {
    ncores <-  1
  } else ncores <- detectCores() - 1


  data <- mclapply(c(1:K), generation, mc.cores = ncores)

  X <- transpose(lapply(data, function(x) x$X))
  Sigma <- lapply(data, function(x) x$Sigma)
  Omega <- lapply(data, function(x) x$Omega)
  remove(data)

  if (verbose) message("3/3 Completed.")
  
  if (length(n) == 1) {
    X <- X[[1]]
  }

  return(list(trueGraph = G, data = X, Sigma = Sigma))
}
