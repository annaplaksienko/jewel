#' Generate a set of scale-free graphs and corresponding datasets (using the graphs as their Gaussian graphical models)
#'
#' Function first generates \code{K} scale-free graphs with \code{p} vertices. They have the same order and degree distribution and share most of the edges, but some edges may vary (user can control how many).
#' Function then generates corresponding precision and covariance matrices, all of the size \code{p} by \code{p} (see the paper for the details of the procedure).
#' Then for each \code{l}-th element of vector \code{n} it generates \code{K} data matrices, each of the size \code{n_l} by \code{p}, 
#' i.e., for the same underlying graphs we can generate several sets of \code{K} datasets with different sample sizes.
#'
#' @param K number of graphs/data matrices.
#' @param p number of nodes in the true graphs.
#' @param n a numerical vector of the sample sizes for each desired set of 
#' \code{K} data matrices. Can be a vector of one element if the user wishes to 
#' obtain only one dataset of \code{K} matrices.
#' @param power a number, power of preferential attachment for the Barabasi-Albert 
#' algorithm for the generation of the scale-free graph. Bigger number means 
#' more connected hubs. The default value is 1.
#' @param m number of edges to add at each step of Barabasi-Albert algorithm 
#' for generation of the scale-free graph. The default value is 1.
#' @param perc a number, tuning parameter for the difference between graphs. 
#' Number of trials to perform in the rewiring procedure of the first graph is 
#' \code{p * perc}. Bigger the number, more different are the graphs.
#' @param int a vector of two numbers, \code{a} and \code{b}. Entries of 
#' precision matrices are sampled from the uniform distribution on the interval 
#' \code{[-b, -a] + [a, b]}. The default values are \code{a = 0.2, b = 0.8}.
#' @param ncores number of cores to use in parallel data generation. 
#' If \code{NULL}, set to \eqn{\#physical cores - 1}.
#' @param makePlot If makePlot = FALSE, plotting of the generated graphs is 
#' disabled. The default value is TRUE.
#' @param verbose If verbose = FALSE, tracing information printing is disabled. 
#' The default value is TRUE.
#'
#' @importFrom igraph barabasi.game as_adjacency_matrix graph_from_adjacency_matrix plot.igraph intersection gsize rewire ecount keeping_degseq
#' @importFrom MASS mvrnorm
#' @importFrom stats runif
#' @importFrom Matrix forceSymmetric
#' @importFrom parallel detectCores mclapply
#' @importFrom purrr transpose
#' @importFrom matrixcalc is.positive.definite
#' 
#' @return The following list is returned
#'  \itemize{
#'    \item \code{Graphs} – a list of adjacency matrices of the \code{K} generated graphs.
#'    \item \code{CommomGraph} - a matrix, common part (intersection) of the \code{K} generated graphs.
#'    \item \code{Data} - a list of lists, for each sample size of the input vector \code{n} one obtains \code{K} data matrices, each of the size \code{n_l} by \code{p}.
#'    \item \code{Sigma} - a list of \code{K} covariance matrices of the size \code{p} by \code{p}.
#' }
#' @export
#' 
#' @examples
#' data <- generateData_rewire(K = 3, p = 50, n = 20, ncores = 1, verbose = FALSE)

generateData_rewire <- function (K, p, n,
                                 power = 1, m = 1, perc = 0.05, 
                                 int = NULL,
                                 ncores = NULL, makePlot = TRUE, verbose = TRUE) {
  
  generation <- function(g) {
    
    #change all non-zero entries to samples from uniform distribution
    #on the interval [d, c] united [a, b]
    
    nn <- ecount(g)
    samp1 <- runif(nn, a, b)
    samp2 <- runif(nn, d, c)
    new_entries <- sample(c(samp1, samp2), size = nn)
    
    #precision matrices
    Omega <- as_adjacency_matrix(g)
    Omega[Omega == 1] <- new_entries
    Omega <- forceSymmetric(Omega)
    
    #ensure positive definiteness
    #diagonal elements = |lambda_min(Omega_k)| + 0.1
    diag(Omega) <- abs(min(eigen(Omega)$values)) + 0.1;
    
    #check if matrices are positive definite
    #package matrixcalc
    SymPosDefCheck <- is.positive.definite(as.matrix(Omega))
    if (SymPosDefCheck == 0) {
      return("Some matrix is not symmmetric positive definite")
    }
    
    #covariance matrices
    OmegaInv <- solve(as.matrix(Omega))
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
  
  if (verbose) message("1/3 Constructing true graphs...")
  
  G_list <- vector(mode = "list", length = K)
  G_list[[1]] <- barabasi.game(n = p, power = power, m = m,
                          directed = FALSE)
  if (makePlot) {
    plot.igraph(G_list[[1]], vertex.label = NA, vertex.size = 1,
                main = "G^(1)")
  }
  

  if (K > 1) {
    for (k in 2:K) {
      G_list[[k]] <- rewire(G_list[[1]], 
                            with = keeping_degseq(niter = gsize(G_list[[1]]) * perc))
    }
  }
  
  G <- as.matrix(as_adjacency_matrix(do.call(intersection, G_list)))
  
  size <- m * p - 2 * m + 1
  common_size <- round(sum(G) / 2)
  diff <- 100 - round(common_size / size * 100)
  
  if (verbose) message("Size of each graph is ", size, " edges.")
  if (verbose) message("Size of the common part is ", common_size, " edges. Difference is ", diff, "%.")
  
  if (verbose) message("2/3 Completed. Constructing the data...")
  
  #boundaries for sampling from uniform distribution
  # d < c < a < b
  if (is.null(int)) {
    a <- 0.2
    b <- 0.8
  } else {
    if (length(int) != 2) {
      stop("int must be a vector of two numbers.")
    } else {
      if (int[1] >= int[2]) {
        stop("Left interval must be less than the right one.")
      } else {
        a <- int[1]
        b <- int[2]
      }
    }
  }
  c <- -a;
  d <- -b;
  
  if (.Platform$OS.type == "windows") {
    ncores <- 1
  } else if (is.null(ncores)) {
    ncores <- detectCores(logical = FALSE) - 1 
  } 
  
  data <- mclapply(G_list, generation, mc.cores = ncores)
  
  X <- transpose(lapply(data, function(x) x$X))
  Sigma <- lapply(data, function(x) x$Sigma)
  remove(data)
  
  G_list <- lapply(G_list, function(x) as.matrix(as_adjacency_matrix(x)))
  
  if (length(n) == 1) {
    X <- X[[1]]
    names(X) <- sapply(1:K, function(i) sprintf("X%i", i))
  }
  
  names(Sigma) <- sapply(1:K, function(i) sprintf("Sigma%i", i))
  names(G_list) <- sapply(1:K, function(i) sprintf("G%i", i))
  
  if (verbose) message("3/3 Completed.")
  return(list(Graphs = G_list,
              CommonGraph = G, 
              Data = X, 
              Sigma = Sigma))
}
