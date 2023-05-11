#' Estimate Gaussian graphical models from multiple datasets
#'
#' This function estimates Gaussian graphical models (i.e. networks of conditional dependencies, direct connections between variables) given multiple datasets. 
#' We assume that datasets contain measurements of the same variables collected under different conditions (different equipment, locations, even sub-types of disease).
#'
#' @param X a list of \code{K} numeric data matrices of \code{n_k} samples and 
#' \code{p} variables (\code{n_k} can be different for each matrix).
#' @param lambda1 a number, first regularization parameter (of the common penalty).
#' @param lambda2 an optional number, second regularization parameter 
#' (of the class-specific penalty). If NULL, set to \code{lambda_2 = lambda_1 * 1.4}
#' @param Theta an optional list of \code{K} regression coefficient matrices 
#' of the size \code{p} by \code{p}. User-provided initialization can be used 
#' for warm-start procedures. If NULL, initialized as all zeros.
#' @param W an optional list of \code{K} weights matrices of the size 
#' \code{p} by \code{p}. User-provided initialization can be used when some 
#' vertices are believed to be hubs. If NULL, initialized as all ones.
#' @param tol an optional number, convergence threshold controlling the relative 
#' error between iterations. The default value is 0.01.
#' @param maxIter an optional number, maximum allowed number of iterations. 
#' The default value is 10 000.
#' @param stability if stability = TRUE, stability selection procedure to reduce 
#' the number of false positives will be applied. \code{n_k / 2} samples are 
#' randomly chosen in each dataset \code{stability_nsubsets} times and then 
#' __jewel__ method is applied to each subset. In the final estimate, we include 
#' only the edges that appear in at least \code{stability_frac} proportion of the subsets. 
#' By default this procedure is disabled since it increases the running time.
#' @param stability_nsubsets an optional number, how many times to subsample 
#' datasets and apply __jewel__ for stability selection procedure. The default value is 25.
#' @param stability_frac an optional number, in what proportion of the stability 
#' results on subsampled data an edge has to be present to be included into the 
#' final estimate. The default value is 0.8.
#' @param verbose if verbose = FALSE, tracing information printing is disabled. 
#' The default value is TRUE.
#'
#' @importFrom parallel detectCores makeCluster clusterEvalQ clusterApply stopCluster
#'
#' @return The following list is returned
#' \itemize{
#'   \item \code{CommonG} - an adjacency matrix of the common estimated graph (intersection of \code{K} estimated graphs).
#'   \item \code{G_list} - a list of \code{K} adjacency matrices for each estimated graph.
#'   \item \code{Theta} - a list of \code{K} estimated covariance matrices (when stability selection is disabled).
#'   \item \code{BIC} – a number, value of Bayesian information criterion for resulting graphs (when stability selection is disabled).
#' }
#' 
#' @export
#' @examples
#' {
#' K <- 3
#' p <- 50
#' n <- 20
#' data <- generateData_rewire(K = K, p = p, n = n, ncores = 1, verbose = FALSE)
#' G_list_true <- data$Graphs
#' X <- data$Data
#' true_degrees <- rowSums(G_list_true[[1]])
#' cut <- sort(true_degrees, decreasing = TRUE)[ceiling(p * 0.03)]
#' apriori_hubs <- ifelse(true_degrees >= cut, 10, 1)
#' W <- constructWeights(apriori_hubs, K = K)
#' res <- jewel(X, lambda1 = 0.25, W = W, verbose = FALSE)
#' }

jewel <- function(X, lambda1, lambda2 = NULL, 
                  Theta = NULL, W = NULL,
                  tol = 0.01, maxIter = 10000,
                  stability = FALSE, stability_nsubsets = 25, stability_frac = 0.8,
                  verbose = TRUE) {
  
  if (stability == FALSE) {
    return(jewel_inner(X, lambda1, lambda2, 
                       Theta, W,
                       tol, maxIter, verbose))
  } else {
    
    if (verbose) message("Stability procedure will be run in parallel so printing is disabled. Starting...")
    
    K <- length(X)
    subsets <- vector(mode = "list", length = stability_nsubsets)
    for (s in 1:stability_nsubsets) {
      choose <- lapply(1:K, function(x) sample(1:dim(X[[x]])[1], 
                                               size = dim(X[[x]])[1] / 2))
      subsets[[s]] <- mapply(function(x, y) x[y, ],
                             X, choose, SIMPLIFY = FALSE)
    }
    names(subsets) <- sapply(1:stability_nsubsets, function(i) sprintf("Subset%i", i))
    
    ncores <- detectCores(logical = FALSE) - 1 
    cl <- makeCluster(ncores)
    
    clusterEvalQ(cl, library("jewel"))
    results <- vector(mode = "list", length = stability_nsubsets)
    names(results) <- sapply(1:stability_nsubsets, function(i) sprintf("ResultsSubset%i", i))
    results <- clusterApply(cl, subsets, jewel_inner, 
                            lambda1, lambda2,
                            Theta, W,
                            tol, maxIter, verbose)
    stopCluster(cl)
    
    stab_G_list <- lapply(results, function(x) x$G_list)

    temp_G <- vector(mode = "list", length = K)
    for (k in 1:K) {
      temp_G[[k]] <- lapply(stab_G_list, function(x) x[[k]])
      temp_G[[k]] <- Reduce('+', temp_G[[k]])
    }

    G_list <- vector(mode = "list", length = K)
    G_list <- lapply(temp_G, function(x) x >= stability_nsubsets * stability_frac)
    if(!is.null(names(X))) {
      names(G_list) <- names(X)
    }
    
    CommonG <- Reduce('+', G_list)
    CommonG <- (CommonG == K)
    
    if (verbose) message("Complete.")
    
    return(list(G_list = G_list, CommonG = CommonG))
    }
}