#' Construct weights for _jewel_ minimization problem from prior information on vertices degrees.
#'
#' Function takes a numerical vector of vertices degrees and constructs weights with the rule \code{W_ij = 1 / sqrt(d_i * d_j)} and then the whole matrix is normilized by the maximum. 
#
#'
#' @param d either one numerical vector or a list of \code{K} numerical vectors 
#' of length \code{p} with user-provided degrees of vertices for each class. 
#' If there is only one vector, we assume degrees are the same for all \code{K} classes. 
#' In that case parameter \code{K} (number of classes) must be provided.
#' Note that for successful _jewel_ estimation true degrees are not necessary: 
#' for example, user can provide a vector where known hubs have degree 10 
#' and the rest of the vertices have degree 1.
#' @param K number of classes (i.e. datasets, i.e. desired graphs). By default it is length(d).
#' In length(d) = 1, \code{K} must be provided by the user.
#'
#' @return W - a list of \code{K} numeric matrices of the size \code{p} by \code{p}
#'
#' @export
#' 
#' @examples
#' {
#' K <- 3
#' p <- 50
#' n <- 20
#' data <- generateData_rewire(K = K, p = p, n = n, ncores = 1, verbose = FALSE)
#' G_list_true <- data$Graphs
#' true_degrees <- rowSums(G_list_true[[1]])
#' cut <- sort(true_degrees, decreasing = TRUE)[ceiling(p * 0.03)]
#' apriori_hubs <- ifelse(true_degrees >= cut, 10, 1)
#' W <- constructWeights(apriori_hubs, K = K)
#' }

constructWeights <- function(d, K = NULL) {
  W <- vector(mode = "list", length = K)
  
  if (is.list(d) && length(d) > 1) {
    K <- length(d)
    names(W) <- names(d)
  } else {
    if (is.null(K)) {
      stop("Please provide the number of classes K.")
    } else {
      d_list <- vector(mode = "list", length = K)
      for (k in 1:K) {
        d_list[[k]] <- d
      }
      d <- d_list
    }
  }
  
  W <- vector(mode = "list", length = K)
  W <- lapply(d, function(x) 1 / sqrt(x %*% t(x)))
  W <- lapply(W, function(x) x / max(x))
  names(W) <- sapply(1:K, function(i) sprintf("W%i", i))
  
  if (sum(sapply(d, function(x) is.null(names(x)))) != 0) {
    warning("Some of degree vectors aren't named. This may lead to the wrong match between variables and their degrees later in jewel function. Please check.")
  } else {
    for (k in 1:K) {
      rownames(W[[k]]) <- colnames(W[[k]])
    }
  }
  
  return(W = W)
}
