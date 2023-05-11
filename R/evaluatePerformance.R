#' Evaluation of graph estimation method's performance if the true graph is known.
#'
#' Function compares adjacency matrices of the true and estimated simple graphs and calculates the number of true positives (correctly estimated edges), true negatives (correctly estimated absence of edges), false positives (edges present in the estimator but not in the true graph) and false negatives (failure to identify an edge).
#'
#' @param G true graph's adjacency matrix.
#' @param G_hat estimated graph's adjacency matrix. Must have the same dimensions as \code{G}.
#'
#' @return performance - a numeric vector of length 4 with TP, TN, FP, FN.
#'
#' @export
#' 
#' @examples
#' {
#' K <- 3
#' p <- 50
#' n <- 20
#' data <- generateData_rewire(K = K, p = p, n = n, ncores = 1, verbose = FALSE)
#' G_common_true <- data$CommonGraph
#' X <- data$Data
#' res <- jewel(X, lambda1 = 0.25)
#' G_common_est <- res$CommonG
#' evaluatePerformance(G = G_common_true, G_hat = G_common_est)
#' }

evaluatePerformance <- function (G, G_hat) {

  G <- removeDiagonal(as.matrix(G))
  G_hat <- removeDiagonal(G_hat)

  comparison <- (G == G_hat)

  TP <- sum(G[which(comparison)]) / 2
  TN <- sum(comparison) / 2 - TP
  FP <- sum(G_hat[which(!comparison)]) / 2
  FN <- sum(G[which(!comparison)]) / 2

  performance <- c(TP, TN, FP, FN);
  names(performance) <- c("TP", "TN", "FP", "FN")

  return(performance)
}
