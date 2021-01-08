#' Evaluation of graph estimation methods perfomance if true graph is known
#'
#' Function compares adjacency matrices of true and estimated simple graphs and calculates the number of true positives (correctly estimated edges), true negatives (correctly estimated absence of edges), false positives (edges present in the estimator but not in the true graph) and false negatives (failure to identify an edge).
#'
#' @param G True graph's adjacency matrix
#' @param G_hat Estimated graph's adjacenct matrix. Must have the same dimensions as \code{G}.
#'
#' @return performance - vector of length 4 with TP, TN, FP, FN
#'
#' @export

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
