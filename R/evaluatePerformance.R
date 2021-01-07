#' Evaluation of performance of graph estimation methods if true graph is known
#'
#' Function compares true and estimated simple graphs and calculates number of true positives (correctly estimated edges), true negatives (correctly estimated absence of edges), false positives (edges present in the estimator but not in the true graph) and false negatives (failure to identify an edge).
#'
#' @param G True graph
#' @param G_hat Estimated graph of the same size as \code{G}.
#'
#' @return vector of length 4 with TP, TN, FP, FN
#'
#' @export

evaluatePerformance <- function (G, G_hat) {

  G <- as.matrix(removeDiagonal(G))
  G_hat <- removeDiagonal(G_hat)

  comparison <- as.matrix(G == G_hat)

  TP <- sum(G[which(comparison)]) / 2
  TN <- sum(comparison) / 2 - TP
  FP <- sum(G_hat[which(!comparison)]) / 2
  FN <- sum(G[which(!comparison)]) / 2

  performance <- c(TP, TN, FP, FN);
  names(performance) <- c("TP", "TN", "FP", "FN")

  return(performance)
}
