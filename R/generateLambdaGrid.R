#' Generation of the sequence of regularization parameters
#'
#' Function generates a uniform in logarithmic space grid of regularization parameters. \code{Lambda_max} is \code{max({1 / (n - 1)} max(X^TX))},  \code{lambda_min = lambda_max * eps}.
#'
#' @param X list of \code{K} numeric data matrices of size \code{n_k} by {p} (\code{n_k} can be different for each matrix)
#' @param n desired number of parameters. The default value is 50.
#' @param eps lambda_min = lambda_max * eps. The default value is 0.1
#' @param scale if TRUE and hence data is scaled, then resulting grid is independent of \code{X} and goes from \code{eps} to 1 uniformly in log-scale. The default value is TRUE.
#'
#' @return \code{lambda} - vector of regularization parameters of length \code{n}
#' @export

generateLambdaGrid <- function (X, n = 50, eps = 0.1, scale = TRUE) {

  if (scale == TRUE) X <- lapply(X, scale)
  lambda_max <- max(sapply(X, function(x) (1 / (dim(x)[1] - 1)) * max(t(as.matrix(x)) %*% as.matrix(x))));
  lambda_min <- eps * lambda_max;
  lambda <- exp(seq(log(lambda_min), log(lambda_max), length.out = n));

  return(lambda);

}
