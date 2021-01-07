#' Generate a sequence of lambda parameters
#'
#' Function that generates a uniform in logarithmic space grid of tuning parameters. Lambda_max is max(\code{1 / (n - 1)} max(X^TX)),  lambda_min = lambda_max * eps.
#'
#' @param X list of \code{K} data matrices of size \code{n_k} by {p} (\code{n_k} can be different for each class).
#' @param n desired number of parameters. The default value is 50.
#' @param eps lambda_min = lambda_max * eps. The default value is 0.1
#' @param scale if TRUE and data is scaled, than resulting grid is independent of X and goes from eps to 1 uniformly in log-scale. The default value is TRUE.
#'
#' @return vector of regularization parameters of length \code{n}
#' @export

generateLambdaGrid <- function (X, n = 50, eps = 0.1, scale = TRUE) {

  if (scale == TRUE) X <- lapply(X, scale)
  lambda_max <- max(sapply(X, function(x) (1 / (dim(x)[1] - 1)) * max(t(as.matrix(x)) %*% as.matrix(x))));
  lambda_min <- eps * lambda_max;
  lambda <- exp(seq(log(lambda_min), log(lambda_max), length.out = n));

  return(lambda);

}
