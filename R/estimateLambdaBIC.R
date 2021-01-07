#' Estimation of optimal lambda parameter for jewel method with Bayesian information criterion
#'
#' Given a grid of lambda parameters, function evaluates Baeysian information criterion for each element.
#' Optimal lambda is chosen as a BIC's minimum. Warm start is implemented.
#'
#' @param X list of \code{K} data matrices of size \code{n_k} by \code{p} (\code{n_k} can be different for each class)
#' @param lambda vector of grid lambda parameters for which function evaluates BIC
#'
#' @return The following list is returned
#' \itemize{
#'   \item lambda_opt - a number, optimal value of regularization parameter according to BIC procedure
#'   \item CV_err - a vector of cross-validation erros for each element of input vector \code{lambda}.
#' }
#'
#' @export

estimateLambdaBIC <- function (X, lambda) {

  message("1/4 Starting iterations over lambda. Iteration number...")

  BIC <- rep(NA, length(lambda));

  for (l in 1:length(lambda)) {

    message(paste0("lambda ", l))

    if (l != 1) {
      GGM_result <- jewel(X, lambda[l],
                          Theta = Theta_warm_up,
                          reltol = 0.0001,
                          verbose = FALSE);
    } else GGM_result <- jewel(X, lambda[l],
                                reltol = 0.0001,
                                verbose = FALSE);

    BIC[l] <- GGM_result$BIC;
    Theta_warm_up <- GGM_result$Theta;

  }

  message("2/4 Completed. Estimating the optimal lambda...")

  #if minimum is achieved for two lambdas, we will get first one
  lambda_opt_BIC <- lambda[which.min(BIC)]

  message(paste0("3/4 BIC optimal lambda is ", lambda_opt_BIC, ". Generating the plot..."))

  plot(lambda, BIC, type = "p", col = "blue");
  points(lambda_opt_BIC, BIC[which.min(BIC)], col = "red");

  message("4/4 Completed.")

  return(list(lambda_opt = lambda_opt_BIC,
         BIC = BIC))
}
