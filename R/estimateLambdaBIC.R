#' Estimate the optimal regularization parameter for jewel method with Bayesian information criterion
#'
#' Function evaluates Baeysian information criterion (BIC) for each element of the grid.
#' Optimal lambda is chosen as the one for which BIC's minimum is obtained. Warm start is implemented (matrices computed with the previous parameter are used as a starting point for the current element of the grid).
#'
#' @param X a list of \code{K} numeric data matrices of \code{n_k} samples and \code{p} variables (\code{n_k} can be different for each matrix).
#' @param lambda an optional numeric vector of parameters for which the function evaluates BIC. Users are encouraged to tailor the grid to their specific data. If NULL, the default value is uniform in log sequence from 0.01 to 0.5. 
#' @param tol an optional number, convergence threshold controlling the relative error between iterations. The smaller it is, more precise is BIC. The default value is 0.0001.
#' @param makePlot If makePlot = FALSE, plotting of BIC is disabled. The default value is TRUE.
#'
#' @importFrom graphics points lines
#' @importFrom utils timestamp
#'
#' @return The following list is returned
#' \itemize{
#'   \item \code{lambda_opt} - a number, optimal value of regularization parameter according to BIC procedure;
#'   \item \code{BIC} - a numeric vector of BICs for each element of the grid of \code{lambda}.
#' }
#'
#' @export

estimateLambdaBIC <- function (X, lambda = NULL, tol = 0.0001,
                               makePlot = TRUE) {
  
  if (is.null(lambda)) {
    lambda <- exp(seq(log(0.01), log(0.5), length.out = 50))
  }
  BIC <- rep(NA, length(lambda))
  
  message("1/3 Starting iterations over lambda. Iteration number...")
  
  for (l in 1:length(lambda)) {
    
    timestamp()
    message(paste0("Computing for lambda", l))
    
    if (l != 1) {
      jewel_result <- jewel(X, lambda[l],
                             Theta = Theta_warm_up,
                             tol,
                             verbose = FALSE);
    } else {
      jewel_result <- jewel(X, lambda[l],
                             tol,
                             verbose = FALSE);
    }
    
    BIC[l] <- jewel_result$BIC
    Theta_warm_up <- jewel_result$Theta_list
  }
  
  #if minimum is achieved for two lambdas, we will get first one
  lambda_opt_BIC <- lambda[which.min(BIC)]
  
  message(paste0("2/3 Completed. BIC optimal lambda is ", lambda_opt_BIC, ". Generating the plot..."))
  
  if (makePlot == TRUE) {
    plot(lambda, BIC, type = "p", col = "blue")
    points(lambda_opt_BIC, BIC[which.min(BIC)], col = "red")
    lines(lambda, BIC, type = "p", col = "blue")
  }
  
  message("3/3 Completed.")
  
  return(list(lambda_opt = lambda_opt_BIC,
              BIC_values = BIC))
}
