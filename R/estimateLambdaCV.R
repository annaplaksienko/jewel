#' Estimation of the optimal regularization parameter for jewel method based on cross-validation
#'
#' Given a grid of regulariation parameters, function performs cross-validation and 
#' estimates the optimal parameter as the one for which the cross-validation error's minimum is obtained. 
#' Parallelization over folds and warm start are implemented.
#'
#' @param X list of \code{K} numeric data matrices of size \code{n_k} by \code{p} (\code{n_k} can be different for each matrix)
#' @param lambda vector of parameters over which cross-validation is performed
#' @param k_folds number of folds in which data is divided. The default value is 5.
#' @param verbose If verbose = FALSE, tracing information printing is disabled. The default value is TRUE.
#' @param makePlot If makePlot = FALSE, plotting of CV error is disabled. The default value is TRUE.
#'
#' @importFrom SMUT eigenMapMatMult
#' @import parallel
#'
#' @return The following list is returned
#' \itemize{
#'   \item \code{lambda_opt} - a number, optimal value of regularization parameter according to cross-validation procedure;
#'   \item \code{CV_err} - a vector of cross-validation errors for each element of input vector \code{lambda}
#' }
#'
#' @export

estimateLambdaCV <- function (X, lambda, k_folds = 5, verbose = TRUE, makePlot = TRUE) {

  #parallelization of cross validation over folds with warm start
  cross_validation <- function(fold) {

    tol <- 0.0001

    err_lambda <- rep(NA, length(lambda))

    X_F <- mapply(function(x, y) x[which(y == fold), ], X, folds, SIMPLIFY  = FALSE)
    X_minus_F <- mapply(function(x, y) x[-which(y == fold), ], X, folds, SIMPLIFY  = FALSE)

    for (l in 1:length(lambda)) {

      if (l > 1) {
        Theta_minus_F <- jewel(X = X_minus_F, lambda[l], Theta = Theta_minus_F,
                                          tol = tol, verbose = FALSE)$Theta
      } else Theta_minus_F <- jewel(X = X_minus_F, lambda[l],
                                               tol = tol, verbose = FALSE)$Theta


      err_lambda[l] <- sum(sapply(index, function(c)
        norm(X_F[[c]] - eigenMapMatMult(X_F[[c]], Theta_minus_F[[c]]), type = "F")^2 / n_F[c]))

      if (err_lambda[l] == 0) break;
    }

    names(err_lambda) <- sapply(1:length(lambda), function(i) sprintf("lambda %i", i))

    return(err_lambda)
  }

  #get number of classes
  K <- length(X);
  #get dimensions of each matrix
  n_k <- sapply(X, function(x) dim(x)[1])
  p <- dim(X[[1]])[2]
  index <- c(1:K)

  if (verbose) message("1/4 Constructing folds...")

  #size of each fold for different n_k
  n_F = n_k / k_folds;

  #for each category we create it's own vector of folds
  folds <- rep(list(NA), K)
  folds <- lapply(index, function (c) sample(rep(1:k_folds, each = n_F[c])))

  if (verbose) message("2/4 Completed. Starting cross-validation...")

  #start the cluster
  ncores <-  detectCores(logical = FALSE) - 1;
  cl <- makeCluster(ncores, type = "SOCK")
  clusterEvalQ(cl, library("jewel"))
  clusterEvalQ(cl, library("SMUT"))

  clusterExport(cl, varlist = c("X", "lambda", "folds", "index", "n_F"), envir = environment())
  err <- clusterApply(cl, c(1:k_folds), cross_validation)
  names(err) <- sapply(1:k_folds, function(i) sprintf("Fold %i", i))

  stopCluster(cl)

  err <- matrix(unlist(err), nrow = k_folds, ncol = length(lambda), byrow = TRUE)

  #if minimum is achieved for two lambdas, we will get first one
  CV_err <- colSums(err) / 5
  lambda_opt_CV <- lambda[which.min(CV_err)]

  if (verbose) message(paste0("3/4 Completed. CV optimal lambda is ", lambda_opt_CV, ". Generating the plot..."))

  if (makePlot == TRUE) {
    plot(lambda, CV_err, type = "p", col = "green")
    points(lambda_opt_CV, CV_err[which.min(CV_err)], col = "red");
  }


  if (verbose) message(paste0("4/4 Completed."))

  return(list(lambda_opt = lambda_opt_CV, CV_err = CV_err))

}


