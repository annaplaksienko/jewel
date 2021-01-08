#' Joint node-wise estimation of Gaussian graphical model from multiple datasets
#'
#' Implementation of the jewel method for estimation of the graph of conditional dependencies between the variables given multiple datasets, 
#' i.e. when observations of variables are collected under different conditions.
#'
#' @param X List of \code{K} numeric data matrices of size \code{n_k} by \code{p} (\code{n_k} can be different for each matrix).
#' @param lambda Regularization parameter which controls the sparsity of the resulting graph - bigger it is, less edges one gets.
#' @param Theta List of \code{K} starting regression coefficient matrices of size \code{p} by \code{p}. If not provided, initialized as all zeros.
#' @param reltol Convergence threshold controlling the relative error between iterations. The default value is 0.01.
#' @param maxIter Maximum allowed number of iterations. The default value is 10 000.
#' @param verbose If verbose = FALSE, tracing information printing is disabled. The default value is TRUE.
#'
#' @importFrom SMUT eigenMapMatMult
#'
#' @return The following list is returned
#' \itemize{
#'   \item \code{EstAdjMat} - adjacency matrix of the estimated graph
#'   \item \code{Theta} - list of \code{K} estimated covariance matrices of size \code{p} by \code{p}
#'   \item \code{residual} - list of \code{K} matrices of residuals of size \code{n_k} by \code{p} (\code{n_k} can be different for each matrix)
#'   \item \code{BIC} - value of Bayesian information criterion
#' }
#'
#' @export


jewel <- function (X, lambda, Theta = NULL,
                              reltol = 0.01, maxIter = 10000,
                              verbose = TRUE) {

  #get number of input matrices
  K <- length(X)
  #get dimensions of each matrix
  n_k <- sapply(X, function(x) dim(x)[1])
  p <- dim(X[[1]])[2]

  #normalize the data
  X <- mapply(function(y) scale(y), X, SIMPLIFY = FALSE);

  #assemble a long X matrix
  Xl <- do.call(rbind, X)
  nindex <- rep(1:K, n_k)
  pindex <- rep(1:K, rep(p, K))
  index <- c(1:K)

  #Theta - list of regression coefficient matrices
  #active - matrix indicating the active variables
  #r - list of residuals for each class
  if (is.null(Theta) == 1) {
    Theta <- matrix(0, nrow = p - 1, ncol = p);
    Theta <- do.call(rbind, rep(list(Theta), K))
    active <- matrix(TRUE, nrow = p - 1, ncol = p);
    r <- Xl;
  } else {
    r <- mapply(function(x, y) x - eigenMapMatMult(x, y),
                X, Theta, SIMPLIFY = FALSE);
    r <- do.call(rbind, r);

    Theta <- lapply(Theta, removeDiagonal);
    active <- (Theta[[which.max(sapply(Theta, function(x) sum(x != 0)))]] != 0);
    Theta <- do.call(rbind, Theta)
  }

  #to avoid dividing by zero
  eps = 2.220446e-16

  numIter <- 1;
  check_conv <- 10000;

  if (verbose) message("1/3 Initialization completed. Starting iterations. Iteration number...")

  while (numIter <= maxIter && check_conv > reltol) {

    numIter <- numIter + 1;
    if (verbose) message(paste0(numIter - 1));

    Theta_old <- Theta;

    for (j in 1:(p-1)) {

      A = (j-1) + which(active[j:nrow(active),j], arr.ind = TRUE);
      if (length(A) == 0) break;
      jminus <- setdiff(1:p, j);

      for (l in 1:length(A)) {

        lminus <- setdiff(1:p, A[l] + 1);

        za <- c(NA, K)
        zb <- c(NA, K)

        za <- sapply(index, function (c)
          (1 / n_k[c]) * Xl[nindex == c, jminus[A[l]]] %*% r[nindex == c, j]
          + Theta[(p - 1) * (c - 1) + A[l], j])
        zb <- sapply(index, function (c)
          (1 / n_k[c]) * Xl[nindex == c, lminus[j]] %*% r[nindex == c, A[l] + 1]
          + Theta[(p - 1) * (c - 1)  + j, A[l] + 1])
        z <- c(za, zb);

        thrld <- 1 - lambda * sqrt(2 * K) / (sqrt(sum(z^2)) + eps)

        if (thrld <= 0) {
          z <- z * 0;
          active[j, A[l] + 1] <- FALSE;
          active[A[l], j] <- FALSE;
        } else {
          z <- z * thrld;
        }

        za <- z[1:K];
        zb <- z[(K + 1) : (2 * K)];


        r[ , j] <- unlist(lapply(index, function (c)
          r[nindex == c, j] - Xl[nindex == c, jminus[A[l]]] *
            (za[c] - Theta[(p - 1) * (c - 1) + A[l], j])))

        r[, A[l] + 1] <- unlist(lapply(index, function (c)
          r[nindex == c, A[l] + 1] - Xl[nindex == c, lminus[j]] *
            (zb[c] - Theta[(p - 1) * (c - 1) + j, A[l] + 1])))

        for (c in index) {
          Theta[(p - 1) * (c - 1) + A[l], j] <- za[c]
          Theta[(p - 1) * (c - 1) + j, A[l] + 1] <- zb[c]
        }

      }

    }

    check_conv <- sum(abs(Theta - Theta_old)) / (sum(abs(Theta_old)) + eps)
  }

  message((paste0("jewel: Total number of iterations is ", numIter - 1, " and error is ", check_conv)))

  if (verbose) message("2/3 Iterations completed. Evaluating BIC and assembling the output...")

  #evaluate BIC
  BIC <- sum(sapply(index, function(c) n_k[c] * sum(apply(r[nindex == c, ], 2,
                                                          function (y) log( sum(y^2) ) ) ) ) ) +
          sum(sapply(n_k, log)) * (sum(active) / 2)

  active <- addZeroDiagonal(active)

  Theta_list <- rep(list(NA), K)

  for (i in 1:K) {
    Theta_list[[i]] <- addZeroDiagonal(Theta[((p - 1) * (i - 1) + 1) : ((p - 1) * i),]);
  }

  if (verbose) message("3/3 Completed.")

  return(list(EstAdjMat = active,
              Theta = Theta_list,
              residual = r,
              BIC = BIC));
}
