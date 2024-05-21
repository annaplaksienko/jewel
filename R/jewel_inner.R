#' Estimate Gaussian graphical models from multiple datasets
#'
#' This function estimates Gaussian graphical models (i.e. networks of conditional dependencies, direct connections between variables) given several datasets. 
#' We assume that datasets contain measurements of the same variables collected under different conditions (different equipment, locations, even sub-types of disease).
#'
#' @param X a list of \code{K} numeric data matrices of \code{n_k} samples and \code{p} variables (\code{n_k} can be different for each matrix).
#' @param lambda1 a number, first regularization parameter (of the common penalty).
#' @param lambda2 an optional number, second regularization parameter (of the class-specific penalty). If NULL, set to \code{lambda_2 = lambda_1 * 1.4}
#' @param Theta an optional list of \code{K} regression coefficient matrices of the size \code{p} by \code{p}. User-provided initialization can be used for warm-start procedures. If NULL, initialized as all zeros.
#' @param W an optional list of \code{K} weights matrices of the size \code{p} by \code{p}. User-provided initialization can be used when some vertices are believed to be hubs. If NULL, initialized as all ones.
#' @param tol an optional number, convergence threshold controlling the relative error between iterations. The default value is 0.01.
#' @param maxIter an optional number, maximum allowed number of iterations. The default value is 10 000.
#' @param verbose if verbose = FALSE, tracing information printing is disabled. The default value is TRUE.
#'
#' @importFrom SMUT eigenMapMatMult
#'
#' @return The following list is returned
#' \itemize{
#'   \item \code{CommonG} - an adjacency matrix of the common estimated graph (intersection of \code{K} estimated graphs).
#'   \item \code{G_list} - a list of \code{K} adjacency matrices for each estimated graph.
#'   \item \code{Theta} - a list of \code{K} estimated covariance matrices.
#'   \item \code{BIC} – a number, value of Bayesian information criterion for resulting graphs.
#' }
#' @keywords internal

jewel_inner <- function(X, lambda1, lambda2 = NULL, 
                  Theta = NULL, W = NULL,
                  tol = 0.01, maxIter = 10000,
                  verbose = TRUE) {
  
  #normalize the data
  X <- mapply(function(y) scale(y), X, SIMPLIFY = FALSE)
 
  #get the number of input matrices
  K <- length(X)
  #get the dimensions of each matrix
  n_k <- sapply(X, function(x) dim(x)[1])
  p <- dim(X[[1]])[2]
  
  #check that variables are the same across classes
  vars <- colnames(X[[1]])
  if (is.null(vars)) {
    warning("Colnames of some datasets are empty. Be aware that without variable names there is no check if their order (columns order) is the same between datasets. That can lead to wrong estimation.")
  } else {
    if (sum(sapply(2:K, function(k) sum(colnames(X[[k]]) %in% vars) == p)) != (K - 1) |
        sum(sapply(2:K, function(k) sum(vars %in% colnames(X[[k]])) == p)) != (K - 1)) {
      stop("Variables don't match across classes. Check colnames of each element of X.")
    } else {
      #if variables are the same, make sure they are in the same order
      for (k in 1:K) {
        X[[k]] <- X[[k]][, vars]
      }
    }
  }
  
  
  #assemble a long X matrix (faster than lists)
  Xl <- do.call(rbind, X)
  nindex <- rep(1:K, n_k)
  index <- c(1:K)
  
  #Theta - list of regression coefficient matrices
  #Active - matrix indicating current active variables in the common part
  #Active_K_long - matrices indicating current active variables in each whole graph
  #R_specific and R_common - residuals
  #zero diagonals are dropped
  if (is.null(Theta)) {
    Theta <- Xi <- Gamma <- matrix(0, nrow = (p - 1) * K, ncol = p);
    Active <- matrix(TRUE, nrow = (p - 1), ncol = p);
    Active_K_long <- matrix(TRUE, nrow = (p - 1) * K, ncol = p);
    R_specific <- R_common <- Xl;
    
    #auxilary
    Gamma_list <- rep(list(NA), K)
    Xi_list <- rep(list(NA), K)
  } else {
    
    r <- mapply(function(x, y) x - eigenMapMatMult(x, y),
                X, Theta, SIMPLIFY = FALSE)
    R_specific <- R_common <- do.call(rbind, r)
    remove(r)
     
    Active <- lapply(Theta, function(x) x != 0)
    Active <- Reduce('+', Active)
    Active <- (Active == K)
    
    Xi <- lapply(Theta, function(x) x * Active)
    Gamma <- mapply(function(x, y) x - y,
                    Theta, Xi, SIMPLIFY = FALSE)
    
    Active_K <- lapply(Gamma, function(x) removeDiagonal(x != 0))
    Active_K_long <- do.call(rbind, Active_K)
    remove(Active_K)

    Gamma <- lapply(Gamma, removeDiagonal)
    Gamma <- do.call(rbind, Gamma)
    
    Xi <- lapply(Xi, removeDiagonal)
    Xi <- do.call(rbind, Xi)
    
    Theta <- lapply(Theta, removeDiagonal)
    Theta <- do.call(rbind, Theta)
    
    #auxilary for residual update
    Gamma_list <- rep(list(NA), K)
    Xi_list <- rep(list(NA), K)
  }
  
  
  #if weights are not provided, set them to 1
  if (is.null(W)) {
    W <- matrix(1, nrow = (p - 1) * K, ncol = p)
    colnames(W) <- vars
  } else {
    if (sum(sapply(W, function(x) !is.null(colnames(x)))) != K) {
      warning("Colnames of some weight matrices are empty. Be aware that without variable names there is no check if variable order matches with the one in datasets provided. That can lead to wrong estimation.")
      W <- lapply(W, removeDiagonal)
      W <- do.call(rbind, W)
    } else if (sum(sapply(1:K, function(k) sum(colnames(W[[k]]) %in% vars) == p)) != K |
        sum(sapply(1:K, function(k) sum(vars %in% colnames(W[[k]])) == p)) != K) {
      stop("There are some colnames that do not match between provided X and W. Please check.")
    } else {
      #order weights as variables in X
      W <- lapply(W, function(x) x[vars, vars])
      W <- lapply(W, removeDiagonal)
      W <- do.call(rbind, W)
    }
  }
  
  #if second regularization parameter not provided, 
  #set it to lambda_2 = 1.4 * lambda_1
  if (is.null(lambda2)) {
    lambda2 <- lambda1 * 1.4
  }
  
  #to avoid dividing by zero add eps to the denominator
  eps <- 2.220446e-16
  
  numIter <- 1
  check_conv <- 10000
  
  if (verbose) message("1/3 Initialization completed. Starting iterations...")
  
  while (numIter <= maxIter && check_conv > tol) {
    
    numIter <- numIter + 1
    if (verbose) message("jewel: iteration number ", numIter - 1)
    
    Theta_old <- Theta;
    
    order1 <- sample(1:(p-1), size = (p-1))
    
    #COMMON part (fix Gamma, update XI)
    for (j in order1) {
      
      #choose under-diagonal elements for j-th column
      #get indices of those which are TRUE 
      Act <- (j-1) + which(Active[j:(p - 1), j], arr.ind = TRUE)
      if (length(Act) == 0) break
      jminus <- setdiff(1:p, j)
      
      for (a in 1:length(Act)) {
        
        aminus <- setdiff(1:p, Act[a] + 1)
        
        za <- c(NA, K)
        zb <- c(NA, K)
        
        za <- sapply(index, function (k)
          (1 / n_k[k]) * Xl[nindex == k, jminus[Act[a]]] %*% R_common[nindex == k, j] + 
            Xi[(p - 1) * (k - 1) + Act[a], j])
        zb <- sapply(index, function (k)
          (1 / n_k[k]) * Xl[nindex == k, aminus[j]] %*% R_common[nindex == k, Act[a] + 1] + 
            Xi[(p - 1) * (k - 1)  + j, Act[a] + 1])
        z <- c(za, zb);
        
        W_avr <- mean(sapply(index, function(k) W[(p - 1) * (k - 1)  + j, Act[a] + 1]))
        thrld <- 1 - lambda1 * sqrt(2 * K) * W_avr / (sqrt(sum(z^2)) + eps)
        
        if (thrld <= 0) {
          z <- z * 0;
          Active[j, Act[a] + 1] <-  FALSE;
          Active[Act[a], j] <- FALSE;
        } else {
          z <- z * thrld;
        }
        
        za <- z[1:K];
        zb <- z[(K + 1) : (2 * K)];
        
        R_common[ , j] <- unlist(lapply(index, function (k)
          R_common[nindex == k, j] - Xl[nindex == k, jminus[Act[a]]] *
            (za[k] - Xi[(p - 1) * (k - 1) + Act[a], j])))
        
        R_common[, Act[a] + 1] <- unlist(lapply(index, function (k)
          R_common[nindex == k, Act[a] + 1] - Xl[nindex == k, aminus[j]] *
            (zb[k] - Xi[(p - 1) * (k - 1) + j, Act[a] + 1])))
        
        for (k in index) {
          Xi[(p - 1) * (k - 1) + Act[a], j] <- za[k]
          Xi[(p - 1) * (k - 1) + j, Act[a] + 1] <- zb[k]
        }
      }
    }
    
    #update R_specific
    {
      for (i in 1:K) {
        Gamma_list[[i]] <- addZeroDiagonal(Gamma[((p - 1) * (i - 1) + 1) : ((p - 1) * i),]);
        Xi_list[[i]] <- addZeroDiagonal(Xi[((p - 1) * (i - 1) + 1) : ((p - 1) * i),]);
      }
      
      update <- mapply(function(x, y, z) x - eigenMapMatMult(x, y) - eigenMapMatMult(x, z),
                       X, Gamma_list, Xi_list, SIMPLIFY = FALSE)
      
      R_specific <-  do.call(rbind, update)
    }
    
    order2 <- sample(1:(p-1), size = (p-1))
    
    #SPECIFIC part (fix Xi, update GAMMA)
    for (k in 1:K) {
      for (j in order2) {
        
        #choose under-diagonal elements for j-th column
        #get indices of those which are TRUE 
        Act <- (p - 1) * (k - 1) + (j - 1) + 
          which(Active_K_long[((p - 1) * (k - 1) + j): ((p - 1) * k), j], 
                arr.ind = TRUE);
        if (length(Act) == 0) break;
        
        #for all active variables in the current column
        for (a in 1:length(Act)) {
          
          i <- Act[a] + k - p * (k - 1)
          
          za <- (1 / n_k[k]) * Xl[nindex == k, i] %*% R_specific[nindex == k, j] + 
            Gamma[Act[a], j]
          zb <- (1 / n_k[k]) * Xl[nindex == k, j] %*% R_specific[nindex == k, i] + 
            Gamma[(p - 1) * (k - 1)  + j, i]
          z <- c(za, zb);
          
          thrld <- 1 - lambda2 * sqrt(2) * W[(p - 1) * (k - 1)  + j, i] / (sqrt(sum(z^2)) + eps)
          
          if (thrld <= 0) {
            z <- z * 0;
            Active_K_long[Act[a], j] <- FALSE
            Active_K_long[(p - 1) * (k - 1) + j, i] <- FALSE
          } else {
            z <- z * thrld;
          }
          
          za <- z[1]
          zb <- z[2]
          
          R_specific[nindex == k, j] <- R_specific[nindex == k, j] - Xl[nindex == k, i] *
            (za - Gamma[Act[a], j])
          
          R_specific[nindex == k, i] <- R_specific[nindex == k, i] - Xl[nindex == k, j] *
            (zb - Gamma[(p - 1) * (k - 1) + j, i])
          
          Gamma[Act[a], j] <- za
          Gamma[(p - 1) * (k - 1) + j, i] <- zb
          
        }
      }
    }
    
    #update R_common
    {
      for (i in 1:K) {
        Gamma_list[[i]] <- addZeroDiagonal(Gamma[((p - 1) * (i - 1) + 1) : ((p - 1) * i),]);
        Xi_list[[i]] <- addZeroDiagonal(Xi[((p - 1) * (i - 1) + 1) : ((p - 1) * i),]);
      }
      
      update <- mapply(function(x, y, z) x - eigenMapMatMult(x, y) - eigenMapMatMult(x, z),
                       X, Gamma_list, Xi_list, SIMPLIFY = FALSE)
      
      R_common <-  do.call(rbind, update)
    }
    
    #check convergence
    Theta <- Xi + Gamma
    check_conv <- sum(abs(Theta - Theta_old)) / (sum(abs(Theta_old)) + eps)
  }
  
  if (verbose) message("jewel: total number of iterations is ", numIter - 1, 
                  " and the error is ", round(check_conv, digits = 5))
  
  if (verbose) message("2/3 Iterations completed. Assembling the output...")
  
  #construct a list of matrices from long matrices
  Theta_list <- rep(list(NA), K)
  A_list <- rep(list(NA), K)
  names(A_list) <- names(X)
  
  for (i in 1:K) {
    Theta_list[[i]] <- addZeroDiagonal(Theta[((p - 1) * (i - 1) + 1) : ((p - 1) * i), ])
    A_list[[i]] <- (Theta_list[[i]] != 0)
    colnames(A_list[[i]]) <- rownames(A_list[[i]]) <- colnames(X[[i]])
  }
  
  #common graph is the intersection of K estimated graphs
  A <- Reduce('+', A_list)
  A <- (A == K)
  
  if (verbose) message("3/3 Completed.")
  
  return(list(G_list = A_list,
              CommonG = A,
              Theta = Theta_list));
}
