#' Removing diagonal from a matrix
#'
#' Function removes a diagonal from a square \code{p} by \code{p} matrix
#' and returns \code{(p-1)} by \code{p} matrix
#'
#' @param M a matrix in which you need to remove a diagonal
#' @keywords internal

removeDiagonal <- function(M) {
  matrix(M[lower.tri(M, diag = F) | upper.tri(M, diag = F)], 
         nrow = nrow(M) - 1, ncol = ncol(M))
}

#' Adding zero diagonal to a matrix
#'
#' Function adds zero diagonal to \code{(p-1)} by \code{p} matrix and returns \code{p} by \code{p} matrix
#'
#' @param M a matrix to which you need to add zero diagonal
#' @keywords internal

addZeroDiagonal <- function(M) {

  p <- max(dim(M))
  M_new <- matrix(NA, p, p)
  diag(M_new) <- 0;
  #replacing everything except the diagonal in the auxiliary matrix with Theta
  M_new[which(is.na(M_new), arr.ind=TRUE)] <- as.vector(M);

  return(M_new)
}
