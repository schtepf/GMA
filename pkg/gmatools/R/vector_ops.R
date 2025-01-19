## Operations on Euclidean vectors

##' Unit vector for k-th axis in n-dimensional Euclidean space
##' 
##' Returns a unit vector for the `k`-th axis in `n`-dimensional Euclidean space
##' in the form of a column vector (or row vector if specified by the user).
##' 
##' @param n dimensionality of the Euclidean space
##' @param k which axis vector to return
##' @param column,byrow alternative flags to select a column vector (default) or row vector
##' @returns A \code{n x 1} matrix (default) or \code{1 x n} matrix (with `byrow=TRUE`)
##' 
##' @examples
##' axis.vector(5, 2)
##' axis.vector(5, 2, byrow=TRUE)
##' @export
axis.vector <- function (n, k, column=!byrow, byrow=FALSE) {
  stopifnot(1 <= k, k <= n)
  v <- rep(0, n)
  v[k] <- 1
  if (column) matrix(v, n, 1) else matrix(v, 1, n)
}

##' Gram-Schmidt orthonormalisation of a list of vectors
##' 
##' This function applies Gram-Schmidt orthonormalisation to a list of vectors
##' (given as columns of matrix `M`) via QR decomposition. This algorithm ensures
##' that the span of the first \eqn{k} vectors remains the same after orthonormalisation.
##' 
##' @param M matrix of one or more numeric column vectors
##' @param fail if `fail=FALSE` is specified, the function will silently return `NULL`
##'        if the vectors don't have full rank (rather than aborting with an error).
##' @returns A matrix of the same shape as `M` with the orthonormalised vectors (or `NULL`).
##' @examples
##' tmp <- matrix(1, 5, 3)
##' tmp <- tmp * (row(tmp) <= col(tmp))
##' tmp
##' gma.orthogonalize(tmp)
##' @export
gma.orthogonalize <- function (M, fail=TRUE) {
  stopifnot(is.matrix(M), is.numeric(M))
  d <- ncol(M) # dimensionality of basis
  res <- qr(M)
  if (res$rank < d) {
    if (fail) stop("vectors to be orthogonalized do not have full rank (not a basis)")
    return(NULL)
  }
  Q <- qr.Q(res) # orthogonal basis
  sign.adj <- sign(colSums(M * Q))  # 1 if orthogonal vector has same orientation as corresponding input vector
  .scaleMargins(Q, cols=sign.adj)   # so orthogonal vectors are projections of original basis vectors
}

##' Check orthogonality of a column matrix
##' 
##' This function checks whether a given column matrix is orthogonal, i.e. whether
##' its column vectors are orthonormal.
##' 
##' @param M matrix of one or more numeric column vectors
##' @param fail if TRUE, abort with error message if the matrix is not orthogonal
##' @param tol tolerance of orthogonality test (default: \eqn{10^{-7}})
##' @returns `TRUE` or `FALSE` (unless the function aborts because of `fail=TRUE`).
##' @export
gma.is.orthogonal <- function (M, fail=FALSE, tol=1e-7) {
  stopifnot(is.matrix(M), is.numeric(M))
  d <- ncol(M)
  A <- crossprod(M) # M'M should be diagonal matrix
  D <- A - diag(d)  # matrix of differences
  ok <- all(abs(D) < tol) 
  if (fail && !ok) stop(sprintf("%d x %d matrix is not orthogonal (maximum deviaton = %g)", nrow(M), d, max(abs(D))))
  ok
}
