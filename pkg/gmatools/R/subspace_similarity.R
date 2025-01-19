## Various measures of subspace similarity

##' Compute similarity between two subspaces A and B
##' 
##' Computes various measures of the similarity of two subspaces `A` and `B`,
##' which are meaningful even if the dimensionalities of the two subspaces differ.
##' 
##' @param A column matrix of basis vectors of subspace A
##' @param B column matrix of basis vectors of subspace B
##' @param orthogonalize can be set to FALSE if `A` and `B` are already orthnormal bases
##' @param tol tolerance for testing orthogonality in this case (default: \eqn{10^{-7}})
##' @param basis matched basis vectors for the two subspaces are attached as attributes if `basis=TRUE`
##' @param method similarity measure: number of shared dimensions (`dim`), preserved variance (`R2`), 
##'        or vector of singular values for matched basis vectors (`sigma`)
##' @returns Either a numeric similarity score or a vector of singular values. If `basis=TRUE`, 
##'          attributes `basisA` and `basisB` provide column matrices of matched basis vectors
##' 
##' @details
##' Subspace similarity is based on an SVD of the crossproduct of the two sets of basis vectors,
##' which determines orthogonal transformations of both bases such that corresponding basis 
##' vectors are optimally matched. The singular values are cosine similarities between pairs of
##' matched basis vectors. Depending on `method`, subspace similarity is quantified as:
##' - `dim`: sum of the singular values, which can be understood as the fractional number of
##'          shared dimensions between the subspaces. Its range is \eqn{[0, \min \{\dim A, \dim B\}]}.
##' - `R2`: average proportion of preserved variance if random vectors from B are projected into A.
##'         Its range is \eqn{[0, 1]}.
##' - `sigma`: vector of individual singular values for matched basis vectors
##' If both subspaces have the same dimensionality, then all measures are symmetric. Otherwise,
##' similarity is computed for a projection from B into A.
##' 
##' Further information and some mathematical details are sketched in [SIGIL Unit 7](https://sigil.r-forge.r-project.org/#unit07).
##' @export
gma.similarity <- function (A, B, method=c("dimensions", "R2", "sigma"), orthogonalize=TRUE, tol=1e-7, basis=FALSE) {
  method <- match.arg(method)
  .ensure.matrix(A, "argument A= of gma.similarity()")
  .ensure.matrix(B, "argument B= of gma.similarity()")
  nA <- ncol(A)
  nB <- ncol(B)
  if (nA == 0L || nB == 0L) {
    return( if (method == "sigma") numeric(0) else 0 )
  }
  if (orthogonalize) {
    A <- gma.orthogonalize(A, fail=TRUE)
    B <- gma.orthogonalize(B, fail=TRUE)
  } else {
    if (!gma.is.orthogonal(A, tol=tol)) stop("basis of subspace A is not orthonormal")
    if (!gma.is.orthogonal(A, tol=tol)) stop("basis of subspace B is not orthonormal")
  }
  res <- svd(crossprod(A, B), nu=nA, nv=nB) # A'B = UDV'
  U <- res$u
  sigma <- res$d
  V <- res$v
  retval <- switch(
    method,
    dimensions = sum(sigma),
    R2 = mean(sigma^2),
    sigma = sigma,
    stop("internal error"))
  if (basis) {
    attr(retval, "basisA") <- A %*% U
    attr(retval, "basisB") <- B %*% V
  }
  retval
}
