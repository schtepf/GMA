## Our own implementation of LDA and repeated-measures LDA

##' Repeated-measures Linear Discriminant Analysis
##' 
##' Our own implementation of linear discriminant analysis (LDA), including
##' a novel repeated-measures version in order to remove confounding factors
##' from within-group variance. For a complete description of the algorithm
##' and its mathematical details, see 
##' [Mathematics of GMA](https://github.com/schtepf/GMA/blob/main/doc/gma_maths.pdf).
##'
##' @param x row matrix of data vectors (as in [MASS::lda()])
##' @param grouping factor sepcifying group membership of each data point (as in [MASS::lda()])
##' @param cohorts optional factor specifying cohort membership for each data point. 
##'        If specified, a repeated-measures LDA is carried out that factors cohorts out of within-group variance.
##' @param tol minimum condition number required for covariance matrix.
##'        Unlike [MASS::lda()], we do not yet support the case of a singular covariance matrix.
##' @param balanced if TRUE, all groups make the same contribution to the between-group covariance matrix (regardless of group size)
##' @param R2 choose as many discriminants as needed to cover at least this percentage of relative group separation variance (default: 100%)
##' @param ... any further argument are ignored, but allowed for compatibiltiy with [MASS::lda()]
##' @returns An object of class `lda` which mimics the value returned by [MASS::lda()]. In particular, the object includes components
##' - `scaling`: column matrix of LDA axis vectors; multiply with data matrix to obtain LDA dimension scores
##' - `prior`: the prior probabilities of class membership used by the algorithm (depending on the `balanced` flag)
##' - `means`: the group means in the original space (which can be transformed into LDA scores with `scaling`)
##' - `svd`: singular values for the ratio of between- and within-group variance, from which `R2` values are derived
##' - `counts`: number of data points in each group
##' - `lev`: group labels (the levels of `grouping`)
##' @export
gma.lda <- function (x, grouping, cohorts=NULL, tol=1e-4, balanced=FALSE, R2=100, ...) {
  X <- x # notation and equations follow "Mathematics of GMA" throughout the implementation
  n <- nrow(X)
  d <- ncol(X)
  m <- colMeans(X) # overall mean
  
  if (length(grouping) != n) stop("length of <grouping> doesn't match rows of data matrix <x>")
  gi <- droplevels(as.factor(grouping)) # only nonempty groups
  lev <- levels(gi)
  g <- length(lev)
  if (!(g > 1)) stop("need data points from two different groups at least")
  nj <- as.vector(table(grouping))
  names(nj) <- lev
  
  repeated <- FALSE # carry out repeated-measures LDA
  if (!is.null(cohorts)) {
    if (length(cohorts) != n) stop("length of <cohorts> doesn't match rows of data matrix <x>")
    ci <- droplevels(as.factor(cohorts)) # only nonempty cohorts
    clev <- levels(ci)
    c <- length(clev)
    nk <- as.vector(table(cohorts))
    names(nk) <- clev
    repeated <- TRUE
  }
  
  if (repeated) {
    njk <- as.matrix(table(grouping, cohorts))
    if (any(njk == 0)) stop("all (group, cohort) combinations must occur in the data set")
    cells <- data.frame(
      group = rep(lev, c),
      cohort = rep(clev, each=g),
      njk = as.vector(njk)
    )
    cells$cell <- paste(cells$group, cells$cohort)
    gici <- factor(paste(grouping, cohorts), levels=cells$cell) # (gi, ci) in standard ordering of cells
    M <- rowsum(X, gici) / cells$njk # cell means for each (group, cohort) combination
    MG <- rowsum(X, gi) / nj         # group means to be returned in lda object
    
    if (balanced) {
      prior <- structure(rep(1 / g, g), names=lev)
      C <- rowsum(M, factor(cells$cohort, levels=clev)) / g # cohort means for uniform prior
    } else {
      prior <- nj / n
      C <- rowsum(X, ci) / nk          # cohort means
    } 
    
    MC <- C[cells$cohort, , drop=FALSE] # cohort mean corresponding to each cell mean
    
    XW <- X - M[gici, , drop=FALSE]
    XW.df <- n - c * g
    MB <- (M - MC) * sqrt(prior[cells$group] * nk[cells$cohort])
    MB.df <- c * (g - 1)
  }
  else {
    M <- rowsum(X, gi) / nj # group means orderd according to lev
    MG <- M                 # to be returned in lda object
    if (balanced) {
      prior <- structure(rep(1/g, g), names=lev)
      m <- colMeans(M) # overall mean for balanced class distribution
    } else prior <- nj / n
    
    XW <- X - M[gi, ] # X_W = X - X_M
    XW.df <- n - g
    MB <- scale(M, center=m, scale=FALSE) * sqrt(prior * n) # M_B = diag(n_j)^1/2 * (M - 1 m^T)
    MB.df <- g - 1
  }
  
  XW.svd <- svd(XW, nu=0)
  D12 <- XW.svd$d / sqrt(XW.df) # diagonal of D^1/2
  ok <- length(D12) == d && min(D12) >= tol * max(D12)
  if (!ok) stop("within-group covariance matrix does not have full rank -- not supported by this implementation")
  U <- XW.svd$v
  S <- t(U) / D12 # S = D^-1/2 * U^T
  
  MBp <- tcrossprod(MB, S) # M_B' = M_B * S^T
  MBp.svd <- svd(MBp, nu=0)
  E <- MBp.svd$d^2 / MB.df # diagonal of E
  r <- length(E) # = rank(B') <= g - 1
  if (r < 1) stop("group means are numerically identical and cannot be separated")
  r <- min(r, sum(E >= tol * E[1]))
  if (R2 < 100) r <- min(r, sum(cumsum(E) / sum(E) < R2 / 100) + 1)
  Vr <- MBp.svd$v[, 1:r, drop=FALSE]
  E <- E[1:r]
  
  A <- crossprod(Vr, S) # A = A' * S = Vr^T * S
  colnames(A) <- colnames(X)
  rownames(A) <- sprintf("LD%d", 1:r)
  
  structure(list( # return same object as lda()
    prior=prior, counts=nj, means=MG,
    scaling=t(A), lev=lev, svd=sqrt(E),
    N=n, call=match.call()
  ), class="lda")
}
