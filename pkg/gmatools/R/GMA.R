## GMA: a class for working with orthogonal subspace projections

##' R6 Class for Geometric Multivariate Analysis
##' 
##' @description
##' **TODO**
##' 
##' @returns Unless otherwise specified, methods invisibly return the object so they can be chained.
##' @export
##' @importFrom R6 R6Class
##' @importFrom stats prcomp
GMA <- R6Class("GMA", list(
  #' @field P orthonormal basis of focus space (column matrix)
  P = NULL,
  #' @field Q orthonormal basis of complement space (column matrix)
  Q = NULL,
  #' @field axes basis spanning focus space, as specified by user (column matrix) 
  axes = NULL,
  #' @field data data matrix of feature vectors in original coordinates (row matrix)
  data = NULL,
  #' @field n dimensionality of the feature space
  n = 0L,
  #' @field rotated TRUE if rotation has been applied, which invalidates `axes`
  rotated = FALSE,
  
  #' @description
  #' Create a new GMA object. **TODO**
  #' @param M numeric matrix of feature vectors (as rows) representing the main data set to be analysed
  #' @param n.dim dimensionality of the feature space, to be specified if `M` is not provided
  #' @return A new `GMA` object.
  initialize = function (M=NULL, n.dim=0L) {
    if (!is.null(M)) {
      .ensure.matrix(M, "argument M= of GMA$new()")
      if (n.dim != 0L && n.dim != ncol(M)) stop("n.dim inconsistent with data matrix M")
      n.dim <- ncol(M)
    }
    else {
      n.dim <- as.integer(n.dim)
      M <- diag(n.dim)
    }
    if (n.dim == 0) stop("can't create GMA() object for 0-dimensional space")
    self$n <- n.dim
    self$data <- M
    self$P <- self$axes <- matrix(0, nrow=n.dim, ncol=0L)
    self$Q <- prcomp(M, center=TRUE, scale=FALSE)$rotation
    self$rotated <- FALSE
  },
  
  #' @description
  #' Print a brief description of a GMA object.
  print = function () {
    k <- nrow(self$data)
    d1 <- ncol(self$P)
    cat(sprintf("GMA object representing projection of %d x %d data matrix into %d-dimensional subspace", k, self$n, d1))
    if (isTRUE(self$rotated)) cat(", with rotation")
    cat("\n")
    invisible(self)
  },
  
  #' @description
  #' Retrieve basis of focus space and/or complement
  #' @param space select desired basis vectors to be returned
  #'   - `focus`: orthonormal basis of focus space (default)
  #'   - `complement` orthonormal basis of complement space
  #'   - `both`: orthonormal basis of focus space + complement
  #'   - `axes`: user-specified basis of focus space (need not be orthonormal)
  #' @param dim = return only these dimensions of the selected basis
  #' @return A column matrix of basis vectors.
  basis = function (space=c("focus", "both", "complement", "axes"), dim=NULL) {
    space <- match.arg(space)
    res <- switch(space,
                  "both" = cbind(self$P, self$Q),
                  "focus" = self$P,
                  "complement" = self$Q,
                  "axes" = self$axes)
    if (!is.null(dim)) res[, dim, drop=FALSE] else res
  },

  #' @description
  #' Project data points into subspace and/or complement
  #' @param space select subspace and type of projection
  #'   - `focus`: orthogonal projection into focus space (default)
  #'   - `complement` orthogonal projection into complement space
  #'   - `both`: transform into orthonormal basis of focus space + complement
  #'   - `axes`: focus space coordinates in user-specified basis
  #' @param dim = return only these coordinates of the selected projection
  #' @param M = data matrix of row vectors to project (defaults to internal data matrix)
  #' @param byrow set `byrow=FALSE` if `M` is a matrix of column vectors
  #' @return A matrix of row vectors with the coordinates of data points after projection (or column vectors for `byrow=FALSE`).
  projection = function (space=c("focus", "both", "complement", "axes"), dim=NULL, M=NULL, byrow=TRUE) {
    space <- match.arg(space)
    if (space == "axes") stop("space='axes' is not supported yet")
    if (is.null(M)) {
      M <- self$data
      if (!byrow) M <- t(M)
    } else {
      n. <- if (byrow) ncol(M) else nrow(M)
      if (n. != self$n) stop("dimensionality of specified data matrix M doesn't match GMA space")
    }
    P. <- self$basis(space, dim=dim) # orthonormal basis for the projection
    if (byrow) M %*% P. else crossprod(P., M)
  },
  
  #' @description
  #' Compute proportion of variance (\eqn{R^2}) captured by projection into focus space.
  #' The return value is a vector of \eqn{R^2} contributions (%) for each selected dimension.
  #' Note that additional complement dimensions can be included with the `dim` argument.
  #' @param M row matrix of data points to be projected (defaults to original data matrix)
  #' @param dim = vector of dimensions for which \eqn{R^2} is computed (defaults to complete focus space)
  #' @return A numeric vector specifying partial \eqn{R^2} for each selected dimensions as a percentage.
  R2 = function (M=NULL, dim=NULL) {
    if (is.null(M)) M <- self$data
    if (is.null(dim)) dim <- seq_len(ncol(self$P))
    M <- scale(M, center=TRUE, scale=FALSE) # center data set
    full.ss <- sum(M^2) # denominator (n-1) of variance is constant and can be ignored
    Mproj <- self$projection(space="both", M=M)
    Mproj <- Mproj[, dim, drop=FALSE] # select relevant dimensions
    if (ncol(Mproj) == 0L) return(numeric(0))
    partial.ss <- colSums(Mproj^2) # Mproj is also centered
    100 * partial.ss / full.ss
  },
  
  #' @description
  #' Add one or more basis dimensions to the focus space. The user-specified basis vectors must be 
  #' a linearly independent extension of the existing focus space. They are automatically transformed
  #' into an orthonormal basis extension of the focus space.
  #' @param basis column matrix of basis vectors to be added (a single basis vector can be given as plain vector)
  #' @param space orthonormal coordinate system in which basis vectors are specified
  #'   - `orig`: original coordinates (feature space, the default)
  #'   - `both`: coordinates in orthonormal basis (of focus space + complement)
  #'   - `complement`: coordinates within complement space
  #' @param names optional names for the additional dimensions
  #' @param normalize normalise specified basis vectors to unit length if `normalize=TRUE`
  #' @param tol tolerance for internal orthogonality checks
  add = function (basis, space=c("orig", "both", "complement"), names=colnames(basis), normalize=FALSE, tol=1e-7) {
    space <- match.arg(space)
    basis <- .ensure.matrix(basis, "argument basis= of $add() method", vector.ok=TRUE)
    d1 <- ncol(self$P)  # dimensionality of current focus space
    d2 <- ncol(basis)   # number of basis vectors to be added
    nb <- nrow(basis)   # dimensionality of basis vectors
    n <- self$n
    
    if (space == "complement") {
      if (nb != n - d1) stop(sprintf("basis vectors of dimensionality %d don't match %d-dimensional complement space", nb, n-d1))
      basis <- self$Q %*% basis  # transform basis vectors into original space
    } else if (space == "both") {
      if (nb != n) stop(sprintf("basis vectors of dimensionality %d don't match %d-dimensional orthonormal space", nb, n))
      basis <- cbind(self$P, self$Q) %*% basis # transform into original space
    } else {
      if (nb != n) stop(sprintf("basis vectors of dimensionality %d don't match %d-dimensional original space", nb, n))
    }
    if (normalize) basis <- .scaleMargins(basis, cols=1 / .euclideanNorms(basis, byrow=FALSE))
    
    Pnew <- gma.orthogonalize(cbind(self$P, basis), fail=FALSE) # extend projection subspace
    if (is.null(Pnew)) stop("basis vectors do not form a linearly independent extension of the current subspace")
    if (d1 > 0) {
      if (max(abs(Pnew[, seq_len(d1), drop=FALSE] - self$P)) > tol) stop("internal error -- QR decomposition does not preserve existing orthonormal basis")
    }
    if (!is.null(names)) {
      if (length(names) != d2) stop("length of names= doesn't match number of basis vectors")
      colnames(basis) <- names
    }
    self$axes <- cbind(self$axes, basis) # extend list of original basis vectors
    colnames(Pnew) <- colnames(self$axes)
    self$P <- Pnew # orthonormal basis of new subspace
    
    M <- self$data
    Mcomp <- M - M %*% tcrossprod(Pnew) # M * (I - P P') = projection of row matrix M into new complement space
    Qnew <- prcomp(Mcomp, center=TRUE, scale=FALSE)$rotation # PCA orthogonalization of complement space
    Qnew <- Qnew[, seq_len(n - d1 - d2), drop=FALSE]
    if (!gma.is.orthogonal(cbind(Pnew, Qnew), tol=tol)) stop("internal error -- basis of focus space + complement is not orthonormal")
    self$Q <- Qnew
    
    rownames(self$P) <- rownames(self$Q) <- rownames(self$axes) <- colnames(self$data)
    invisible(self)
  },
  
  #' @description
  #' Apply operation in subspace coordinates and re-transform resulting vectors into original feature space.
  #' Only orthogonal transformations are allowed.
  #' @param M data matrix of row or column vectors to be operated on
  #' @param FUN function applied to projected data matrix
  #' @param ... any additional arguments are passed to `FUN` after the transformed data matrix
  #' @param space select desired subspace (space, complement, both), defaults to "focus"
  #'   - `focus`: orthogonal projection into focus space (default)
  #'   - `complement` orthogonal projection into complement space
  #'   - `both`: transform into orthonormal basis of focus space + complement
  #' @param dim use only these dimensions of the selected subspace
  #' @param byrow.in whether input `M` is a row or column matrix
  #' @param byrow.out whether result of `FUN` is a row or column matrix
  #' @param byrow default value for both `byrow.in` and `byrow.out` (defaults to TRUE)
  #' @return A row (if `byrow.out=TRUE`) or column matrix of vectors returned by `FUN`, transformed back into the original feature space.
  subspace.apply = function (M, FUN, ..., space=c("focus", "both", "complement"), dim=NULL, byrow.in=byrow, byrow.out=byrow, byrow=TRUE) {
    space <- match.arg(space)
    n.in <- if (byrow.in) ncol(M) else nrow(M)
    if (n.in != self$n) stop("dimensionality of data matrix doesn't match GMA space")
    P. <- self$basis(space, dim)
    d <- ncol(P.)
    M1 <- if (byrow.in) M %*% P. else crossprod(P., M)
    res1 <- FUN(M1, ...)
    d.out <- if (byrow.out) ncol(res1) else nrow(res1)
    if (d.out != d) stop("dimensionality of result doesn't match subspace -- did you forget to set byrow.out?")
    if (byrow.out) tcrossprod(res1, P.) else P. %*% res1
  },
  
  #' @description
  #' Identify LDA dimensions in full space or specified subspace.
  #' Note that the dimensions returned are LDA discriminants rather than orthonormal basis vectors.
  #' @param categories vector of categories to be separated by LDA (must match rows of `M`)
  #' @param cohorts optional vector of cohort assignments for repeated-measures LDA (see [`?gma.lda`][gma.lda()])
  #' @param balanced set `balanced=TRUE` to give all groups the same contribution to between-group variance (regardless of group size)
  #' @param R2 what percentage of relative between-group variance needs to be captured by discriminant space
  #' @param M data matrix of row vectors (defaults to internal data matrix of GMA object)
  #' @param space select desired subspace in which LDA is carried out
  #'   - `both`: transform into orthonormal basis of focus space + complement (default)
  #'   - `complement` orthogonal projection into complement space
  #'   - `focus`: orthogonal projection into focus space
  #' @param dim use only these dimensions of the selected subspace
  #' @param idx perform LDA on subset of M selected by `idx`
  #' @return column matrix of LDA discriminants (in coordinates of original feature space)
  discriminant = function (categories, cohorts=NULL, balanced=FALSE, R2=100, M=NULL, space=c("both", "complement", "focus"), dim=NULL, idx=NULL) {
    space <- match.arg(space)
    if (is.null(M)) M <- self$data
    if (nrow(M) != length(categories)) stop("length of categories= must correspond to number of row vectors in M=")
    categories <- as.factor(categories)
    if (!is.null(cohorts)) {
      if (nrow(M) != length(cohorts)) stop("length of cohorts= must correspond to number of row vectors in M=")
      cohorts <- as.factor(cohorts)
    }
    if (!is.null(idx)) {
      M <- M[idx, , drop=FALSE]
      categories <- droplevels(categories[idx])
      if (!is.null(cohorts)) cohorts <- droplevels(cohorts[idx])
    }
    lda.dims <- function (x) gma.lda(x, categories, cohorts=cohorts, balanced=balanced, R2=R2)$scaling
    res <- self$subspace.apply(M, lda.dims, space=space, dim=dim, byrow.out=FALSE)
    rownames(res) <- colnames(M)
    res
  },
  
  #' @description
  #' Extend focus space with LDA dimensions.
  #' @param categories vector of categories to be separated by LDA (must match rows of `M`)
  #' @param cohorts optional vector of cohort assignments for repeated-measures LDA (see [`?gma.lda`][gma.lda()])
  #' @param balanced set `balanced=TRUE` to give all groups the same contribution to between-group variance (regardless of group size)
  #' @param R2 what percentage of relative between-group variance needs to be captured by discriminant space
  #' @param M data matrix of row vectors (defaults to internal data matrix of GMA object)
  #' @param space select desired subspace in which LDA is carried out
  #'   - `complement` orthogonal projection into complement space (default)
  #'   - `both`: transform into orthonormal basis of focus space + complement
  #' @param dim use only these dimensions of the selected subspace
  #' @param idx perform LDA on subset of M selected by `idx`
  #' @param max.dim maximum number of LDA dimensions to add to the focus space
  add.discriminant = function (categories, cohorts=NULL, balanced=FALSE, R2=100, M=NULL, space=c("complement", "both"), idx=NULL, max.dim=Inf) {
    space <- match.arg(space)
    new.dim <- self$discriminant(categories, cohorts=cohorts, balanced=balanced, R2=R2, M=M, space=space, idx=idx)
    if (ncol(new.dim) > max.dim) new.dim <- new.dim[, 1:max.dim, drop=FALSE]
    self$add(new.dim)
  },
  
  #' @description
  #' Compute similarity between the focus spaces of two GMA objects.
  #' @param other a second GMA object based on the same underlying feature space
  #' @param method similarity measure to compute (`dim`, `R2`, or `sigma`); see [?gma.similarity][gma.similarity()] for details
  #' @param add.basis if TRUE, attach matching basis vectors as attributes
  #' @return The selected similarity score, or a vector of cosine similarities between matching dimensions (for `method="sigma"`).
  similarity = function (other, method=c("dimensions", "R2", "sigma"), add.basis=FALSE) {
    method <- match.arg(method)
    if (!inherits(other, "GMA")) stop("other= must be a GMA object")
    if (self$n != other$n) stop("both GMA objects must have the same underlying dimensionality n")
    A <- self$basis()
    B <- other$basis()
    gma.similarity(A, B, method=method, orthogonalize=FALSE, basis=add.basis)
  },
  
  #' @description
  #' Apply a rotation to (selected dimensions of) the focus space. 
  #' In contrast to factor analysis, all GMA rotations are orthogonal transformations that preserve the geometry of the focus space.
  #' @param type selects the rotation to be performed:
  #' - `pca`:    perform a PCA on `M`(or internal data matrix) in the subspace and rotate the dimensions to principal components (default)
  #' - `swap`:   reorder dimensions according to the permutation in `perm` (default: reverse order)
  #' - `flip`:   reverse sign of selected dimensions
  #' - `auto`:   swap and flip dimensions to give best match to those in `basis`
  #' - `turn`:   rotate two basis dimensions of focus space by user-specified angle `phi` (counter-clockwise, note that projected data points rotate clockwise)
  #' - `manual`: rotate basis vectors to match manually specified axes (in `basis`), up to differences in sign
  #'   - `basis` can specify fewer axes than the dimensionality of the subspace, use permutation in `dim` to select which basis vectors are matched
  #'   - axes are projected into focus space (or subspace selected by `dim`) and (re-)orthogonalised
  #' - `align`:  align with subspace in `basis` by finding optimally matching orthonormal basis vectors in both spaces
  #'   - uses [gma.similarity()] algorithm; both subspaces must have the same dimensionality
  #'   - apply twice to align two GMA focus spaces: `A$rotation("align", basis=B)` and `B$rotation("align", basis=A)`
  #' @param dim  subset of the focus space dimensions on which rotation will be performed
  #'             (since the dimensions are orthogonal, rotations on disjoint subsets are completely independent)
  #' @param M optional data matrix for `pca` rotation (defaults to internal data)
  #' @param perm a permutation of the selected focus space dimensions (`1:d`), defaults to `rev(1:d)`
  #' @param basis a column matrix of basis vectors (for `auto`, `manual`, `align`) or a GMA object (whose focus space basis will be used)
  #' @param phi manually specified angle for `turn` rotation (in degrees)
  #' @param tol,check,debug internal use arguments for diagnostics and debugging
  #' @details
  #' CAVEAT: rotations invalidate the user-specified basis vectors (`$axes`), which are replaced with the orthonormal basis vectors `$P`.
  rotation = function (type=c("pca", "swap", "flip", "auto", "turn", "manual", "align"), 
                       dim=NULL, M=NULL, perm=NULL, basis=NULL, phi=0, 
                       tol=1e-6, check=TRUE, debug=FALSE) {
    type <- match.arg(type)
    k <- ncol(self$P) # dimensionality of focus space
    if (k == 0) stop("cannot apply rotation to empty focus space")
    if (is.null(dim)) dim <- seq_len(k)
    stopifnot(length(dim) >= 1, all(1 <= dim & dim <= k), !any(duplicated(dim)))
    P1 <- self$P[, dim, drop=FALSE] # subspace for the rotation (of dimensionality d)
    d <- length(dim)
    if (!is.null(basis) && inherits(basis, "GMA")) {
      basis <- basis$basis("focus")
      if (ncol(basis) > d) basis[, 1:d, drop=FALSE]
    }
    if (type == "pca") {
      if (is.null(M)) M <- self$data
      M1 <- M %*% P1  # projection of M into rotation subspace
      R <- prcomp(M1)$rotation # PCA rotation in the subspace
      if (ncol(R) < d) stop(sprintf("M is singular in %d-dim subspace [%s], PCA returned only %d dimensions", d, paste(dims, collapse=", "), ncol(R)))
      P1.new <- P1 %*% R # principal components in original coordinates form new basis vectors
      ## sign matching against old basis vectors
      cos.mat <- crossprod(P1.new, P1) # find closest original dim to each new dim
      sign.vec <- apply(cos.mat, 1, function (x) if (x[which.max(abs(x))] < 0) -1 else 1) # flip sign if best match has opposite direction
      P1.new <- .scaleMargins(P1.new, cols=sign.vec)
    }
    else if (type == "swap") {
      if (is.null(perm)) perm <- rev(seq_len(d))
      if (!all.equal(sort(perm), seq_len(d))) stop("perm= is not a suitable permutation of the selected dimensions")
      P1.new <- P1[, perm, drop=FALSE]
    }
    else if (type == "flip") {
      P1.new <- -P1
    }
    else if (type == "auto") {
      if (is.null(basis)) stop("basis= must be specified for type='auto'")
      stopifnot(ncol(basis) == d)
      stopifnot(nrow(basis) == nrow(P1))
      P1.new <- matrix(0, nrow=nrow(P1), ncol=d)
      for (i in seq_len(d)) {
        cos.vec <- crossprod(basis[, i, drop=FALSE], P1) # cosine between basis[i] and remaining P1
        j <- which.max(abs(cos.vec)) # find best match 
        sign.adj <- if (cos.vec[j] < 0) -1 else 1
        P1.new[, i] <- sign.adj * P1[, j]
        P1 <- P1[, -j, drop=FALSE]
      }
    }
    else if (type == "turn") {
      if (d != 2) stop("dim= must select exactly two dimensions of focus space for type='turn'")
      P1.new <- P1 %*% rbind(c(cospi(phi/180), -sinpi(phi/180)), c(sinpi(phi/180), cospi(phi/180)))
    }
    else if (type == "manual") {
      if (is.null(basis)) stop("basis= must be specified for type='manual'")
      m <- ncol(basis)
      stopifnot(1 <= m && m <= d)
      stopifnot(nrow(basis) == nrow(P1))
      basis <- P1 %*% t(P1) %*% basis # project basis into target subsapce
      basis <- gma.orthogonalize(basis, fail=FALSE) # may need to re-orthonormalise after projection (also checks rank)
      if (is.null(basis)) stop("basis= is not linearly independent after projection into focus space")
      P1.new <- P1
      for (k in seq_len(min(m, d - 1))) {
        # use planar rotation to align P1.new[, k] with basis[, k], applied to P1.new[, k:d]
        u <- P1.new[, k, drop=FALSE] # u, v will become orthonormal basis for rotation plane
        v <- basis[, k, drop=FALSE]
        cos.phi <- sum(v * u) # angle between basis vector k and user axis k
        phi <- 180 * acos(pmin(1, pmax(-1, cos.phi))) / pi # for debugging output only
        if (debug) cat(sprintf("%d) rotation angle phi = %.2f deg\n", k, phi))
        if ((1 - abs(cos.phi)) < tol) {
          warning(sprintf("(axis, basis) pair #%d too close to collinearity at %.2f degrees, can't rotate", k, phi))
          P1.new[, k] <- sign(cos.phi) * v # force basis[k] to axis[k] (or -axis[k] to avoid sign flip)
          P1.new[, k:d] <- gma.orthogonalize(P1.new[, k:d]) # re-orthogonalise basis vectors instead of rotation
        }
        else {
          v <- (v - cos.phi * u)    # Gram-Schmidt orthogonalization
          v <- v / sqrt(sum(v * v)) # re-normalise, same as v / sin(phi)
          if (check) gma.is.orthogonal(cbind(u, v), fail=TRUE)
          sin.phi <- sin(acos(cos.phi)) # sin.phi <- sqrt(1 - cos.phi^2) is almost identical
          R <- diag(nrow(P1)) + sin.phi * (tcrossprod(v, u) - tcrossprod(u, v)) + (cos.phi - 1) * (tcrossprod(u) + tcrossprod(v)) # rotation matrix R
          if (check) gma.is.orthogonal(R, fail=TRUE)
          P1.new[, k:d] <- R %*% P1.new[, k:d]
          if (check || debug) {
            u <- P1.new[, k, drop=FALSE]  # compare basis and axis vector after rotation (should be identical)
            v <- basis[, k, drop=FALSE]
            d.uv <- min(sum((u - v)^2), sum((u + v)^2)) # we expect u = v or u = -v (if rotation could not be performed)
            if (debug) cat(sprintf("   | b[%d] - a[%d] |^2 = %.6f\n", k, k, d.uv))
            if (check) stopifnot(d.uv < tol)
          }
        }
        if (check || debug) {
          P1.sim <- gma.similarity(P1.new, P1) # ensure that focus space hasn't been changed by rotation
          if (debug) cat(sprintf("   preservation of focus space: lost %g dims\n", d - P1.sim))
          if (check) stopifnot((d - P1.sim) < tol)
        }
      }
    }
    else if (type == "align") {
      if (is.null(basis)) stop("basis= must be specified for type='align'")
      m <- ncol(basis)
      if (m != d) stop("basis= must have excactly same dimensionality as selected subspace for type='align'")
      res <- gma.similarity(P1, basis, method="sigma", orthogonalize=TRUE, basis=TRUE)
      P1.new <- attr(res, "basisA")
      if (debug) {
        cat(sprintf("Match quality of aligned basis:\n"))
        cat(sprintf(" - basis[%d]:  sigma = %.3f\n", dim, res), sep="")
      }
      if (check) {
        P1.sim <- gma.similarity(P1.new, P1) # ensure that focus space hasn't been changed
        if (debug) cat(sprintf(" - preservation of focus space: lost %g dims\n", d - P1.sim))
        if (check) stopifnot((d - P1.sim) < tol)
      }
    }
    else {
      stop("not yet implemented")
    }
    self$P[, dim] <- P1.new
    self$axes <- self$P # user-specified basis is invalidated, overwrite with orthonormal
    self$rotated <- TRUE
    invisible(self)
  }
  
))

