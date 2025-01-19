## Internal helper functions not meant for export

## Inefficient substitute for scaleMargins() in order to avoid dependency on 'wordspace' package
.scaleMargins <- function (M, rows=NULL, cols=NULL) {
  if (!(is.matrix(M) && is.numeric(M))) stop("M= must be a numeric matrix")
  nr <- nrow(M)
  nc <- ncol(M)
  if (!is.null(rows)) {
    if (length(rows) == 1) rows <- rep(rows, nr)
    if (length(rows) != nr) stop("rows= must either be a scalar or conformable with the rows of M=")
    M <- M * rows
  }
  if (!is.null(cols)) {
    if (length(cols) == 1) cols <- rep(cols, nc)
    if (length(cols) != nc) stop("cols= must either be a scalar or conformable with the columns of M=")
    M <- t(t(M) * cols) # shouldn't be much worse than using sweep()
  }
  M
}

## Check that argument is a numeric matrix, optionally promoting vectors to column vectors
.ensure.matrix <- function (M, message="argument", vector.ok=FALSE, minR=0, minC=0) {
  if (is.vector(M) && vector.ok) M <- t(t(M)) # promote to column vector, preserving names
  if (!(is.matrix(M) && is.numeric(M))) stop(sprintf("%s must be a numeric matrix", message))
  if (nrow(M) < minR) stop(sprintf("%s must have >= %d rows", message, minR))
  if (ncol(M) < minC) stop(sprintf("%s must have >= %d columns", message, minC))
  invisible(M)
}

## Load optional packages, aborting if they're not available
.ensure.library <- function (pkgs, message="some 'gmatools' function") {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly=TRUE)) stop(sprintf("%s requires package '%s' to be installed", message, p))
  }
}
