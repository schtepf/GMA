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
