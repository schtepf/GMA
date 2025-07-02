## Scatterplot matrix display optimised for GMA studies

#' Scatterplot Matrix as Visualisation for GMA
#'
#' **TODO:**
#' write a description
#'
#' @param M a row matrix of data points to be plotted
#' @param dims which dimensions of `M` to include in the scatterplot matrix (default: all)
#' @param Meta a data frame of metadata corresponding to the data points in `Meta`
#' @param select expression to select a subset of data points to be plotted, either as logical vector or vector of row numbers.
#'        The expression is evaluated within `Meta`, so its columns can directly be used as variables. Note that `select` only
#'        reduces data points to a subset but does not reorder them.
#' @param alpha.select if set to a positive value in the range \eqn{(0, 1]}, data points omitted by a `select` expression are not
#'        hidden completely, but their opacity is reduced to the specified value
#' @param pch expression to determine plot symbols for the data points, usually based on columns of `Meta` (which can be accessed directly as variables)
#' @param col expression to determine colours for the data points, usually based on columns of `Meta` (which can be accessed directly as variables)
#' @param pch.vals plot symbols used for the levels of the expression `pch` (values are looked up by level index, not by level name)
#' @param pch.cols option colours for the levels of the expression `pch`, so plot symbols can appear in different colours (cannot be combined with `col`)
#' @param col.vals point colours used for the levels of the expression `col` (values are looked up by level index, not by level name)
#' @param cex global magnification factor for plot symbols and text
#' @param legend.cex magnification factor for legend boxes (not relative to `cex`, but default is `1.4 * cex`)
#' @param randomize if TRUE (the default), data points are plotted in random order so no categories are hidden.
#'        Specify a positive integer to set a random seed for reproducibility, or set to FALSE to keep original plotting order.
#' @param gap distance between panels of the scatter plot matrix, in margin lines (default: 0.5)
#' @param oma adjusted outer margins of individual panels for a more compact display (internal use)
#' @param iso if TRUE, adjust ranges for all axes so that the scatterplots are isometric (provided that the overall plot has square aspect).
#'        That is, ranges have the same width on each axis, but are shifted to center the data.
#' @param compact if TRUE, only the upper triangle of scatterplots is shown and legends are moved to the lower triangle, for a more compact display
#' @param lim use `lim=c(min, max)` to specify a fixed range for all coordinate axes,
#'        or choose individual ranges as an \eqn{n_{\text{dim}} \times 2}{n_dim x 2} matrix (overrides `iso=TRUE`)
#' @param ... further graphics parameters are passed through to the underlying plot functions
#' @export
#' @importFrom corpora alpha.col corpora.palette
#' @importFrom graphics Axis box legend mtext par plot points strwidth text
#' @importFrom grDevices dev.flush dev.hold
gma.pairs <- function (M, dims=NULL, Meta=NULL, select=NULL, alpha.select=0, pch=NULL, col=NULL,
                       pch.vals=1:10, pch.cols=NULL, col.vals=corpora.palette("simple"),
                       cex=1, legend.cex=1.4*cex, randomize=TRUE, gap=.5,
                       oma=c(2,2,2,2), iso=FALSE, compact=FALSE, lim=NULL, ...) {
  if (is.null(dims)) dims <- seq_len(ncol(M))
  n.dim <- length(dims)
  stopifnot(n.dim >= 2)
  if (!all(dims %in% seq_len(ncol(M)))) stop("invalid dimensions selected")

  if (!is.null(lim)) {
    if (is.matrix(lim)) {
      if (nrow(lim) != n.dim || ncol(lim) != 2) stop(sprintf("lim= must be a %d x 2 matrix or a vector c(min, max)", n.dim))
    }
    else {
      if (length(lim) != 2) stop(sprintf("lim= must be a %d x 2 matrix or a vector c(min, max)", n.dim))
      lim <- cbind(rep(lim[1], n.dim), rep(lim[2], n.dim))
    }
  }

  select.expr <- substitute(select) # evaluate arguments in context of metadata
  pch.expr <- substitute(pch)
  col.expr <- substitute(col)

  if (is.null(Meta)) {
    item.ids <- rownames(M)
    if (is.null(item.ids)) item.ids <- sprintf("%d04d", 1:nrow(M))
    Meta <- data.frame(id=item.ids) # provide dummy metadata
  } else {
    if (nrow(Meta) != nrow(M)) stop("metadata table Meta must have the same number of rows in the same order as the data matrix M")
  }

  select <- eval(select.expr, Meta, parent.frame())
  if (is.numeric(select)) {
    ## ensure that subset index is given as a logical vector (not as line numbers)
    stopifnot(all(1 <= select), all(select <= nrow(Meta)))
    select <- seq_len(nrow(Meta)) %in% select
  }
  pch <- eval(pch.expr, Meta, parent.frame())
  if (length(pch) == 1) pch <- rep(pch, nrow(M))  # NULL has length 0 and is not affected
  col <- eval(col.expr, Meta, parent.frame())
  if (length(col) == 1) col <- rep(col, nrow(M))

  if (!is.null(pch)) {
    stopifnot(length(pch) == nrow(M))
    if (is.character(pch)) pch <- factor(pch)
  }
  if (!is.null(col)) {
    stopifnot(length(col) == nrow(M))
    if (is.character(col)) col <- factor(col)
  }
  if (!is.null(select) && !(alpha.select > 0)) {
    ## reduce to selected subset of points (unless we just fade out non-selected points)
    M <- M[select, , drop=FALSE]
    if (!is.null(pch)) pch <- pch[select]
    if (!is.null(col)) col <- col[select]
  }
  labels <- if (compact || is.null(colnames(M))) rep("", n.dim) else colnames(M)[dims]

  k <- if (compact && n.dim > 2) 2 else 1
  if (is.null(pch)) {
    pch <- rep(pch.vals[1], nrow(M))
  } else {
    pch <- if (is.numeric(pch)) factor(pch, levels=1:max(pch)) else as.factor(pch)
    pch.levels <- levels(pch)
    pch.codes <- as.integer(pch)
    pch.n <- length(pch.levels)
    pch.ncol <- if (pch.n > 5) 2 else 1
    stopifnot(length(pch.vals) >= pch.n)
    pch <- pch.vals[pch.codes]
    if (!is.null(pch.cols)) {
      stopifnot(length(pch.cols) >= pch.n)
    } else {
      pch.cols <- rep("black", pch.n)
    }
    labels[k] <- ":pch:"; k <- k+1
  }

  if (is.null(col)) {
    col <- rep(col.vals[1], nrow(M))
  } else {
    col <- if (is.numeric(col)) factor(col, levels=1:max(col)) else as.factor(col)
    col.levels <- levels(col)
    col.codes <- as.integer(col)
    col.n <- length(col.levels)
    col.ncol <- ceiling(sqrt(col.n / 5)) # 1 col up to 5, 2 cols up to 20, 3 cols up to 45
    stopifnot(length(col.vals) >= col.n)
    col <- col.vals[col.codes]
    if (compact && n.dim == 3 && k > 2) {
      labels[2] <- ":pch:col:" # special case for compact layout
    } else {
      labels[k] <- ":col:"; k <- k+1
    }
  }

  if (!is.null(select) && alpha.select > 0) {
    ## fade out non-selected points by reducing their opacity
    col[!select] <- alpha.col(col[!select], alpha.select)
  }

  my.panel <- function (x, y, label, ...) {
    if (label == ":pch:") legend(x, y, xjust=0.5, yjust=0.5, legend=pch.levels, pch=pch.vals[1:pch.n], col=pch.cols[1:pch.n], pt.cex=1.2, pt.lwd=1.2, cex=legend.cex, ncol=pch.ncol, bty="n") else
      if (label == ":col:") legend(x, y, xjust=0.5, yjust=0.5, legend=col.levels, fill=col.vals[1:col.n], border=col.vals[1:col.n], cex=legend.cex, ncol=col.ncol, bty="n") else
        if (label == ":pch:col:") {
          # legend(x, y, xjust=0.5, yjust=1.1, legend=pch.levels, pch=pch.vals[1:pch.n], col=pch.cols[1:pch.n], pt.cex=1.2, pt.lwd=1.2, cex=legend.cex, ncol=pch.ncol)
          # legend(x, y, xjust=0.5, yjust=-0.1, legend=col.levels, fill=col.vals[1:col.n], border=col.vals[1:col.n], cex=legend.cex, ncol=col.ncol)
          legend("top", inset=.05, legend=pch.levels, pch=pch.vals[1:pch.n], col=pch.cols[1:pch.n], pt.cex=1.2, pt.lwd=1.2, cex=legend.cex, ncol=pch.ncol, bty="n")
          legend("bottom", inset=.05, legend=col.levels, fill=col.vals[1:col.n], border=col.vals[1:col.n], cex=legend.cex, ncol=col.ncol, bty="n")
        } else
          text(x, y, label, adj=c(0.5, 0.5), cex=legend.cex)
  }

  if (randomize) {
    if (is.numeric(randomize)) set.seed(randomize)
    nR <- nrow(M)
    idx <- sample.int(nR)
    pch <- pch[idx]
    col <- col[idx]
    M <- M[idx, , drop=FALSE]
  }

  upper.panel <- points
  lower.panel <- if (compact) function (x, y, ...) {} else upper.panel
  .pairsCompact(M[, dims], pch=pch, col=col, cex=cex,
                upper.panel=upper.panel, lower.panel=lower.panel, text.panel=my.panel,
                labels=labels, gap=gap, oma=oma, compact=compact, iso=iso, lim=lim, ...)
}

## modified version of built-in pairs() function for scatterplot matrix
##  - row1attop = FALSE is not supported
##  - new option compact=TRUE to omit leftmost column and bottom row from plot
##  - new option iso=TRUE to force isometric scaling of all axes (but perhaps with different center)
##  - new option lim=matrix(ndim, 2) to specify explicit axis limits (overrides iso=TRUE)
.pairsCompact <- function (x, labels, panel = points, ...,
                           horInd = 1:nc, verInd = 1:nc,
                           lower.panel = panel, upper.panel = panel,
                           diag.panel = NULL, text.panel = textPanel,
                           label.pos = 0.5 + has.diag/3, line.main = 3,
                           cex.labels = NULL, font.labels = 1,
                           row1attop = TRUE, gap = 1, log = "",
                           compact=FALSE, iso=FALSE, lim=NULL)
{
  if (!missing(horInd) || !missing(verInd)) stop("this implementation doesn't support modified horInd or verInd")
  if(doText <- missing(text.panel) || is.function(text.panel))
    textPanel <-
      function(x = 0.5, y = 0.5, txt, cex, font)
        text(x, y, txt, cex = cex, font = font)

  localAxis <- function(side, x, y, xpd, bg, col=NULL, main, oma, ...) {
    ## Explicitly ignore any color argument passed in as
    ## it was most likely meant for the data points and
    ## not for the axis.
    xpd <- NA
    if(side %% 2L == 1L && xl[j]) xpd <- FALSE
    if(side %% 2L == 0L && yl[i]) xpd <- FALSE
    if(side %% 2L == 1L) Axis(x, side = side, xpd = xpd, ...)
    else Axis(y, side = side, xpd = xpd, ...)
  }

  localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main)
    lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main)
    upper.panel(...)

  localDiagPanel <- function(..., main, oma, font.main, cex.main)
    diag.panel(...)

  dots <- list(...); nmdots <- names(dots)
  if (!is.matrix(x)) {
    x <- as.data.frame(x)
    for(i in seq_along(names(x))) {
      if(is.factor(x[[i]]) || is.logical(x[[i]]))
        x[[i]] <- as.numeric(x[[i]])
      if(!is.numeric(unclass(x[[i]])))
        stop("non-numeric argument to 'pairs'")
    }
  } else if (!is.numeric(x)) stop("non-numeric argument to 'pairs'")
  panel <- match.fun(panel)
  if((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))
    lower.panel <- match.fun(lower.panel)
  if((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))
    upper.panel <- match.fun(upper.panel)
  if((has.diag  <- !is.null( diag.panel)) && !missing( diag.panel))
    diag.panel <- match.fun( diag.panel)

  if(row1attop) {
    tmp <- lower.panel; lower.panel <- upper.panel; upper.panel <- tmp
    tmp <- has.lower; has.lower <- has.upper; has.upper <- tmp
  }

  nc <- ncol(x)
  if (nc < 2L) stop("only one column in the argument to 'pairs'")
  if(!all(horInd >= 1L & horInd <= nc))
    stop("invalid argument 'horInd'")
  if(!all(verInd >= 1L & verInd <= nc))
    stop("invalid argument 'verInd'")
  if(doText) {
    if (missing(labels)) {
      labels <- colnames(x)
      if (is.null(labels)) labels <- paste("var", 1L:nc)
    }
    else if(is.null(labels)) doText <- FALSE
  }
  if (compact) {
    verInd <- if (row1attop) verInd[1:(nc-1)] else verInd[-1]
    horInd <- horInd[-1]
  }
  oma <- if("oma" %in% nmdots) dots$oma
  main <- if("main" %in% nmdots) dots$main
  if (is.null(oma))
    oma <- c(4, 4, if(!is.null(main)) 6 else 4, 4)
  opar <- par(mfrow = c(length(horInd), length(verInd)),
              mar = rep.int(gap/2, 4), oma = oma)
  on.exit(par(opar))
  dev.hold(); on.exit(dev.flush(), add = TRUE)

  if (iso) axis.width <- max(apply(x, 2, function (y) diff(range(y)))) # maximal width of range on an axis
  if (!is.null(lim) && (nrow(lim) != nc || ncol(lim) != 2)) stop("lim= must be a n_dim x 2 matrix")

  xl <- yl <- logical(nc)
  if (is.numeric(log)) xl[log] <- yl[log] <- TRUE
  else {xl[] <- grepl("x", log); yl[] <- grepl("y", log)}
  for (i in if(row1attop) verInd else rev(verInd))
    for (j in horInd) {
      l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y", ""))
      if (iso || !is.null(lim)) {
        if (!is.null(lim)) {
          iso.xlim <- lim[j, ] # explicit axis limits take precedence
          iso.ylim <- lim[i, ]
        }
        else {
          iso.xlim <- mean(range(x[, j])) + c(-0.5, 0.5) * axis.width # determine isometric axis limits
          iso.ylim <- mean(range(x[, i])) + c(-0.5, 0.5) * axis.width
        }
        localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, type = "n", xlim=iso.xlim, ylim=iso.ylim, ..., log = l)
      } else {
        ## plot with automatic axis limits if neither iso=TRUE nor lim= have been specified
        localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, type = "n", ..., log = l)
      }
      if(i == j || (i < j && has.lower) || (i > j && has.upper) ) {
        box()
        if(i == min(verInd)  && (!(j %% 2L) || !has.upper || !has.lower ))
          localAxis(1L + 2L*row1attop, x[, j], x[, i], ...)
        if(i == max(verInd) && (  j %% 2L  || !has.upper || !has.lower ))
          localAxis(3L - 2L*row1attop, x[, j], x[, i], ...)
        if(j == min(horInd)  && (!(i %% 2L) || !has.upper || !has.lower ))
          localAxis(2L, x[, j], x[, i], ...)
        if(j == max(horInd) && (  i %% 2L  || !has.upper || !has.lower ))
          localAxis(4L, x[, j], x[, i], ...)
        mfg <- par("mfg")
        if(i == j) {
          if (has.diag) localDiagPanel(as.vector(x[, i]), ...)
          if (doText) {
            par(usr = c(0, 1, 0, 1))
            if(is.null(cex.labels)) {
              l.wid <- strwidth(labels, "user")
              cex.labels <- max(0.8, min(2, .9 / max(l.wid)))
            }
            xlp <- if(xl[i]) 10^0.5 else 0.5
            ylp <- if(yl[j]) 10^label.pos else label.pos
            text.panel(xlp, ylp, labels[i],
                       cex = cex.labels, font = font.labels)
          }
        } else if(i < j)
          localLowerPanel(as.vector(x[, j]), as.vector(x[, i]), ...)
        else
          localUpperPanel(as.vector(x[, j]), as.vector(x[, i]), ...)
        if (any(par("mfg") != mfg))
          stop("the 'panel' function made a new plot")
      } else par(new = FALSE)

    }
  if (!is.null(main)) {
    font.main <- if("font.main" %in% nmdots) dots$font.main else par("font.main")
    cex.main <- if("cex.main" %in% nmdots) dots$cex.main else par("cex.main")
    mtext(main, 3, line.main, outer=TRUE, at = 0.5, cex = cex.main, font = font.main)
  }
  invisible(NULL)
}
