## 3D scatterplot optimised for GMA studies

## - must provide matrix M and corresponding data frame with metadata (Meta)
## - parameters select, pch, col, size are evaluated within data frame Meta
## - if Meta=NULL, a dummy data frame without information is provided
## - expressions for connect.select, connect.lwd and connect.col can access metadata of start and end points as start$<var> and end$<var>
## - if connect.arrow=TRUE, points are connected by arrows rather than simple lines (expensive)
## - arrow width is automatically selected, but can be overridden by setting connect.arrow=<width>; in either case, it is further scaled by connect.lwd
## - if legend=TRUE, display legend on 2D graphics device, or write directly to PDF file if legend.pdf= is specified
#' @export
gma.pairs3d <- function (M, dims=c(2,1,3), Meta=NULL, select=NULL, pch=NULL, col=NULL, 
                         pch.vals=c("circles","triangles","squares"), pch.cols=NULL, 
                         col.vals=corpora.palette("simple"), linecol.vals=rev(corpora.palette("simple")), 
                         size=.05, legend.cex=3, legend.magnify=2, legend=FALSE, bbox=TRUE,
                         connect=NULL, connect.select=NULL, connect.lwd=1, connect.col=NULL, connect.arrow=FALSE, ...) {
  .ensure.library(c("rgl"), "gma.pairs3d()")
  if (length(dims) != 3) stop("need exactly 3 dimensions for 3D visualization")
  size <- rep(size, length.out=nrow(M))
  pt.size <- max(size) # assumed point size for adjustments
  
  select.expr <- substitute(select) # will evaluate arguments in context of metadata
  pch.expr <- substitute(pch)
  col.expr <- substitute(col)
  connect.select.expr <- substitute(connect.select)
  connect.col.expr <- substitute(connect.col)
  
  if (is.null(Meta)) {
    item.ids <- rownames(M)
    if (is.null(item.ids)) item.ids <- sprintf("%d04d", 1:nrow(M))
    Meta <- data.frame(id=item.ids) # provide dummy metadata
  } else {
    if (nrow(Meta) != nrow(M)) stop("metadata table Meta must have the same number of rows in the same order as the data matrix M")
  }
  
  select <- eval(select.expr, Meta, parent.frame())
  if (is.logical(select)) select <- which(select) # ensure that subset index is given as line numbers (NULL is not affected)
  pch <- eval(pch.expr, Meta, parent.frame())
  if (length(pch) == 1) pch <- rep(pch, nrow(M))  # NULL has length 0 and is not affected
  col <- eval(col.expr, Meta, parent.frame())
  if (length(col) == 1) col <- rep(col, nrow(M))
  
  if (!is.null(connect)) {
    if (!( is.matrix(connect) && ncol(connect) == 2 )) stop("connect= argument must be 2-column matrix")
    meta.envir <- list(from=Meta[connect[,1], ], to=Meta[connect[,2], ])
    connect.select <- eval(connect.select.expr, meta.envir, parent.frame())
    connect.col <- eval(connect.col.expr, meta.envir, parent.frame())
    if (length(connect.col) == 1) connect.col <- rep(connect.col, nrow(connect)) # NULL has length 0 and is not affected
    if (is.character(connect.col)) connect.col <- factor(connect.col)
    
    if (!is.null(connect.select)) {
      ## reduce to selected subset of connections (and adjust vector of colour values)
      connect <- connect[connect.select, , drop=FALSE]
      if (!is.null(connect.col)) connect.col <- connect.col[connect.select]
    }
    
    ## look up coordinates of start and end points before row indices may be changed by select
    P1 <- M[connect[,1], dims, drop=FALSE]
    P2 <- M[connect[,2], dims, drop=FALSE]
  }
  
  if (!is.null(pch)) {
    stopifnot(length(pch) == nrow(M))
    if (is.character(pch)) pch <- factor(pch)
  }
  if (!is.null(col)) {
    stopifnot(length(col) == nrow(M))
    if (is.character(col)) col <- factor(col)
  }
  if (!is.null(select)) {
    ## reduce to selected subset of points (adjusting point symbols and colour values)
    M <- M[select, , drop=FALSE]
    if (!is.null(pch)) pch <- pch[select]
    if (!is.null(col)) col <- col[select]
    
    if (!is.null(connect)) {
      ## drop connections if their start or end points are no longer displayed
      idx.keep <- connect[,1] %in% select & connect[,2] %in% select
      P1 <- P1[idx.keep, , drop=FALSE]
      P2 <- P2[idx.keep, , drop=FALSE]
      connect <- connect[idx.keep, , drop=FALSE]
      if (!is.null(connect.col)) connect.col <- connect.col[idx.keep]
    }
  }
  
  par3d(skipRedraw=TRUE)
  clear3d()
  bg3d() # reset background

  leg <- list(pch=FALSE, col=FALSE, connect=FALSE)
  if (is.null(pch)) {
    pch <- rep(pch.vals[1], nrow(M))
  } else if (!is.factor(pch)) {
    pch <- pch.vals[pch]
  } else {
    pch <- droplevels(pch)
    pch.levels <- levels(pch)
    pch.codes <- as.integer(pch)
    pch.n <- length(pch.levels)
    if (pch.n > length(pch.vals)) stop(sprintf("too few point shapes supplied (need at least %d)", pch.n))
    pch <- pch.vals[pch.codes]
    leg$pch <- TRUE
  }
  
  if (is.null(col)) {
    col <- rep(col.vals[1], nrow(M))
  } else if (!is.factor(col)) {
    col <- col.vals[col]
  } else {
    col <- droplevels(col)
    col.levels <- levels(col)
    col.codes <- as.integer(col)
    col.n <- length(col.levels)
    if (col.n > length(col.vals)) stop(sprintf("too few point colours supplied (need at least %d)", col.n))
    col <- col.vals[col.codes]
    leg$col <- TRUE
  }
  
  for (sym in levels(factor(pch))) {
    sym.idx <- pch == sym
    .draw.3dpoints(M[sym.idx, dims], size=size[sym.idx], shape=sym, color=col[sym.idx])
  }
  
  if (!is.null(connect)) {
    ## add connections between points (lines or arrows)
    
    if (is.null(connect.col)) {
      connect.col <- rep(linecol.vals[1], nrow(connect))
    } else if (!is.factor(connect.col)) {
      connect.col <- linecol.vals[connect.col]
    } else {
      connect.col <- droplevels(connect.col)
      connect.levels <- levels(connect.col)
      connect.codes <- as.integer(connect.col)
      connect.n <- length(connect.levels)
      if (connect.n > length(linecol.vals)) stop(sprintf("too few line colours supplied (need at least %d)", connect.n))
      connect.col <- linecol.vals[connect.codes]
      leg$connect <- TRUE
    }
    
    if (connect.arrow) {
      ## connect points by arrows (width scaled by line width)
      d <- sqrt(rowSums((P1 - P2)^2)) # arrow lengths
      arrow.width <- (2 * pt.size) / 3   # auto-scale arrow width based on point size (width of shaft = 1/3 diameter of point)
      if (is.numeric(connect.arrow)) arrow.width <- connect.arrow
      arrow.width <- arrow.width * connect.lwd
      for (i in 1:nrow(connect)) {
        p1 <- P1[i,]
        p2 <- P2[i,]
        .arrow3d((p1 + p2)/2, p2 - p1, d[i] - 2 * pt.size, width=arrow.width, adj=0.5, col=connect.col[i], tipex=3, sides=32) # NB: arrows shortened by point size (radius) on either end
      }
    } else {
      ## connect points by lines with specified line width
      P.pairs <- matrix(as.vector(rbind(t(P1), t(P2))), ncol=3, byrow=TRUE) # interleave P1 and P2 by rows (must have 3 columns each)
      segments3d(P.pairs, lwd=connect.lwd, color=rep(connect.col, each=2), line_antialias=TRUE)
    }
  }
  
  if (legend) {
    bgplot3d({
      par(mar=c(0, 0, 0, 0))
      plot(0, 0, type="n", xlim=0:1, ylim=0:1, xaxs="i", yaxs="i", axes=FALSE, bty="n")
      if (leg$pch) {
        sym.names <- c("circles", "triangles", "squares")
        sym.codes <- c(16, 17, 15)
        pch.legend <- sym.codes[ pmatch(pch.vals[1:pch.n], sym.names, duplicates.ok=TRUE) ]
        col.legend <- if (!is.null(pch.cols)) pch.cols[1:pch.n] else "black"
        legend("topleft", bty="n", legend=pch.levels, pch=pch.legend, cex=legend.cex, col=col.legend)
      }
      if (leg$col) {
        legend("topright", bty="n", legend=col.levels, fill=col.vals[1:col.n], cex=legend.cex)
      }
      if (leg$connect) {
        legend("bottomright", bty="n", legend=connect.levels, col=linecol.vals[1:connect.n], lwd=4, cex=legend.cex)
      }
    }, magnify=legend.magnify, leg=leg)    
  }
  
  if (bbox) {
    bbox3d(color=c("grey30", "black"), shininess=100, emission="grey60")
    for (side in c("z-", "z+", "y-", "y+")) grid3d(side, col="#999999", lwd=.5)
    for (side in c("x-", "x+")) grid3d(side, col="#777777", lwd=.5)
  } else {
    axes3d(c("x", "y", "z"))
    if (!is.null(colnames(M))) {
      axisnames <- colnames(M)[dims]
      title3d(xlab=axisnames[1], ylab=axisnames[2], zlab=axisnames[3])
    }
    ## decorate3d(box=FALSE, axes=TRUE) # -- open box, similar to standard 2D plots
  }
  
  aspect3d("iso")
  par3d(skipRedraw=FALSE)
}

## ----- helper functions -----

## 3-dimensional scatterplot with standard shapes; only single shape allowed per call;
## all other grqphics parameters (such as col=) are passed throuh
.draw.3dpoints <- function (xyz, shape=c("circles", "triangles", "squares"), size=1, ...) {
  shape <- match.arg(shape)
  if (shape == "circles") {
    spheres3d(xyz, radius=size, ...)
  } else if (shape == "triangles") {
    shapelist3d(tetrahedron3d(), xyz, size=size, ...)  
  } else { # shape == "squares"
    shapelist3d(cube3d(), xyz, size=size, ...)  
  }
}

## add a single solid arrow to the 3D plot
##  - arguments: loc(ation), dir(ection), len(gth), width, adj(ustment of anchor point)
.arrow3d <- function (loc, dir, length, width=length/40, adj=0.5, col="grey40", tipex=2, sides=32, ignoreExtent=TRUE) {
  dir <- as.vector(dir)
  dir.length <- sqrt(sum(dir^2))
  dir <- dir / dir.length # normalized direction vector
  tip.width <- width * tipex
  tip.length <- width / sin(pi/9) # 20 degrees from center line
  p1 <- loc - adj * length * dir  # tail of arrow
  p2 <- loc + ((1 - adj) * length - tip.length) * dir # start of tip
  p3 <- loc + (1 - adj) * length * dir # end of tip
  if (ignoreExtent) {
    par.save <- par3d(ignoreExtent=TRUE)
    on.exit(par3d(par.save), add=TRUE)
  }
  shade3d(cylinder3d(rbind(p1, p2), radius=width/2, sides=sides, closed=-1), col=col)
  shade3d(cylinder3d(rbind(p2, p3), radius=c(tip.width/2, 1e-6), sides=sides, closed=-1), col=col)  
}

##  - use this variant to connect start and end point (loc, dir and length are computed automatically)
.connect3d <- function (from, to, width=NULL, ...) {
  dir <- to - from
  l <- sqrt(sum(dir^2))
  if (is.null(width)) width <- l/40
  arrow3d(from, dir, l, width=width, adj=0, ...)
}

