## Barplot of feature weights for focus space dimensions

##' Barplot of feature weights for one or more focus space dimensions
##' 
##' **TODO DESCRIPTION**
##' The **ggplot2** package must be installed in order to use this function.
##' 
##' @param basis column matrix of basis vectors for the focus space dimensions
##' @param dim dimensions to be shown in plot (default: all). If more than one dimension is selected, the barplots will be stacked vertically.
##' @param idx optional index vector selecting a subset of feature weights to be displayed
##' @param names,feature.names = labels for dimensions and features in the plot. Labels must be provided for selected dimensions only,
##'        but for all features (before selection with `idx`).
##' @param ylim display range on y-axis (NULL or two-element vector)
##' @param zlim range of values mapped to the colour scale, with out-of-band values clamped to maximal colours.
##'        It is always ensured that the colour range is symmetric around 0 so that zero weights are white.
##' @param ylab label for y-axis of the plot
##' @param main optional main title for the plot
##' @param bw if TRUE, produce a black & white version of the plot
##' @param base_size basic font size for the plot (see [ggplot2::theme_bw()] documentation)
##' @returns A `ggplot` object that has to be printed to display the graph.
##' @export
gma.plot.weights <- function (basis, dim=seq_len(ncol(basis)), idx=NULL, names=NULL, feature.names=NULL, 
                              ylab="normalized feature weights", main=NULL, bw=FALSE, base_size=12, ylim=NULL, zlim=NULL) {
  .ensure.library(c("ggplot2", "scales", "rlang"), "gma.plot.weights()")
  basis <- .ensure.matrix(basis, vector.ok=TRUE)
  n <- nrow(basis)
  k <- length(dim)
  stopifnot(all(dim >= 1 & dim <= ncol(basis)))
  if (is.null(names)) names <- if (!is.null(colnames(basis))) colnames(basis)[dim] else paste("Dim", dim)
  if (is.null(feature.names)) feature.names <- if (!is.null(rownames(basis))) rownames(basis) else sprintf("feature #%d", 1:n)
  stopifnot(length(names) == k)
  stopifnot(length(feature.names) == n)
  if (!is.null(idx)) {
    basis <- basis[idx, , drop=FALSE]
    feature.names <- feature.names[idx]
    n <- nrow(basis)
  }
  
  .data <- NULL # silence R CHECK warnings without listing rlang in Imports:
  info.tbl <- data.frame(feature=factor(rep(feature.names, k), levels=feature.names), 
                         dimension=factor(rep(names, each=n), levels=names), 
                         weight=as.vector(basis[, dim]))
  bp <- ggplot2::ggplot(info.tbl, ggplot2::aes(x=.data$feature, y=.data$weight)) + 
    ggplot2::facet_grid(dimension ~ .) + 
    ggplot2::geom_bar(ggplot2::aes(fill=.data$weight), stat="identity", position="identity")
  bp <- bp + ggplot2::ylab(ylab) + ggplot2::xlab("")
  if (!is.null(main)) bp <- bp + ggplot2::ggtitle(main)
  bp <- bp + 
    ggplot2::theme_bw(base_size=base_size) + 
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1))
  if (bw) {
    bp <- bp + ggplot2::scale_fill_gradient2(low="black", mid="white", high="black", limits=zlim, oob=scales::squish) 
  } else {
    bp <- bp + ggplot2::scale_fill_gradient2(low="#CC0000", mid="white", high="#009900", limits=zlim, oob=scales::squish)
  }
  if (!is.null(ylim)) bp <- bp + ggplot2::lims(y=ylim)
  bp
}
