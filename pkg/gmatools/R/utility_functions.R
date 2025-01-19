## Small convenience utilities exported by the package

##' Signed logarithmic transformation that is smooth at 0
##'
##' This function applies a signed logarithmic transformation to the argument
##' using the mapping \deqn{x \mapsto \mathop{\text{sgn}}(x) \cdot \log_b(1 + |x|)}{x -> sgn(x) * log_b(1 + |x|)}
##' 
##' @param x a numeric vector or matrix containing the values to be transformed
##' @param base base \eqn{b} of the logarithm to apply (defaults to \eqn{e})
##' @returns A numeric vector or matrix of the same shape with the logarithmic transformation applied.
##' 
##' @examples
##' signed.log(seq(-3, 3), base=2)
##' \dontrun{
##' plot(signed.log, from=-4, to=4, xlim=c(-4, 4), ylim=c(-4, 4),
##'      lwd=3, col=2, xaxs="i", yaxs="i")
##' abline(h=0, v=0, lwd=2); abline(0, 1, col=4)}
##' @export
signed.log <- function (x, base=exp(1)) {
  sign(x) * log(1 + abs(x), base=base)
}

##' Expand range by a specified proportion on either side
##' 
##' This function computes the range of the specified numerical values,
##' then expands it on either side by a specified proportion of the full range.
##' Its main use is to pick a suitable range for scatterplots, ensuring a margin around all data points.
##' 
##' @param x a numeric vector or matrix whose range is to be extended
##' @param by proportion by which to expand the range on either side (default: 20%)
##' @returns A numeric vector of length 2 with the expanded range.
##' 
##' @examples
##' expand.range(1:10, by=.1)
##' @export
expand.range <- function(x, by=.2) {
  rg <- range(x)
  rg + diff(rg) * by * c(-1, 1)
}

##' Compute Cohen's d as scale-adjusted effect size of a t-test
##'
##' This function computes Cohen's \eqn{d} as a scale-adjusted measure of effect size
##' for a t-test of two independent samples `x` and `y`, 
##' according to https://en.wikipedia.org/wiki/Effect_size#Cohen's_d
##' 
##' @param x,y two numeric vectors giving the data points to be compared
##' @param conf.level if specified, return confidence interval of \eqn{d} at the specified
##'        confidence level in range \eqn{[0, 1]} rather than point estimate
##' @returns Point estimate or confidence interval for Cohen's \eqn{d}.
##' 
##' @details
##' Effect size is based on the difference \eqn{\mu_x - \mu_y} in order to be compatible with
##' [t.test()]. Positive effect size thus indicates that `x` values are larger than `y` values on average.
##' @export
##' @importFrom stats t.test var
cohen.d <- function (x, y, conf.level=NULL) {
  nx <- length(x)
  ny <- length(y)
  s2 <- ((nx-1) * var(x) + (ny-1) * var(y)) / (nx + ny - 2) # within-group variance
  if (is.null(conf.level)) {
    (mean(x) - mean(y)) / sqrt(s2)
  }
  else {
    t.test(x, y, conf.level=conf.level)$conf.int / sqrt(s2)
  }
}
