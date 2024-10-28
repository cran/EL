\name{EL.plot}
\alias{EL.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draws plots using the smoothed two-sample empirical likelihood method
}
\description{
Draws P-P and Q-Q plots, ROC curves, quantile differences (qdiff) and CDF
differences (ddiff) and their respective confidence bands (pointwise or simultaneous) using the empirical
likelihood method.
}
\usage{
EL.plot(method, X, Y, bw = bw.nrd0, conf.level = NULL,
        simultaneous = FALSE, bootstrap.samples = 300,
        more.warnings = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
 "pp", "qq", "roc", "qdiff" or "fdiff".
}
  \item{X}{
a vector of data values.
}
  \item{Y}{
a vector of data values.
}
  \item{bw}{
a function taking a vector of values and returning the corresponding
bandwidth or a vector of two values corresponding to the respective
bandwidths of X and Y.
}
  \item{conf.level}{
confidence level for the intervals. A number between 0 and 1 or NULL when
no confidence bands should be calculated. Depending on the value of
'simultaneous' either pointwise intervals or simultaneous confidence bands
will be drawn.
}
  \item{simultaneous}{
if this is TRUE, simultaneous confidence bands will be constructed, using
a nonparametric bootstrap procedure to select the level of confidence
bands. The default is FALSE, in which case simple pointwise confidence bands are
calculated.
}
  \item{bootstrap.samples}{
the number of samples used to bootstrap the simultaneous confidence bands
when 'simultaneous = TRUE'.
}
  \item{more.warnings}{
if this is FALSE (the default) a single warning will be produced if
there is any problem calculating the estimate or the confidence
bands. If this is set to TRUE a warning will be produced for every
point at which there was a problem.
}
  \item{\dots}{
further arguments passed to plot.
}
}
\details{
  The plotting interval for P-P plots, ROC curves and differences of
  quantile functions is [0, 1] (where these functions are defined). The
  Q-Q plot is drawn from the minimum to the maximum of 'Y'. Finally, for
  the plot of distribution function differences the interval from max(min(X),
  min(Y)) to min(max(X), max(Y)) is used.

  Confidence bands are drawn only if 'conf.level' is not 'NULL'.

When constructing simultaneous confidence bands, the plot is drawn on an
interval that is narrowed by 5\% on both sides, since the procedure is
usually sensitive at the end-points, which can result in large bands.  The
confidence level for the simultaneous confidence bands is bootstrapped
using 50 evenly spaced points in this interval. If the default interval
produces too large confidence bands, use the function 'EL.smooth' where
the intervals are specified manually. Note that calculation of
simultaneous confidence bands can take a long time.

 }

\value{
none.
}
\references{
J. Valeinis, E. Cers. Extending the two-sample empirical likelihood. To be published.
Preprint available at http://home.lanet.lv/~valeinis/lv/petnieciba/EL_TwoSample_2011.pdf.

P. Hall and A. Owen (1993). Empirical likelihood bands in density estimation.
Journal of Computational and Graphical statistics, 2(3), 273-289.
}
\author{
E. Cers, J. Valeinis
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{EL.smooth}}
\code{\link{EL.statistic}}
}
\examples{
## The examples showcase all available graphs

X1 <- rchisq(100, 2.5)
X2 <- rnorm(100, 0, 1)

\donttest{
X1 <- rchisq(100, 2.5)
X2 <- rnorm(100, 0, 1)

# Intro
xlim <- c(min(X1, X2) - 0.5, max(X1, X2) + 0.5)
D1 <- density(X1)
D2 <- density(X2)
ylim <- c(min(D1$y, D2$y), max(D1$y, D2$y))
df <- data.frame(x1 = D1$x, y1 = D1$y, x2 = D2$x, y2 = D2$y)
p1 <- ggplot2::ggplot(data = df) +
    ggplot2::geom_line(ggplot2::aes(x=x2, y=y2, color=paste0('X2 (bw=', round(D2$bw, 2), ')'))) +
    ggplot2::geom_line(ggplot2::aes(x=x1, y=y1, color=paste0('X1 (bw=', round(D1$bw, 2), ')'))) +
    ggplot2::guides(color = ggplot2::guide_legend(title = NULL)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::labs(x="X", y="Density")
p1

# CDF differences
p2 <- EL.plot("fdiff", X1, X2, main="F difference", conf.level=0.95)
tt <- seq(max(c(min(X1), min(X2))), min(c(max(X1), max(X2))), length=30)
ee <- ecdf(X2)(tt) - ecdf(X1)(tt)
p2 <- p2 + ggplot2::geom_point(data=data.frame(tt = tt, ee = ee), ggplot2::aes(x=tt, y=ee))
p2

# Quantile differences
p3 <- EL.plot("qdiff", X1, X2, main="Quantile difference", conf.level = 0.95)
tt <- seq(0.01, 0.99, length=30)
ee <- quantile(X2, tt) - quantile(X1, tt)
p3 <- p3 + ggplot2::geom_point(data=data.frame(tt = tt, ee = ee), ggplot2::aes(x=tt, y=ee))
p3

# Q-Q plot
p4 <- EL.plot("qq", X1, X2, main="Q-Q plot", conf.level=0.95)
tt <- seq(min(X2), max(X2), length=30)
ee <- quantile(X1, ecdf(X2)(tt))
p4 <- p4 + ggplot2::geom_point(data=data.frame(tt = tt, ee = ee), ggplot2::aes(x=tt, y=ee))
p4

# P-P plot
p5 <- EL.plot("pp", X1, X2, main="P-P plot", conf.level=0.95, ylim=c(0,1))
tt <- seq(0.01, 0.99, length=30)
ee <- ecdf(X1)(quantile(X2, tt))
p5 <- p5 + ggplot2::geom_point(data=data.frame(tt = tt, ee = ee), ggplot2::aes(x=tt, y=ee))
p5

# ROC curve
p6 <- EL.plot("roc", X1, X2, main="ROC curve", conf.level=0.95, ylim=c(0,1))
tt <- seq(0.01, 0.99, length=30)
ee <- 1- ecdf(X1)(quantile(X2, 1-tt))
p6 <- p6 + ggplot2::geom_point(data=data.frame(tt = tt, ee = ee), ggplot2::aes(x=tt, y=ee))
p6

# Showing all plots at once is outside of scope from
# these examples but to do so run the following:
# require(cowplot)
# cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~nonparametric }
\keyword{ ~smooth }
\keyword{ ~hplot }
