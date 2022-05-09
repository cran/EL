\name{EL.smooth}
\alias{EL.smooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Smooth estimates and confidence intervals (or simultaneous bands) using the smoothed two-sample EL method
}
\description{
Calculates estimates and pointwise confidence intervals (or simultaneous bands) for P-P and Q-Q plots, ROC curves, quantile differences (qdiff) and CDF differences
(ddiff) using the smoothed empirical likelihood method.
}
\usage{
EL.smooth(method, X, Y, t, bw = bw.nrd0,
          conf.level = NULL, simultaneous = FALSE,
          bootstrap.samples = 300, more.warnings = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
 "pp", "qq", "roc", "qdiff" or "ddiff".
}
  \item{X}{
a vector of data values.
}
  \item{Y}{
a vector of data values.
}
  \item{t}{
a vector of points for which to calculate the estimates and confidence
intervals.
}
  \item{conf.level}{
confidence level for the intervals. A number between 0 and 1 or NULL when no
confidence bands should be calculated. Depending on the value of 'simultaneous'
either pointwise intervals or simultaneous confidence bands will be calculated.
}
  \item{simultaneous}{
if this is TRUE, simultaneous confidence bands will be constructed, using a
nonparametric bootstrap procedure to select the level of confidence bands. The
default is FALSE, in which case simple pointwise confidence bands are calculated.
}
  \item{bootstrap.samples}{
the number of samples used to bootstrap the simultaneous confidence bands
when 'simultaneous = TRUE'.
}
  \item{bw}{
a function taking a vector of values and returning the corresponding
bandwidth or a vector of two values corresponding to the respective
bandwidths of X and Y.
}
  \item{more.warnings}{
if this is FALSE (the default) a single warning will be produced if
there is any problem calculating the estimate or the confidence
bands. If this is set to TRUE a warning will be produced for every
point at which there was a problem.
}
}
\details{
Confidence bands are drawn only if 'conf.level' is not 'NULL'.

  When constructing simultaneous confidence bands, it is advisable to
check whether the chosen range of 't' values does not produce too large
bands (for example, for the P-P plot in the example below the
interval [0.05, 0.95] was a sensible choice). This has to be checked for
each data sample separately by hand.  Note that the calculation of
simultaneous confidence bands can take a long time.
}
\value{
\item{estimate }{ the estimated values at points 't'.}
\item{conf.int }{ a two column matrix where each row represents the
lower and upper bounds of the confidence bands corresponding to the
values at points 't'.}
\item{simultaneous.conf.int}{ will be a true value if simultaneous
confidence bands are constructed.
}
\item{bootstrap.crit}{ the critical value from the bootstrapped -2 *
  log-likelihood statistic for simultaneous confidence
  bands using the confidence level 'conf.level'. Only calculated when
  'conf.level' is not NULL and 'simultaneous' is TRUE.
}
}
\references{
J. Valeinis and E. Cers. Extending the two-sample empirical likelihood. To be published.
Preprint available at http://home.lanet.lv/~valeinis/lv/petnieciba/EL_TwoSample_2011.pdf

P. Hall and A. Owen (1993). Empirical likelihood bands in density estimation. Journal of
Computational and Graphical statistics, 2(3), 273-289.
}
\author{
E. Cers, J. Valeinis
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{EL.plot}}
\code{\link{EL.statistic}}
}
\examples{
\donttest{
#### Simultaneous confidence bands for a P-P plot
X1 <- rnorm(200)
X2 <- rnorm(200, 1)

x <- seq(0.05, 0.95, length=19)
y <- EL.smooth("pp", X1, X2, x, conf.level=0.95,
               simultaneous=TRUE, bw=c(0.3, 0.3))

## Plot the graph with both pointwise and simultaneous confidence bands
EL.plot("pp", X1, X2, conf.level=0.95, bw=c(0.3, 0.3))
lines(x, y$conf.int[1,], lty="dotted")
lines(x, y$conf.int[2,], lty="dotted")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~hplot }
\keyword{ ~nonparametric }
\keyword{ ~smooth }
