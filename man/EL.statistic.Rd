\name{EL.statistic}
\alias{EL.statistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The two-sample empirical likelihood statistic}

\description{ 
Calculates -2 times the log-likelihood ratio statistic when the function
of interest (either of P-P or Q-Q plot, ROC curve, difference of quantile or
distribution functions) at some point 't' is equal to 'd'.
}
\usage{
EL.statistic(method, X, Y, d, t, bw = bw.nrd0)
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
  \item{d}{
a number
}
  \item{t}{
a number.
}
  \item{bw}{
a function taking a vector of values and returning the corresponding
bandwidth or a vector of two values corresponding to the respective
bandwidths of X and Y.
}
}
\value{
-2 times the logarithm of the two-sample empirical likelihood ratio.
}
\references{
J.Valeinis, E.Cers. Extending the two-sample empirical likelihood. To be published.
Preprint available at http://home.lanet.lv/%7Evaleinis/lv/petnieciba/EL_TwoSample_2011.pdf
}
\author{
E. Cers, J. Valeinis
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{EL.smooth}}
}
\examples{

EL.statistic("pp", rnorm(100), rnorm(100), 0.5, 0.5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~nonparametric }
\keyword{ ~smooth }% __ONLY ONE__ keyword per line
\keyword{ ~htest }
