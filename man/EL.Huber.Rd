\name{EL.Huber}
\alias{EL.Huber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Empirical likelihood test for the difference of smoothed Huber estimators
}
\description{
Empirical likelihood inference for the difference of smoothed Huber estimators. 
This includes a test for the null hypothesis for a constant difference of smoothed 
Huber estimators, confidence interval and EL estimator.
}
\usage{
EL.Huber(X, Y, mu = 0, conf.level = 0.95, 
         scaleX=1, scaleY=1, VX = 2.046, VY = 2.046, k = 1.35)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a vector of data values.
}
  \item{Y}{
a vector of data values.
}
  \item{mu}{
a number specifying the null hypothesis.
}
  \item{conf.level}{
confidence level of the interval.
}
  \item{scaleX}{
the scale estimate of sample 'X'.
}
  \item{scaleY}{
the scale estimate of sample 'Y'.
}
  \item{VX}{
the asymptotic variance of initial (nonsmooth) Huber estimator for the sample 'X'.
}
  \item{VY}{
the asymptotic variance of initial (nonsmooth) Huber estimator for the sample 'Y'.
}
  \item{k}{
tuning parameter for the Huber estimator.
}
}
\value{
A list of class 'htest' containing the following components:
\item{estimate }{the empirical likelihood estimate for the difference of two smoothed Huber estimators.}
\item{conf.int }{a confidence interval for the difference of two smoothed Huber estimators.}
\item{p.value }{the p-value for the test.}
\item{statistic }{the value of the test statistic.}
\item{method }{the character string 'Empirical likelihood smoothed Huber estimator difference test'.}
\item{null.value }{the specified hypothesized value of the mean difference 'mu' under the null hypothesis.}
\item{data.name }{a character string giving the names of the data.}
}
\details{
A common choice for a robust scale estimate (parameters scaleX and scaleY) is 
the mean absolute deviation (MAD).
}
\references{
J. Valeinis, E. Cers. Extending the two-sample empirical likelihood. To be published.
Preprint available at http://home.lanet.lv/~valeinis/lv/petnieciba/EL_TwoSample_2011.pdf.

F. Hampel, C. Hennig and E. A. Ronchetti (2011).  A smoothing principle for the Huber and other location M-estimators, Computational Statistics & Data Analysis, 55(1), 324-337.
}
\author{
E. Cers, J. Valeinis
}

\seealso{
\code{\link{EL.means}}
}
\examples{
X <- rnorm(100)
Y <- rnorm(100)
t.test(X, Y)
EL.means(X, Y)
EL.Huber(X, Y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~nonparametric }
\keyword{ ~smooth }
\keyword{ ~htest }
