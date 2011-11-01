\name{EL.means}
\alias{EL.means}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Empirical likelihood test for the difference of two sample means
}
\description{
Empirical likelihood inference for the difference of two sample means. 
This includes a test for the null hypothesis for a constant difference 
of mean difference, confidence interval and EL estimator.
}
\usage{
EL.means(X, Y, mu = 0, conf.level = 0.95)
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
}
\value{
A list of class 'htest' containing the following components:
\item{estimate }{the empirical likelihood estimate of the mean difference.}
\item{conf.int }{a confidence interval for the mean difference.}
\item{p.value }{the p-value for the test.}
\item{statistic }{the value of the test statistic.}
\item{method }{the character string 'Empirical likelihood mean difference test'.}
\item{null.value }{the specified hypothesized value of mean differences 'mu' under the null hypothesis.}
\item{data.name }{a character string giving the names of the data.}
}
\references{
J. Valeinis, E. Cers and J. Cielens (2010). Two-sample problems in statistical data modelling. 
Mathematical modelling and analysis, 15(1), 137-151.

J. Valeinis, E. Cers. Extending the two-sample empirical likelihood. To be published.
Preprint available at http://home.lanet.lv/~valeinis/lv/petnieciba/EL_TwoSample_2011.pdf.
}
\author{
E. Cers, J. Valeinis
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{EL.Huber}}
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
\keyword{ ~htest }% __ONLY ONE__ keyword per line
