% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_generics.R
\name{print.qgcompfit}
\alias{print.qgcompfit}
\title{Default printing method for a qgcompfit object}
\usage{
\method{print}{qgcompfit}(x, showweights = TRUE, ...)
}
\arguments{
\item{x}{"qgcompfit" object from \code{qgcomp}, \code{qgcomp.noboot} or \code{qgcomp.boot}
function}

\item{showweights}{logical: should weights be printed, if estimated?}

\item{...}{unused}
}
\description{
Gives variable output depending on whether \code{qgcomp.noboot} or \code{qgcomp.boot}
is called. For \code{qgcomp.noboot} will output final estimate of joint exposure
effect (similar to the 'index' effect in weighted quantile sums), as well
as estimates of the 'weights' (standardized coefficients). For \code{qgcomp.boot},
the marginal effect is given, but no weights are reported since this approach
generally incorporates non-linear models with interaction terms among exposures,
which preclude weights with any useful interpretation.
}
\examples{
set.seed(50)
dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
obj1 <- qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
obj2 <- qgcomp.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=10, seed=125)
# does not need to be explicitly called, but included here for clarity
print(obj1)
print(obj2)
}
\seealso{
\code{\link[qgcomp]{qgcomp.noboot}}, \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}}
}
\concept{variance mixtures}
