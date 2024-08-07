% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_zi.R
\name{qgcomp.zi.noboot}
\alias{qgcomp.zi.noboot}
\title{Quantile g-computation for zero-inflated count outcomes under linearity/additivity}
\usage{
qgcomp.zi.noboot(
  f,
  data,
  expnms = NULL,
  q = 4,
  breaks = NULL,
  id = NULL,
  weights,
  alpha = 0.05,
  bayes = FALSE,
  ...
)
}
\arguments{
\item{f}{R style formula using syntax from 'pscl' package: depvar ~ indvars_count | indvars_zero}

\item{data}{data frame}

\item{expnms}{character vector of exposures of interest}

\item{q}{NULL or number of quantiles used to create quantile indicator variables
representing the exposure variables. If NULL, then gcomp proceeds with un-transformed
version of exposures in the input datasets (useful if data are already transformed,
or for performing standard g-computation)}

\item{breaks}{(optional) NULL, or a list of (equal length) numeric vectors that
characterize the minimum value of each category for which to
break up the variables named in expnms. This is an alternative to using 'q'
to define cutpoints.}

\item{id}{(optional) NULL, or variable name indexing individual units of
observation (only needed if analyzing data with multiple observations per
id/cluster)}

\item{weights}{"case weights" - passed to the "weight" argument of
\code{\link[pscl]{zeroinfl}}.}

\item{alpha}{alpha level for confidence limit calculation}

\item{bayes}{not yet implemented}

\item{...}{arguments to zeroinfl (e.g. dist)}
}
\value{
a qgcompfit object, which contains information about the effect
measure of interest (psi) and associated variance (var.psi), as well
as information on the model fit (fit) and information on the
weights/standardized coefficients in the positive (pos.weights) and
negative (neg.weights) directions.
}
\description{
This function estimates a linear dose-response parameter representing a one quantile
increase in a set of exposures of interest for zero-inflated count outcomes. This function is
limited to linear and additive
effects of individual components of the exposure. This model estimates the parameters of a marginal
structural zero-inflated model (MSM) based on g-computation with quantized exposures.
Note: this function is valid only under linear and additive effects of individual components of the exposure, but when
these hold the model can be fit with very little computational burden.
}
\details{
A zero-inflated version of quantile g-computation based on the implementation in the
'pscl' package. A zero-inflated distribution is a mixture distribution in which one of the
distributions is a point mass at zero (with probability given by a logistic model), and the
other distribution is a discrete or continuous distribution.
This estimates the effect of a joint increase in all exposures on 1) the odds
of belonging to the "zero" vs. "count" portions of the distribution and/or 2) the rate parameter
for the "count" portion of the distribution.
}
\examples{
set.seed(50)
n=100
dat <- data.frame(y=rbinom(n, 1, 0.5)*rpois(n, 1.2), x1=runif(n), x2=runif(n), z=runif(n))

# poisson count model, mixture in both portions
qgcomp.zi.noboot(f=y ~ z + x1 + x2 | x1 + x2, expnms = c('x1', 'x2'), 
    data=dat, q=2, dist="poisson")
    
# negative binomial count model, mixture and covariate in both portions
qgcomp.zi.noboot(f=y ~ z + x1 + x2 | z + x1 + x2, expnms = c('x1', 'x2'), 
   data=dat, q=2, dist="negbin")  
qgcomp.zi.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), 
   data=dat, q=2, dist="negbin") # equivalent
   
# negative binomial count model, mixture only in the 'count' portion of the model
qgcomp.zi.noboot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=2, dist="negbin")

# weighted analysis
dat$w = runif(n)*5
qgcomp.zi.noboot(f=y ~ z + x1 + x2 | x1 + x2, expnms = c('x1', 'x2'), 
    data=dat, q=2, dist="poisson", weights=w)
# Expect this:     
# Warning message:
# In eval(family$initialize) : non-integer #successes in a binomial glm!

}
\seealso{
Other qgcomp_methods: 
\code{\link{qgcomp.cch.noboot}()},
\code{\link{qgcomp.cox.boot}()},
\code{\link{qgcomp.cox.noboot}()},
\code{\link{qgcomp.glm.boot}()},
\code{\link{qgcomp.glm.ee}()},
\code{\link{qgcomp.glm.noboot}()},
\code{\link{qgcomp.hurdle.boot}()},
\code{\link{qgcomp.hurdle.noboot}()},
\code{\link{qgcomp.multinomial.boot}()},
\code{\link{qgcomp.multinomial.noboot}()},
\code{\link{qgcomp.partials}()},
\code{\link{qgcomp.zi.boot}()}
}
\concept{qgcomp_methods}
\concept{variance mixtures}
