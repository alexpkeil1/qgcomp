% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_experimental.R
\name{qgcomp.cch.noboot}
\alias{qgcomp.cch.noboot}
\title{Quantile g-computation for survival outcomes in a case-cohort design under linearity/additivity}
\usage{
qgcomp.cch.noboot(
  f,
  data,
  subcoh = NULL,
  id = NULL,
  cohort.size = NULL,
  expnms = NULL,
  q = 4,
  breaks = NULL,
  weights,
  cluster = NULL,
  alpha = 0.05,
  ...
)
}
\arguments{
\item{f}{R style survival formula, which includes \code{\link[survival]{Surv}}
in the outcome definition. E.g. \code{Surv(time,event) ~ exposure}. Offset
terms can be included via \code{Surv(time,event) ~ exposure + offset(z)}}

\item{data}{data frame}

\item{subcoh}{(From \code{\link[survival]{cch}} help) Vector of indicators for subjects sampled as part of the sub-cohort. Code 1 or TRUE for members of the sub-cohort, 0 or FALSE for others. If data is a data frame then subcoh may be a one-sided formula.}

\item{id}{(From \code{\link[survival]{cch}} help) Vector of unique identifiers, or formula specifying such a vector.}

\item{cohort.size}{(From \code{\link[survival]{cch}} help) Vector with size of each stratum original cohort from which subcohort was sampled}

\item{expnms}{character vector of exposures of interest}

\item{q}{NULL or number of quantiles used to create quantile indicator variables
representing the exposure variables. If NULL, then gcomp proceeds with un-transformed
version of exposures in the input datasets (useful if data are already transformed,
or for performing standard g-computation)}

\item{breaks}{(optional) NULL, or a list of (equal length) numeric vectors that
characterize the minimum value of each category for which to
break up the variables named in expnms. This is an alternative to using 'q'
to define cutpoints.}

\item{weights}{Not used here}

\item{cluster}{not yet implemented}

\item{alpha}{alpha level for confidence limit calculation}

\item{...}{arguments to glm (e.g. family)}
}
\value{
a qgcompfit object, which contains information about the effect
measure of interest (psi) and associated variance (var.psi), as well
as information on the model fit (fit) and information on the
weights/standardized coefficients in the positive (pos.weights) and
negative (neg.weights) directions.
}
\description{
This function performs quantile g-computation in a survival
setting. The approach estimates the covariate-conditional hazard ratio for
a joint change of 1 quantile in each exposure variable specified in expnms
parameter
}
\details{
For survival outcomes (as specified using methods from the
survival package), this yields a conditional log hazard ratio representing
a change in the expected conditional hazard (conditional on covariates)
from increasing every exposure by 1 quantile. In general, this quantity
quantity is not equivalent to g-computation estimates. Hypothesis test
statistics and 95\% confidence intervals are based on using the delta
estimate variance of a linear combination of random variables.
}
\examples{
set.seed(50)
N=200
dat <- data.frame(time=(tmg <- pmin(.1,rweibull(N, 10, 0.1))), 
                d=1.0*(tmg<0.1), x1=runif(N), x2=runif(N), z=runif(N))
expnms=paste0("x", 1:2)
f = survival::Surv(time, d)~x1 + x2
(fit1 <- survival::coxph(f, data = dat))
(obj <- qgcomp.cox.noboot(f, expnms = expnms, data = dat))
\dontrun{

# weighted analysis
dat$w = runif(N)
qdata = quantize(dat, expnms=expnms)
(obj2 <- qgcomp.cox.noboot(f, expnms = expnms, data = dat, weight=w))
obj2$fit
survival::coxph(f, data = qdata$data, weight=w)

# not run: bootstrapped version is much slower
(obj2 <- qgcomp.cox.boot(f, expnms = expnms, data = dat, B=200, MCsize=20000))
}
}
\seealso{
Other qgcomp_methods: 
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
\code{\link{qgcomp.zi.boot}()},
\code{\link{qgcomp.zi.noboot}()}
}
\concept{qgcomp_methods}
\concept{variance mixtures}
