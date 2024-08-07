% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\name{qgcomp.glm.boot}
\alias{qgcomp.glm.boot}
\alias{qgcomp.boot}
\title{Quantile g-computation for continuous and binary outcomes}
\usage{
qgcomp.glm.boot(
  f,
  data,
  expnms = NULL,
  q = 4,
  breaks = NULL,
  id = NULL,
  weights,
  alpha = 0.05,
  B = 200,
  rr = TRUE,
  degree = 1,
  seed = NULL,
  bayes = FALSE,
  MCsize = nrow(data),
  parallel = FALSE,
  parplan = FALSE,
  ...
)

qgcomp.boot(
  f,
  data,
  expnms = NULL,
  q = 4,
  breaks = NULL,
  id = NULL,
  weights,
  alpha = 0.05,
  B = 200,
  rr = TRUE,
  degree = 1,
  seed = NULL,
  bayes = FALSE,
  MCsize = nrow(data),
  parallel = FALSE,
  parplan = FALSE,
  ...
)
}
\arguments{
\item{f}{R style formula}

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
id/cluster). Note that qgcomp.glm.noboot will not produce cluster-appropriate
standard errors. qgcomp.glm.boot can be used for this, which will use bootstrap
sampling of clusters/individuals to estimate cluster-appropriate standard
errors via bootstrapping.}

\item{weights}{"case weights" - passed to the "weight" argument of
\code{\link[stats]{glm}} or \code{\link[arm]{bayesglm}}}

\item{alpha}{alpha level for confidence limit calculation}

\item{B}{integer: number of bootstrap iterations (this should typically be >=200,
though it is set lower in examples to improve run-time).}

\item{rr}{logical: if using binary outcome and rr=TRUE, qgcomp.glm.boot will
estimate risk ratio rather than odds ratio}

\item{degree}{polynomial bases for marginal model (e.g. degree = 2
allows that the relationship between the whole exposure mixture and the outcome
is quadratic (default = 1).}

\item{seed}{integer or NULL: random number seed for replicable bootstrap results}

\item{bayes}{use underlying Bayesian model (\code{arm} package defaults). Results
in penalized parameter estimation that can help with very highly correlated
exposures. Note: this does not lead to fully Bayesian inference in general,
so results should be interpreted as frequentist.}

\item{MCsize}{integer: sample size for simulation to approximate marginal
zero inflated model parameters. This can be left small for testing, but should be as large
as needed to reduce simulation error to an acceptable magnitude (can compare psi coefficients for
linear fits with qgcomp.glm.noboot to gain some intuition for the level of expected simulation
error at a given value of MCsize). This likely won't matter much in linear models, but may
be important with binary or count outcomes.}

\item{parallel}{use (safe) parallel processing from the future and future.apply packages}

\item{parplan}{(logical, default=FALSE) automatically set future::plan to plan(multisession) (and set to existing plan, if any, after bootstrapping)}

\item{...}{arguments to glm (e.g. family)}
}
\value{
a qgcompfit object, which contains information about the effect
measure of interest (psi) and associated variance (var.psi), as well
as information on the model fit (fit) and information on the
marginal structural model (msmfit) used to estimate the final effect
estimates.
}
\description{
This function estimates a dose-response parameter representing a one quantile
increase in a set of exposures of interest. This model estimates the parameters of a marginal
structural model (MSM) based on g-computation with quantized exposures. Note: this function
allows non-linear and non-additive effects of individual components of the exposure, as well as
non-linear joint effects of the mixture via polynomial basis functions, which increase the
computational computational burden due to the need for non-parametric bootstrapping.
qgcomp.boot is an equivalent function (slated for deprecation)
}
\details{
Estimates correspond to the average expected change in the
(log) outcome per quantile increase in the joint exposure to all exposures
in `expnms'. Test statistics and confidence intervals are based on
a non-parametric bootstrap, using the standard deviation of the bootstrap
estimates to estimate the standard error. The bootstrap standard error is
then used to estimate Wald-type confidence intervals. Note that no bootstrapping
is done on estimated quantiles of exposure, so these are treated as fixed
quantities
}
\examples{
set.seed(30)
# continuous outcome
dat <- data.frame(y=rnorm(100), x1=runif(100), x2=runif(100), z=runif(100))
# Conditional linear slope
qgcomp.glm.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4, family=gaussian())
# Marginal linear slope (population average slope, for a purely linear,
#  additive model this will equal the conditional)
 \dontrun{
qgcomp.glm.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4,
  family=gaussian(), B=200) # B should be at least 200 in actual examples
# no intercept model
qgcomp.glm.boot(f=y ~ -1+z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4,
  family=gaussian(), B=200) # B should be at least 200 in actual examples

 # Note that these give different answers! In the first, the estimate is conditional on Z,
 # but in the second, Z is marginalized over via standardization. The estimates
 # can be made approximately the same by centering Z (for linear models), but
 # the conditional estimate will typically have lower standard errors.
 dat$z = dat$z - mean(dat$z)

# Conditional linear slope
qgcomp.glm.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4, family=gaussian())
# Marginal linear slope (population average slope, for a purely linear,
#  additive model this will equal the conditional)

qgcomp.glm.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4,
  family=gaussian(), B=200) # B should be at least 200 in actual examples

# Population average mixture slope which accounts for non-linearity and interactions
qgcomp.glm.boot(y ~ z + x1 + x2 + I(x1^2) + I(x2*x1), family="gaussian",
 expnms = c('x1', 'x2'), data=dat, q=4, B=200)

# generally non-linear/non-addiive underlying models lead to non-linear mixture slopes
qgcomp.glm.boot(y ~ z + x1 + x2 + I(x1^2) + I(x2*x1), family="gaussian",
 expnms = c('x1', 'x2'), data=dat, q=4, B=200, deg=2)

# binary outcome
dat <- data.frame(y=rbinom(50,1,0.5), x1=runif(50), x2=runif(50), z=runif(50))

# Conditional mixture OR
qgcomp.glm.noboot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'),
  data=dat, q=2)

#Marginal mixture OR (population average OR - in general, this will not equal the
# conditional mixture OR due to non-collapsibility of the OR)
qgcomp.glm.boot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'),
  data=dat, q=2, B=3, rr=FALSE)

# Population average mixture RR
qgcomp.glm.boot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'),
  data=dat, q=2, rr=TRUE, B=3)

# Population average mixture RR, indicator variable representation of x2
# note that I(x==...) operates on the quantile-based category of x,
# rather than the raw value
res = qgcomp.glm.boot(y ~ z + x1 + I(x2==1) + I(x2==2) + I(x2==3),
  family="binomial", expnms = c('x1', 'x2'), data=dat, q=4, rr=TRUE, B=200)
res$fit
plot(res)

# now add in a non-linear MSM
res2 = qgcomp.glm.boot(y ~ z + x1 + I(x2==1) + I(x2==2) + I(x2==3),
  family="binomial", expnms = c('x1', 'x2'), data=dat, q=4, rr=TRUE, B=200,
  degree=2)
res2$fit
res2$msmfit  # correct point estimates, incorrect standard errors
res2  # correct point estimates, correct standard errors
plot(res2)
# Log risk ratio per one IQR change in all exposures (not on quantile basis)
dat$x1iqr <- dat$x1/with(dat, diff(quantile(x1, c(.25, .75))))
dat$x2iqr <- dat$x2/with(dat, diff(quantile(x2, c(.25, .75))))
# note that I(x>...) now operates on the untransformed value of x,
# rather than the quantized value
res2 = qgcomp.glm.boot(y ~ z + x1iqr + I(x2iqr>0.1) + I(x2>0.4) + I(x2>0.9),
  family="binomial", expnms = c('x1iqr', 'x2iqr'), data=dat, q=NULL, rr=TRUE, B=200,
  degree=2)
res2
# using parallel processing

qgcomp.glm.boot(y ~ z + x1iqr + I(x2iqr>0.1) + I(x2>0.4) + I(x2>0.9),
  family="binomial", expnms = c('x1iqr', 'x2iqr'), data=dat, q=NULL, rr=TRUE, B=200,
  degree=2, parallel=TRUE, parplan=TRUE)


# weighted model
N=5000
dat4 <- data.frame(id=seq_len(N), x1=runif(N), x2=runif(N), z=runif(N))
dat4$y <- with(dat4, rnorm(N, x1*z + z, 1))
dat4$w=runif(N) + dat4$z*5
qdata = quantize(dat4, expnms = c("x1", "x2"), q=4)$data
# first equivalent models with no covariates
qgcomp.glm.noboot(f=y ~ x1 + x2, expnms = c('x1', 'x2'), data=dat4, q=4, family=gaussian())
qgcomp.glm.noboot(f=y ~ x1 + x2, expnms = c('x1', 'x2'), data=dat4, q=4, family=gaussian(),
              weights=w)

set.seed(13)
qgcomp.glm.boot(f=y ~ x1 + x2, expnms = c('x1', 'x2'), data=dat4, q=4, family=gaussian(),
            weights=w)
# using the correct model
set.seed(13)
qgcomp.glm.boot(f=y ~ x1*z + x2, expnms = c('x1', 'x2'), data=dat4, q=4, family=gaussian(),
            weights=w, id="id")
(qgcfit <- qgcomp.glm.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat4, q=4,
                       family=gaussian(), weights=w))
qgcfit$fit
summary(glm(y ~ z + x1 + x2, data = qdata, weights=w))
}
}
\seealso{
Other qgcomp_methods: 
\code{\link{qgcomp.cch.noboot}()},
\code{\link{qgcomp.cox.boot}()},
\code{\link{qgcomp.cox.noboot}()},
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
