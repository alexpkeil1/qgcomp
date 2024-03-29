% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_zi.R
\name{zimsm_fit.control}
\alias{zimsm_fit.control}
\title{Control of fitting parameters for zero inflated MSMs}
\usage{
zimsm_fit.control(predmethod = c("components", "catprobs"))
}
\arguments{
\item{predmethod}{character in c("components", "catprobs"). "components" simulates from the
model parameters directly while "catprobs" simulates outcomes from the category specific
probabilities, which is output from predict.zeroinfl. The former is slightly
more flexible and stable, but the latter is preferred in zero inflated negative bionomial models.}
}
\description{
this is an internal function called by
\code{\link[qgcomp]{qgcomp.zi.boot}}, but is documented here for clarity.
Generally, users will not need to call this function directly.
}
\details{
Provides fine control over zero inflated MSM fitting
}
