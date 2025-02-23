
#' @title Quantile g-computation for left-censored outcomes
#'  
#' @description This function yields population average effect estimates for 
#'   (possibly left censored)  outcomes
#'  #'
#' @param f an r formula representing the conditional model for the outcome, given all
#' exposures and covariates. Interaction terms that include exposure variables
#' should be represented via the \code{\link[base]{AsIs}} function
#' @param data data frame
#' @param expnms character vector of exposures of interest
#' @param q NULL or number of quantiles used to create quantile indicator variables
#' representing the exposure variables. If NULL, then gcomp proceeds with un-transformed
#' version of exposures in the input datasets (useful if data are already transformed,
#' or for performing standard g-computation)
#' @param breaks (optional) NULL, or a list of (equal length) numeric vectors that 
#' characterize the minimum value of each category for which to 
#' break up the variables named in expnms. This is an alternative to using 'q'
#' to define cutpoints.
#' @param id Not used
#' @param weights "case weights" - passed to the "weight" argument of AER::tobit. Note this must be a standalone vector so that it might accept, for example, an argument like data$weights but will not accept "weights" if there is a variable "weights" in the dataset representing the weights.
#' \code{\link[AER]{tobit}}
#' @param cluster Passed to AER::tobit (Optional variable that identifies groups of subjects, used in computing the robust variance. Like model variables, this is searched for in the dataset pointed to by the data argument). Note: this will result in robust standard errors being returned unless robust=FALSE is used in the function call.
#' @param alpha alpha level for confidence limit calculation
#' @param left Passed to AER::tobit (From tobit docs: left limit for the censored dependent variable y. If set to -Inf, y is assumed not to be left-censored.)
#' @param right Passed to AER::tobit (From tobit docs: right limit for the censored dependent variable y. If set to Inf, the default, y is assumed not to be right-censored.)
#' @param ... arguments to AER::tobit (e.g. dist), and consequently to survival::survreg
#'
#' @return a qgcompfit object, which contains information about the effect
#'  measure of interest (psi) and associated variance (var.psi), as well
#'  as information on the model fit (fit) and information on the
#'  weights/standardized coefficients in the positive (pos.weights) and
#'  negative (neg.weights) directions.
#' @export
#' @importFrom AER tobit
#'
#' @examples
#' # First, a demonstration that the tobit model will return identical results 
#' # (with variation due to numerical algorithms)
#' # in the case of no censored values
#' set.seed(50)
#' # linear model
#' dat <- data.frame(y=runif(50,-1,1), x1=runif(50), x2=runif(50), z=runif(50))
#' qgcomp.glm.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian())
#' qgcomp.tobit.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
#' # not intercept model
#' qgcomp.glm.noboot(f=y ~-1+ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, 
#'     q=2, family=gaussian())
#' qgcomp.tobit.noboot(f=y ~-1+ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, 
#'     q=2, dist="gaussian", left = -Inf, right = Inf)
#' # Next, a demonstration that the tobit model address censoring
#' uncens_dat = qgcomp:::.dgm_quantized(N=1000, b0=0, coef=c(.33,.33,.33,.0))
#' uncens_dat$ycens = ifelse(uncens_dat$y>0, uncens_dat$y, 0) # censor at zero
#' uncens_dat$ycens1 = ifelse(uncens_dat$y>0, uncens_dat$y, 1) # censor at higher value
#' # q set to NULL because the exposures are already quantized
#' # again, first check the uncensored outcome for agreement
#' qgcomp.glm.noboot(f= y ~ x1+x2+x3+x4, expnms = c('x1', 'x2', 'x3', 'x4'), 
#'     data=uncens_dat, q=NULL, family=gaussian())
#' qgcomp.tobit.noboot(f= y ~x1+x2+x3+x4, expnms = c('x1', 'x2', 'x3', 'x4'), 
#'     data=uncens_dat, q=NULL, dist="gaussian", left = -Inf, right = Inf)
#' # Next, after censoring
#' ft_std <- qgcomp.glm.noboot(f= ycens ~ x1+x2+x3+x4, 
#'     expnms = c('x1', 'x2', 'x3', 'x4'), data=uncens_dat, q=NULL, 
#'     family=gaussian())
#' ft_tobit <- qgcomp.tobit.noboot(f= ycens ~x1+x2+x3+x4, 
#'     expnms = c('x1', 'x2', 'x3', 'x4'), data=uncens_dat, q=NULL, 
#'     dist="gaussian", left = 0, right = Inf)
#' # note the tobit regression will be closer to the estimates given when 
#' # fitting with the uncensored outcome (which will typically not be available)
#' summary(ft_std) 
#' summary(ft_tobit)
#' ft_std1 <- qgcomp.glm.noboot(f= ycens1 ~ x1+x2+x3+x4, 
#'     expnms = c('x1', 'x2', 'x3', 'x4'), data=uncens_dat, q=NULL, 
#'     family=gaussian())
#' ft_tobit1 <- qgcomp.tobit.noboot(f= ycens1 ~x1+x2+x3+x4, 
#'     expnms = c('x1', 'x2', 'x3', 'x4'), data=uncens_dat, q=NULL, 
#'     dist="gaussian", left = 1, right = Inf)
#' # the bias from standard methods is more extreme at higher censoring levels
#' summary(ft_std1) 
#' summary(ft_tobit1)
qgcomp.tobit.noboot <- function (f, data, expnms = NULL, q = 4, breaks = NULL, id = NULL, 
                                  weights = NULL, cluster = NULL, alpha = 0.05, left = -Inf, right = Inf, ...) 
{
  #weights <- data[weights_name]
  newform <- terms(f, data = data)
  hasintercept = as.logical(attr(newform, "intercept"))
  class(newform) <- "formula"
  nobs = nrow(data)
  origcall <- thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("f", "data", "weights", "offset"), names(thecall), 
             0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, 
                      FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE
  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame())
  if (hasweights) {
    data$weights <- weights
  }
  else data$weights = rep(1, nobs)
  if (is.null(expnms)) {
    message("Including all model terms as exposures of interest")
    expnms <- attr(newform, "term.labels")
  }
  if (!is.null(q) | !is.null(breaks)) {
    ql <- quantize(data, expnms, q, breaks)
    qdata <- ql$data
    br <- ql$breaks
  }
  else {
    qdata <- data
    br <- breaks
  }
  
  environment(newform) <- list2env(list(qdata = qdata))
  
  fit <- AER::tobit(formula = newform, data = qdata, weights = weights, left = left,
                    right = right, ...)
  fit[["family"]] = tobit()
  
  mod <- summary(fit)
  fit$data <- qdata
  
  # ##
  estb <- sum(mod$coefficients[expnms, 1])
  covMat = fit$var
  class(covMat) <- "matrix"
  colnames(covMat) <-  rownames(covMat) <- names(mod$coefficients[, 1])
  seb <- se_comb(expnms, covmat = covMat)
  if(hasintercept){
    estb <- c(fit$coefficients[1], estb)
    seb <- c(sqrt(covMat[1,1,drop=FALSE]), seb)
  }
  psiidx = 1+hasintercept
  names(estb)[psiidx] <- "psi1"
  tstat <- estb/seb
  pvalz <- 2 - 2 * pnorm(abs(tstat))
  ci <- cbind(estb + seb * qnorm(alpha/2), estb + seb * qnorm(1 - alpha/2))
  wcoef <- fit$coefficients[expnms]
  names(wcoef) <- gsub("_q", "", names(wcoef))
  poscoef <- which(wcoef > 0)
  negcoef <- which(wcoef <= 0)
  pos.weights <- abs(wcoef[poscoef])/sum(abs(wcoef[poscoef]))
  neg.weights <- abs(wcoef[negcoef])/sum(abs(wcoef[negcoef]))
  pos.psi <- sum(wcoef[poscoef])
  neg.psi <- sum(wcoef[negcoef])
  qx <- qdata[, expnms]
  names(qx) <- paste0(names(qx), "_q")
  cnms = "psi1"
  if(hasintercept){
    covmat.coef=vc_comb(aname="(Intercept)", expnms=expnms, covmat = covMat)
    cnms = c("(intercept)", cnms)
  }
  if(!hasintercept)
    covmat.coef=as.matrix(seb^2,nrow=1,ncol=1)
  colnames(covmat.coef) <- rownames(covmat.coef) <- names(estb) <- cnms

  res <- .qgcomp_object(qx = qx, fit = fit, psi = estb, var.psi = seb^2,
                        covmat.psi = seb^2, ci = ci, coef = estb, var.coef = seb^2,
                        covmat.coef = seb^2, ci.coef = ci, expnms = expnms, q = q,
                        breaks = br, degree = 1, pos.psi = pos.psi, neg.psi = neg.psi,
                        pos.weights = sort(pos.weights, decreasing = TRUE), neg.weights = sort(neg.weights,
                                                                                               decreasing = TRUE), pos.size = sum(abs(wcoef[poscoef])),
                        neg.size = sum(abs(wcoef[negcoef])), 
                        tstat = tstat, 
                        pval = pvalz,
                        alpha = alpha, call = origcall, 
                        hasintercept = hasintercept
  )
  # attr(res, "class") <- c("qgcompfit",attr(res, "class"))

  res
}

