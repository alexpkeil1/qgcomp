
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
#' @param weights_name 
#' @param cluster not yet implemented
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
#' set.seed(50)
#' # linear model
#' dat <- data.frame(y=runif(50,-1,1), x1=runif(50), x2=runif(50), z=runif(50))
#' qgcomp.glm.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian())
#' qgcomp.tobit.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
#' # not intercept model
#' qgcomp.glm.noboot(f=y ~-1+ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian())
#' qgcomp.tobit.noboot(f=y ~-1+ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, dist="gaussian", left = -Inf, right = Inf)
qgcomp.tobit.noboot <- function (f, data, expnms = NULL, q = 4, breaks = NULL, id = NULL, 
                                  weights_name = NULL, cluster = NULL, alpha = 0.05, left = -Inf, right = Inf, ...) 
{
  weights <- data[weights_name]
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
  #tob = list(family = "gaussian", link = "identity", linkfun = identity)
  #class(tob) = "family"
  #fit[["family"]] = tob
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
  names(estb)[psiidx] <- c("psi1")
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

