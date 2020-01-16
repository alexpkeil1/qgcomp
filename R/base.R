
se_comb <- function(expnms, covmat, grad=NULL){
  #' @title Calculate standard error of weighted linear combination of random variables
  # given a vector of weights and a covariance matrix (not exported)
  #' @description This function uses the Delta method to calculate standard errors of linear
  #' functions of variables (similar to `lincom` in Stata). Generally, users will not need to 
  #' call this function directly.
  #' @details This function takes inputs of a set of exposure names (character vector)
  #' and a covariance matrix (with colnames/rownames that contain the full set
  #' of exposure names), as well as a possible `grad` parameter to calculate
  #' the variance of a weighted combination of the exposures in `expnms`, where the
  #' weights are based off of `grad` (which defaults to 1, so that this function
  #' yields the variance of a sum of all variables in `expnms`)
  #' 
  #' Here is simple version of the delta method for a linear combination of 
  #' three model coefficients:
  #' 
  #'  \eqn{f(\beta) = \beta_1 + \beta_2 + \beta_3}
  #' given gradient vector 
  #' \deqn{G = [d(f(\beta))/d\beta_1 = 1,
  #'   d(f(\beta))/d\beta_2 = 1,
  #'   d(f(\beta))/d\beta_3 = 1]}
  #' \eqn{t(G) Cov(\beta)  G} = delta method variance, where t() is the transpose operator
  #' 
  #' @param expnms a character vector with the names of the columns to be
  #' of interest in the covariance matrix for a which a standard error will be
  #' calculated (e.g. same as expnms in qgcomp fit)
  #' @param covmat covariance matrix for parameters, e.g. from a model or 
  #' bootstrap procedure
  #' @param grad the "weight" vector for calculating the contribution of each variable
  #' in expnms to the final standard error. For a linear combination, this is equal 
  #' to a vector of ones (and is set automatically). Or can be calculated via the 
  #' grad.poly procedure, in the case of coming up with proper weights when the combination
  #' of expnms derives from a polynomial function (as in qgcomp.boot with degree>1).
  #
  #' @examples
  #' vcov = rbind(c(1.2, .9),c(.9, 2.0))
  #' colnames(vcov) <- rownames(vcov) <- expnms <- c("x1", "x2")
  #' qgcomp:::se_comb(expnms, vcov, c(1, 0))^2 # returns the given variance
  #' qgcomp:::se_comb(expnms, vcov, c(1, 1)) # default linear MSM fit: all exposures
  #' # have equal weight
  #' qgcomp:::se_comb(expnms, vcov, c(.3, .1)) # used when one exposure contributes
  #'   # to the overall fit more than others  = d(msmeffect)/dx

  if(!is.matrix(covmat)) {
    nm <- names(covmat)
    covmat = matrix(covmat)
    colnames(covmat) <- nm
  }
  weightvec <- rep(0, dim(covmat)[1])
  # eventual extension: allow non-zero 'weights' such that the intervention
  # could correspond to 1 unit increases in some variables, and < 1 unit increases
  # in others
  if(is.null(grad)) weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- 1
  if(!is.null(grad)) weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- grad
  var <- weightvec %*% covmat %*% weightvec # delta method
  #var <- sum(wcovmat)
  sqrt(var)[1,,drop=TRUE] # should be a scalar
}

grad.poly <- function(intvals, degree){
  # returns matrix with each column referring
  if(degree==1){
    mat <- matrix(1, nrow=length(intvals), ncol=1)
  }else{
    mat <- matrix(1, nrow=length(intvals), ncol=degree)
    for(d in 2:degree){
      mat[,d] <- d*poly(intvals, degree = d-1, simple = TRUE, raw = TRUE)[,d-1]
    }
  }
  mat
}


quantize <- function (data, expnms, q=4, breaks=NULL) {
  #' @title create variables representing indicator functions with cutpoints defined
  #' by quantiles
  #' @description This function creates categorical variables in place of the
  #' exposure variables named in 'expnms.' For example, a continuous exposure
  #' 'x1' will be replaced in the output data by another 'x1' that takes on values
  #' 0:(q-1), where, for example, the value 1 indicates that the original x1 value
  #' falls between the first and the second quantile.
  #' @details This function is a vectorized version of `quantile_f` from the `gWQS` 
  #' package that also allows the use of externally defined breaks
  #' @param data a data frame
  #' @param expnms a character vector with the names of  the columns to be
  #' quantized
  #' @param q integer, number of quantiles used in creating quantized variables
  #' @param breaks (optional) list of (equal length) numeric vectors that 
  #' characterize the minimum value of each category for which to 
  #' break up the variables named in expnms. This is an alternative to using 'q'
  #' to define cutpoints.
  #' @concept variance mixtures
  #' @import stats
  #' @export
  #' @examples
  #' set.seed(1232)
  #' dat = data.frame(y=runif(100), x1=runif(100), x2=runif(100), z=runif(100))
  #' qdata = quantize(data=dat, expnms=c("x1", "x2"), q=4)
  #' table(qdata$data$x1)
  #' table(qdata$data$x2)
  #' summary(dat[c("y", "z")]);summary(qdata$data[c("y", "z")]) # not touched
  #' dat = data.frame(y=runif(100), x1=runif(100), x2=runif(100), z=runif(100))
  #' # using 'breaks' requires specifying min and max (the qth quantile)
  #' # example with theoretical quartiles (could be other relevant values)
  #' qdata2 = quantize(data=dat, expnms=c("x1", "x2"), 
  #'    breaks=list(c(-1e64, .25, .5, .75, 1e64), 
  #'                c(-1e64, .25, .5, .75, 1e64)
  #'                ))
  #' table(qdata2$data$x1)
  #' table(qdata2$data$x2)
    e <- new.env()
    e$retbr <- list()
    qt <- function(i){
      # not exported
        datmat <- as.numeric(unlist(data[, expnms[i]]))
        if(!is.null(breaks)){
          # prioritize breaks if given by user
          br  <- breaks[[i]]
          e$retbr[[i]] <- breaks[[i]]
        }else{
          br <- unique(quantile(datmat, probs = seq(0, 1, by = 1 / q), na.rm = TRUE))
          br[1] <- -1e64
          br[length(br)] <- 1e64
          e$retbr[[i]] <- br 
        }
        cut(datmat, breaks = br, labels = FALSE,
             include.lowest = TRUE) - 1
    }
    if(length(expnms)==1){
      data[, expnms] <- qt(1)
    }else{
      data[, expnms] <- sapply(1:length(expnms), qt)
    }
    return(list(data=data, breaks=e$retbr))
}

checknames <- function(terms){
  #' @title check for valid model terms if 'expnms' parameter not given
  #' @description This is an internal function called by \code{\link[qgcomp]{qgcomp}},
  #'  \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp.noboot}},
  #'  but is documented here for clarity. Generally, users will not need to call
  #'  this function directly. This function tries to determine whether there are 
  #'  non-linear terms in the underlying model, which helps infer whether the 
  #'  appropriate function is called, and whether more explicit function calls
  #'  are needed.
  #' @param terms model terms from attr(terms(modelfunction, data), "term.labels")
  nonlin <- ifelse(sum(grep("\\(|\\:|\\^", terms))>0, TRUE, FALSE)
  if(nonlin){
    return(FALSE)
  }else{
    return(TRUE)
  }
}


msm.fit <- function(f, qdata, intvals, expnms, rr=TRUE, main=TRUE, degree=1, id=NULL, bayes, ...){
  #' @title fitting marginal structural model (MSM) based on g-computation with
  #' quantized exposures
  #' @description This is an internal function called by \code{\link[qgcomp]{qgcomp}},
  #'  \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp.noboot}},
  #'  but is documented here for clarity. Generally, users will not need to call
  #'  this function directly.
  #' @details This function first computes expected outcomes under hypothetical
  #' interventions to simultaneously set all exposures to a specific quantile. These
  #' predictions are based on g-computation, where the exposures are `quantized',
  #' meaning that they take on ordered integer values according to their ranks,
  #' and the integer values are determined by the number of quantile cutpoints used.
  #' The function then takes these expected outcomes and fits an additional model
  #' (a marginal structural model) with the expected outcomes as the outcome and
  #' the intervention value of the exposures (the quantile integer) as the exposure.
  #' Under causal identification assumptions and correct model specification,
  #' the MSM yields a causal exposure-response representing the incremental
  #' change in the expected outcome given a joint intervention on all exposures.
  #' @param f an r formula representing the conditional model for the outcome, given all
  #' exposures and covariates. Interaction terms that include exposure variables
  #' should be represented via the \code{\link[base]{AsIs}} function
  #' @param qdata a data frame with quantized exposures
  #' @param expnms a character vector with the names of the columns in qdata that represent
  #' the exposures of interest (main terms only!)
  #' @param intvals sequence, the sequence of integer values that the joint exposure 
  #' is 'set' to for estimating the msm. For quantile g-computation, this is just 
  #' 0:(q-1), where q is the number of quantiles of exposure.
  #' @param rr logical, estimate log(risk ratio) (family='binomial' only)
  #' @param main logical, internal use: produce estimates of exposure effect (psi)
  #'  and expected outcomes under g-computation and the MSM
  #' @param degree polynomial basis function for marginal model (e.g. degree = 2
  #'  allows that the relationship between the whole exposure mixture and the outcome
  #'  is quadratic. Default=1)
  #' @param id (optional) NULL, or variable name indexing individual units of 
  #' observation (only needed if analyzing data with multiple observations per 
  #' id/cluster)
  #' @param bayes use underlying Bayesian model (`arm` package defaults). Results
  #' in penalized parameter estimation that can help with very highly correlated 
  #' exposures. Note: this does not lead to fully Bayesian inference in general, 
  #' so results should be interpereted as frequentist.
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @concept variance mixtures
  #' @import stats arm
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(200), x1=runif(200), x2=runif(200), z=runif(200))
  #' X <- c('x1', 'x2')
  #' qdat <- quantize(dat, X, q=4)$data
  #' mod <- qgcomp:::msm.fit(f=y ~ z + x1 + x2 + I(x1*x2), 
  #'         expnms = c('x1', 'x2'), qdata=qdat, intvals=1:4, bayes=FALSE)
  #' summary(mod$fit) # outcome regression model
  #' summary(mod$msmfit) # msm fit (variance not valid - must be obtained via bootstrap)
    if(is.null(id)) {
      id <- "id__"
      qdata$id__ <- 1:dim(qdata)[1]
    }
    # conditional outcome regression fit
    nidx = which(!(names(qdata) %in% id))
    if(!bayes) fit <- glm(f, data = qdata[,nidx,drop=FALSE], ...)
    if(bayes){
      requireNamespace("arm")
      fit <- bayesglm(f, data = qdata[,nidx,drop=FALSE], ...)
    } 
    if(fit$family$family=="gaussian") rr=FALSE
    ### 
    # get predictions (set exposure to 0,1,...,q-1)
    if(is.null(intvals)){
      intvals <- (1:length(table(qdata[expnms[1]]))) - 1
    }
    predit <- function(idx){
      newdata <- qdata
      newdata[,expnms] <- idx
      suppressWarnings(predict(fit, newdata=newdata, type='response'))
    }
    predmat = lapply(intvals, predit)
    # fit MSM using g-computation estimates of expected outcomes under joint 
    #  intervention
    nobs <- dim(qdata)[1]
    msmdat <- data.frame(
      Ya = unlist(predmat),
      psi = rep(intvals, each=nobs))
    # to do: allow functional form variations for the MSM via specifying the model formula
    if(bayes){
      if(!rr) suppressWarnings(msmfit <- bayesglm(Ya ~ poly(psi, degree=degree, raw=TRUE), data=msmdat,...))
      if(rr)  suppressWarnings(msmfit <- bayesglm(Ya ~ poly(psi, degree=degree, raw=TRUE), data=msmdat, family=binomial(link='log'), start=rep(-0.0001, degree+1)))
    }
    if(!bayes){
      if(!rr) suppressWarnings(msmfit <- glm(Ya ~ poly(psi, degree=degree, raw=TRUE), data=msmdat,...))
      if(rr)  suppressWarnings(msmfit <- glm(Ya ~ poly(psi, degree=degree, raw=TRUE), data=msmdat, family=binomial(link='log'), start=rep(-0.0001, degree+1)))
    }
    res <- list(fit=fit, msmfit=msmfit)
    if(main) {
      res$Ya <- msmdat$Ya   # expected outcome under joint exposure, by gcomp
      res$Yamsm <- predict(msmfit, type='response')
      res$A <- msmdat$psi # joint exposure (0 = all exposures set category with 
       # upper cut-point as first quantile)
    }
    res
}


qgcomp.noboot <- function(f, data, expnms=NULL, q=4, breaks=NULL, id=NULL, alpha=0.05, bayes=FALSE, ...){
  #' @title estimating the parameters of a marginal structural model (MSM) based on 
  #' g-computation with quantized exposures
  #'
  #' @description This function mimics the output of a weighted quantile sums regression in 
  #' large samples. 
  #' 
  #' @details For continuous outcomes, under a linear model with no 
  #' interaction terms, this is equivalent to g-computation of the effect of
  #' increasing every exposure by 1 quantile. For binary outcomes
  #' outcomes, this yields a conditional log odds ratio representing the 
  #' change in the expected conditional odds (conditional on covariates)
  #' from increasing every exposure by 1 quantile. In general, the latter 
  #' quantity is not equivalent to g-computation estimates. Hypothesis test
  #' statistics and 95% confidence intervals are based on using the delta
  #' estimate variance of a linear combination of random variables.
  #' 
  #' @param f R style formula
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
  #' @param id (optional) NULL, or variable name indexing individual units of 
  #' observation (only needed if analyzing data with multiple observations per 
  #' id/cluster)
  #' @param alpha alpha level for confidence limit calculation
  #' @param bayes use underlying Bayesian model (`arm` package defaults). Results
  #' in penalized parameter estimation that can help with very highly correlated 
  #' exposures. Note: this does not lead to fully Bayesian inference in general, 
  #' so results should be interpereted as frequentist.
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and information on the 
  #'  weights/standardized coefficients in the positive (pos.weights) and 
  #'  negative (nweight) directions.
  #' @concept variance mixtures
  #' @import stats arm
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  if (is.null(expnms)) {
    expnms <- attr(terms(f, data = data), "term.labels")
    message("Including all model terms as exposures of interest\n")      
  }
  lin = checknames(expnms)
  if(!lin) stop("Model appears to be non-linear: use qgcomp.boot instead")
  if (!is.null(q) | !is.null(breaks)){
      ql <- quantize(data, expnms, q, breaks)
      qdata <- ql$data
      br <- ql$breaks
    } else{
      qdata <- data
      br <- breaks
      }
    if(is.null(id)) {
      # not yet implemented
      id = "id__"
      qdata$id__ = 1:dim(qdata)[1]
    }
    if(!bayes) fit <- glm(f, data = qdata[,!(names(qdata) %in% id), drop=FALSE], ...)
    if(bayes){
      requireNamespace("arm")
      fit <- bayesglm(f, data = qdata[,!(names(qdata) %in% id), drop=FALSE], ...)
    }
    mod <- summary(fit)
    if(length(setdiff(expnms, rownames(mod$coefficients)))>0){
      stop("Model aliasing occurred, likely due to perfectly correlated quantized exposures. 
           Try one of the following:
             1) set 'bayes' to TRUE in the qgcomp function (recommended)
             2) set 'q' to a higher value in the qgcomp function (recommended)
             3) check correlation matrix of exposures, and drop all but one variable in each highly correlated set  (not recommended)
           ")
    }
    estb <- c(fit$coefficients[1], sum(mod$coefficients[expnms,1, drop=TRUE]))
    seb <- c(sqrt(mod$cov.scaled[1,1]), se_comb(expnms, covmat = mod$cov.scaled))
    tstat <- estb / seb
    df <- mod$df.null - length(expnms)
    pval <- 2 - 2 * pt(abs(tstat), df = df)
    pvalz <- 2 - 2 * pnorm(abs(tstat))
    ci <- cbind(estb + seb * qnorm(alpha / 2), estb + seb * qnorm(1 - alpha / 2))
    # 'weights'
    wcoef <- fit$coefficients[expnms]
    names(wcoef) <- gsub("_q", "", names(wcoef))
    pos.coef <- which(wcoef > 0)
    neg.coef <- which(wcoef <= 0)
    pos.weights <- abs(wcoef[pos.coef]) / sum(abs(wcoef[pos.coef]))
    neg.weights <- abs(wcoef[neg.coef]) / sum(abs(wcoef[neg.coef]))
    # 'post-hoc' positive and negative estimators 
    # similar to constrained gWQS
    pos.psi <- sum(wcoef[pos.coef])
    neg.psi <- sum(wcoef[neg.coef])
    #nmpos <- names(pos.weights)
    #nmneg <- names(neg.weights)
    #se.pos.psi <- se_comb(nmpos, covmat = mod$cov.scaled)
    #se.neg.psi <- se_comb(nmneg, covmat = mod$cov.scaled)
    qx <- qdata[, expnms]
    names(qx) <- paste0(names(qx), "_q")
    res <- list(
      qx = qx, fit = fit, 
      psi = estb[-1], var.psi = seb[-1] ^ 2, covmat.psi=c('psi1' = seb[-1]^2), ci = ci[-1,],
      coef = estb, var.coef = seb ^ 2, covmat.coef=c('(Intercept)' = seb[1]^2, 'psi1' = seb[2]^2), 
      ci.coef = ci,
      expnms=expnms, q=q, breaks=br, degree=1,
      pos.psi = pos.psi, neg.psi = neg.psi,
      pos.weights = sort(pos.weights, decreasing = TRUE),
      neg.weights = sort(neg.weights, decreasing = TRUE), 
      pos.size = sum(abs(wcoef[pos.coef])),
      neg.size = sum(abs(wcoef[neg.coef])),
      bootstrap=FALSE,
      cov.yhat=NULL
    )
      if(fit$family$family=='gaussian'){
        res$tstat <- tstat
        res$df <- df
        res$pval <- pval
      }
      if(fit$family$family=='binomial'){
        res$zstat <- tstat
        res$pval <- pvalz
      }
    attr(res, "class") <- "qgcompfit"
    res
}

#TODO: explain (log) better - here and in the noboot
qgcomp.boot <- function(f, data, expnms=NULL, q=4, breaks=NULL, id=NULL, alpha=0.05, B=200, 
                        rr=TRUE, degree=1, seed=NULL, bayes=FALSE, parallel=FALSE, ...){
  #' @title estimating the parameters of a marginal structural model (MSM) based on 
  #' g-computation with quantized exposures
  #'  
  #' @description This function yields population average effect estimates for 
  #'   both continuous and binary outcomes
  #'  
  #' @details Estimates correspond to the average expected change in the
  #'  (log) outcome per quantile increase in the joint exposure to all exposures 
  #'  in `expnms'. Test statistics and confidence intervals are based on 
  #'  a non-parametric bootstrap, using the standard deviation of the bootstrap
  #'  estimates to estimate the standard error. The bootstrap standard error is 
  #'  then used to estimate Wald-type confidence intervals. Note that no bootstrapping
  #'  is done on estimated quantiles of exposure, so these are treated as fixed
  #'  quantities
  #'
  #' @param f R style formula
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
  #' @param id (optional) NULL, or variable name indexing individual units of 
  #' observation (only needed if analyzing data with multiple observations per 
  #' id/cluster)
  #' @param alpha alpha level for confidence limit calculation
  #' @param B integer: number of bootstrap iterations (this should typically be
  #' >=200, though it is set lower in examples to improve run-time).
  #' @param rr logical: if using binary outcome and rr=TRUE, qgcomp.boot will 
  #'   estimate risk ratio rather than odds ratio
  #' @param degree polynomial basis function for marginal model (e.g. degree = 2
  #'  allows that the relationship between the whole exposure mixture and the outcome
  #'  is quadratic.
  #' @param seed integer or NULL: random number seed for replicable bootstrap results
  #' @param bayes use underlying Bayesian model (`arm` package defaults). Results
  #' in penalized parameter estimation that can help with very highly correlated 
  #' exposures. Note: this does not lead to fully Bayesian inference in general, 
  #' so results should be interpereted as frequentist.
  #' @param parallel use (safe) parallel processing from the future and future.apply packages
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.noboot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and information on the 
  #'  marginal structural model (msmfit) used to estimate the final effect
  #'  estimates.
  #' @concept variance mixtures
  #' @import stats
  #' @export
  #' @examples
  #' set.seed(30)
  #' # continuous outcome
  #' dat <- data.frame(y=rnorm(100), x1=runif(100), x2=runif(100), z=runif(100))
  #' # Conditional linear slope
  #' qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4, family=gaussian())
  #' # Marginal linear slope (population average slope, for a purely linear, 
  #' #  additive model this will equal the conditional)
  #'  \donttest{
  #' qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4, 
  #'   family=gaussian(), B=200) # B should be at least 200 in actual examples
  #'   
  #' # Population average mixture slope which accounts for non-linearity and interactions
  #' qgcomp.boot(y ~ z + x1 + x2 + I(x1^2) + I(x2*x1), family="gaussian", 
  #'  expnms = c('x1', 'x2'), data=dat, q=4, B=200)
  #'  
  #' # binary outcome
  #' dat <- data.frame(y=rbinom(50,1,0.5), x1=runif(50), x2=runif(50), z=runif(50))
  #' 
  #' # Conditional mixture OR
  #' qgcomp.noboot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), 
  #'   data=dat, q=2)
  #'   
  #' #Marginal mixture OR (population average OR - in general, this will not equal the 
  #' # conditional mixture OR due to non-collapsibility of the OR)
  #' qgcomp.boot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), 
  #'   data=dat, q=2, B=3)
  #'   
  #' # Population average mixture RR
  #' qgcomp.boot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), 
  #'   data=dat, q=2, rr=TRUE, B=3)
  #'   
  #' # Population average mixture RR, indicator variable representation of x2
  #' # note that I(x==...) operates on the quantile-based category of x,
  #' # rather than the raw value
  #' res = qgcomp.boot(y ~ z + x1 + I(x2==1) + I(x2==2) + I(x2==3), 
  #'   family="binomial", expnms = c('x1', 'x2'), data=dat, q=4, rr=TRUE, B=200)
  #' res$fit  
  #' plot(res)
  #' 
  #' # now add in a non-linear MSM
  #' res2 = qgcomp.boot(y ~ z + x1 + I(x2==1) + I(x2==2) + I(x2==3), 
  #'   family="binomial", expnms = c('x1', 'x2'), data=dat, q=4, rr=TRUE, B=200, 
  #'   degree=2)
  #' res2$fit  
  #' res2$msmfit  
  #' plot(res2)
  #' # Log risk ratio per one IQR change in all exposures (not on quantile basis)
  #' dat$x1iqr <- dat$x1/with(dat, diff(quantile(x1, c(.25, .75))))
  #' dat$x2iqr <- dat$x2/with(dat, diff(quantile(x2, c(.25, .75))))
  #' # note that I(x>...) now operates on the untransformed value of x,
  #' # rather than the quantized value
  #' res2 = qgcomp.boot(y ~ z + x1iqr + I(x2iqr>0.1) + I(x2>0.4) + I(x2>0.9), 
  #'   family="binomial", expnms = c('x1iqr', 'x2iqr'), data=dat, q=NULL, rr=TRUE, B=200, 
  #'   degree=2)
  #' res2
  #' # using parallel processing
  #'
  #' qgcomp.boot(y ~ z + x1iqr + I(x2iqr>0.1) + I(x2>0.4) + I(x2>0.9), 
  #'   family="binomial", expnms = c('x1iqr', 'x2iqr'), data=dat, q=NULL, rr=TRUE, B=200, 
  #'   degree=2, parallel=TRUE)
  #' }
  # character names of exposure mixture components
    if(is.null(seed)) seed = round(runif(1, min=0, max=1e8))
    if (is.null(expnms)) {
      expnms <- attr(terms(f, data = data), "term.labels")
      message("Including all model terms as exposures of interest\n")      
    }
    lin = checknames(expnms)
    if(!lin) stop("Model appears to be non-linear and I'm having trouble parsing it: 
                  please use `expnms` parameter to define the variables making up the exposure")
    if (!is.null(q) & !is.null(breaks)){
      # if user specifies breaks, prioritize those
      q <- NULL
    }
    if (!is.null(q) | !is.null(breaks)){
      ql <- quantize(data, expnms, q, breaks)
      qdata <- ql$data
      br <- ql$breaks
      if(is.null(q)){
        # rare scenario with user specified breaks and q is left at NULL
        nvals <- length(br[[1]])-1
      } else{
        nvals <- q
      }
      intvals <- (1:nvals)-1
    } else {
      # if( is.null(breaks) & is.null(q)) # also includes NA
      qdata <- data
      # if no transformation is made (no quantiles, no breaks given)
      # then draw distribution values from quantiles of all the exposures
      # pooled together
      # TODO: allow user specification of this
      message("\nNote: using quantiles of all exposures combined in order to set 
          proposed intervention values for overall effect (25th, 50th, 75th %ile)\n")
      intvals = as.numeric(quantile(unlist(data[,expnms]), c(.25, .5, .75)))
      br <- NULL
    }
    if(is.null(id)) {
      id <- "id__"
      qdata$id__ <- 1:dim(qdata)[1]
    }
    ###
    msmfit <- msm.fit(f, qdata, intvals, expnms, rr, main=TRUE,degree=degree, id=id, bayes, ...)
    # main estimate  
    #estb <- as.numeric(msmfit$msmfit$coefficients[-1])
    estb <- as.numeric(msmfit$msmfit$coefficients)
    #bootstrap to get std. error
    nobs <- dim(qdata)[1]
    nids <- length(unique(qdata[,id, drop=TRUE]))
    starttime = Sys.time()
    psi.only <- function(i=1, f=f, qdata=qdata, intvals=intvals, expnms=expnms, rr=rr, degree=degree, nids=nids, id=id, ...){
      if(i==2){
        timeiter = as.numeric(Sys.time() - starttime)
        if((timeiter*B/60)>0.5) message(paste0("Expected time to finish: ", round(B*timeiter/60, 2), " minutes \n"))
      }
      bootids <- data.frame(temp=sort(sample(unique(qdata[,id, drop=TRUE]), nids, replace = TRUE)))
      names(bootids) <- id
      qdata_ <- merge(qdata,bootids, by=id, all.x=FALSE, all.y=TRUE)
      ft = msm.fit(f, qdata_, intvals, expnms, rr, main=FALSE, degree, id, bayes,
              ...)
      yhatty = data.frame(yhat=predict(ft$msmfit), psi=ft$msmfit$data[,"psi"])
      as.numeric(
        c(ft$msmfit$coefficients, with(yhatty, tapply(yhat, psi, mean)))
        #msm.fit(f, qdata_, intvals, expnms, rr, main=FALSE, degree, id, bayes,
        #        ...)$msmfit$coefficients[-1]
      )
    }
    set.seed(seed)
    if(parallel){
      Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
      future::plan(strategy = future::multiprocess)
      bootsamps <- future.apply::future_sapply(X=1:B, FUN=psi.only,f=f, qdata=qdata, intvals=intvals, 
                          expnms=expnms, rr=rr, degree=degree, nids=nids, id=id, ...)
      
      future::plan(future::sequential)
    }else{
      bootsamps <- sapply(X=1:B, FUN=psi.only,f=f, qdata=qdata, intvals=intvals, 
                          expnms=expnms, rr=rr, degree=degree, nids=nids, id=id, ...)
      
    }
    #if(is.null(dim(bootsamps))) {
    #  seb <- sd(bootsamps)
    #  covmat <- var(bootsamps)
    #  names(covmat) <- 'psi1'
    #}else{
    hats = t(bootsamps[-c(1:(degree+1)),])
    cov.yhat = cov(hats)
    bootsamps = bootsamps[1:(degree+1),]
    seb <- apply(bootsamps, 1, sd)
      covmat <- cov(t(bootsamps))
      colnames(covmat) <- rownames(covmat) <- c("(intercept)", paste0("psi", 1:(nrow(bootsamps)-1)))
    #}
    tstat <- estb / seb
    df <- nobs - length(attr(terms(f, data = data), "term.labels")) - 1 - degree # df based on obs - gcomp terms - msm terms
    pval <- 2 - 2 * pt(abs(tstat), df = df)
    pvalz <- 2 - 2 * pnorm(abs(tstat))
    ci <- cbind(estb + seb * qnorm(alpha / 2), estb + seb * qnorm(1 - alpha / 2))
    # 'weights' not applicable in this setting, generally (i.e. if using this function for non-linearity, 
    #   then weights will vary with level of exposure)
    qx <- qdata[, expnms]
    res <- list(
      qx = qx, fit = msmfit$fit, msmfit = msmfit$msmfit, 
      psi = estb[-1], var.psi = seb[-1] ^ 2, covmat.psi=covmat[-1,-1, drop=FALSE], ci = ci[-1,],
      coef = estb, var.coef = seb ^ 2, covmat.coef=covmat, ci.coef = ci,
      expnms=expnms, q=q, breaks=br, degree=degree,
      pos.psi = NULL, neg.psi = NULL, 
      pos.weights = NULL, neg.weights = NULL, pos.size = NULL,neg.size = NULL, bootstrap=TRUE,
      y.expected=msmfit$Ya, y.expectedmsm=msmfit$Yamsm, index=msmfit$A,
      bootsamps = bootsamps,
      cov.yhat=cov.yhat
    )
      if(msmfit$fit$family$family=='gaussian'){
        res$tstat <- tstat
        res$df <- df
        res$pval <- pval
      }
      if(msmfit$fit$family$family=='binomial'){
        res$zstat <- tstat
        res$pval <- pvalz
      }
    attr(res, "class") <- "qgcompfit"
    res
}


qgcomp <- function(f,data=data,family=gaussian(),rr=TRUE,...){
  #' @title estimation of quantile g-computation fit
  #' 
  #' @description This function automatically selects between qgcomp.noboot, qgcomp.boot,
  #'  qgcomp.cox.noboot, and qgcomp.cox.boot
  #'  for the most efficient approach to estimate the average expected 
  #'  change in the (log) outcome per quantile increase in the joint 
  #'  exposure to all exposures in `expnms', given the underlying model. For example,
  #'  if the underlying model (specified by the formula `f`) is a linear model with all
  #'  linear terms for exposure, then `qgcomp.noboot`` will be called to fit the model. Non-linear
  #'  terms or requesting the risk ratio for binomial outcomes will result in the `qgcomp.boot`
  #'  function being called. For a given linear model, boot and noboot versions will give identical
  #'  inference, though when using survival outcomes, the `boot` version uses simulation based 
  #'  inference, which can vary from the `nonboot` version due to simulation error (which can 
  #'  be minimized via setting the MCsize parameter very large - see 
  #'  \code{\link[qgcomp]{qgcomp.cox.boot}} for details).
  #'
  #' @param f R style formula (may include survival outcome via \code{\link[survival]{Surv}})
  #' @param data data frame
  #' @param family `gaussian()`, `binomial()`, `cox()`
  #' @param rr logical: if using binary outcome and rr=TRUE, qgcomp.boot will 
  #' estimate risk ratio rather than odds ratio. Note, to get population average 
  #' effect estimates for a binary outcome, set rr=TRUE (default: ORs are generally not
  #' of interest as population average effects, so if rr=FALSE then a conditional
  #' OR will be estimated, which cannot be interpreted as a population average
  #' effect
  #' @param ... arguments to qgcomp.noboot or qgcomp.boot (e.g. q)
  #' @seealso \code{\link[qgcomp]{qgcomp.noboot}}, \code{\link[qgcomp]{qgcomp.boot}}, 
  #'  \code{\link[qgcomp]{qgcomp.cox.noboot}} and \code{\link[qgcomp]{qgcomp.cox.boot}}
  #'  (\code{\link[qgcomp]{qgcomp}} is just a wrapper for these functions)
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and possibly information on the 
  #'  marginal structural model (msmfit) used to estimate the final effect
  #'  estimates (qgcomp.boot, qgcomp.cox.boot only). 
  #'  If appropriate, weights are also reported, which represent the proportion
  #'  of a directional (positive/negative) effect that is accounted for by
  #'  each exposure.
  #' @concept variance mixtures
  #' @import stats
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' qgcomp.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=10, seed=125)
  #' # automatically selects appropriate method
  #' qgcomp(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' # note for binary outcome this will 
  #' dat <- data.frame(y=rbinom(50, 1, 0.5), x1=runif(50), x2=runif(50), z=runif(50))
  #' qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=binomial())
  #' qgcomp.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=10, seed=125, 
  #'   family=binomial())
  #' # automatically selects appropriate method
  #' qgcomp(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=binomial())
  #' qgcomp(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=binomial(), rr=TRUE)
  # f = y ~ factor(x1) + x2
  #' 
  #' #survival objects
  #' set.seed(50)
  #' N=200
  #' dat <- data.frame(time=(tmg <- pmin(.1,rweibull(N, 10, 0.1))), 
  #'                 d=1.0*(tmg<0.1), x1=runif(N), x2=runif(N), z=runif(N))
  #' expnms=paste0("x", 1:2)
  #' f = survival::Surv(time, d)~x1 + x2
  #' qgcomp(f, expnms = expnms, data = dat)
  #' # note that in the survival models, MCsize should be set to a large number
  #' #  such that results are repeatable (within an error tolerance such as 2 significant digits)
  #' # if you run them under different  seed values
  #' f = survival::Surv(time, d)~x1 + x2 + x1:x2
  #' qgcomp(f, expnms = expnms, data = dat, B=10, MCsize=100)
  requireNamespace("survival")
  issurv = survival::is.Surv(eval(attr(terms(f, data = data), "variables")[[2]], envir = data))
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    if(issurv){
      message("Survival model\n")
    }else{
      print(family)
      stop("'family' not recognized")
    }
  }
  if(!(family$family %in% c("binomial")) & rr) {
    #warning("'rr=TRUE' is for bimomial family only, setting rr=FALSE")
    rr = FALSE
  }
  terms <- attr(terms(f,data=data), 'term.labels')
  doboot <- !checknames(terms)
  if(issurv & doboot ){
    res <- qgcomp.cox.boot(f=f,data=data,...)
  } else if(issurv & !doboot ){
    res <- qgcomp.cox.noboot(f=f,data=data,...)
  } else if(rr | doboot){
    res <- qgcomp.boot(f=f,data=data,family=family,rr=rr,...)
  }else{
    res <- qgcomp.noboot(f=f,data=data,family=family,...)
  }
  res
}

coef.qgcompfit <- function(object, ...){
  #' @importFrom stats coef
  #' @export
  object$coef
}

print.qgcompfit <- function(x, ...){
  #' @title default printing method for a qgcompfit object
  #' 
  #' @description Gives variable output depending on whether `qgcomp.noboot` or `qgcomp.boot`
  #' is called. For `qgcomp.noboot` will output final estimate of joint exposure
  #' effect (similar to the 'index' effect in weighted quantile sums), as well
  #' as estimates of the 'weights' (standardized coefficients). For `qgcomp.boot`,
  #' the marginal effect is given, but no weights are reported since this approach
  #' generally incorporates non-linear models with interaction terms among exposures,
  #' which preclude weights with any useful interpretation.
  #' 
  #' @param x "qgcompfit" object from `qgcomp`, `qgcomp.noboot` or `qgcomp.boot` 
  #' function
  #' @param ... unused
  #' @seealso \code{\link[qgcomp]{qgcomp.noboot}}, \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @concept variance mixtures
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' obj1 <- qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' obj2 <- qgcomp.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=10, seed=125)
  #' # does not need to be explicitly called, but included here for clarity
  #' print(obj1)
  #' print(obj2)
  fam <- x$fit$family$family
  if(is.null(fam)){
    printZI(x)
    return(invisible(x))
  }
  if(!is.null(x$pos.size)) {
    cat(paste0("Scaled effect size (positive direction, sum of positive coefficients = ", signif(x$pos.size, 3) , ")\n"))
    if (length(x$pos.weights) > 0) {
      print(x$pos.weights, digits = 3)
    } else cat("None\n")
    cat("\n")
  }
  if(!is.null(x$neg.size)) {
    cat(paste0("Scaled effect size (negative direction, sum of negative coefficients = ", signif(-x$neg.size, 3) , ")\n"))
    if (length(x$neg.weights) > 0) {
      print(x$neg.weights, digits = 3)
    } else cat("None\n")
    cat("\n")
  }
  if (fam == "binomial"){
    estimand <- 'OR'
    if(x$bootstrap && x$msmfit$family$link=='log') estimand = 'RR'
    cat(paste0("Mixture log(",estimand,")", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "Z"
    rnm = c("(Intercept)", c(paste0('psi',1:max(1, length(coef(x))-1))))
  }
  if (fam == "gaussian"){
    cat(paste0("Mixture slope parameters", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "t"
    rnm = c("(Intercept)", c(paste0('psi',1:max(1, length(coef(x))-1))))
  }
  if (fam == "cox"){
    cat(paste0("Mixture log(hazard ratio)", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "Z"
    rnm = c(paste0('psi',1:max(1, length(coef(x)))))
  }
  if(is.null(dim(x$ci.coef))){
    pdat <- cbind(Estimate=coef(x), "Std. Error"=sqrt(x$var.coef), "Lower CI"=x$ci.coef[1], "Upper CI"=x$ci.coef[2], "test"=x$zstat, "Pr(>|z|)"=x$pval)
    colnames(pdat)[5] = eval(paste(testtype, "value"))
    rownames(pdat) <- rnm
    printCoefmat(pdat,has.Pvalue=TRUE,tst.ind=5L,signif.stars=FALSE, cs.ind=1L:2)
  } else{
    pdat <- cbind(Estimate=coef(x), "Std. Error"=sqrt(x$var.coef), "Lower CI"=x$ci.coef[,1], "Upper CI"=x$ci.coef[,2], "test"=x$zstat, "Pr(>|z|)"=x$pval)
    colnames(pdat)[5] = eval(paste(testtype, "value"))
    rownames(pdat) <- rnm
    printCoefmat(pdat,has.Pvalue=TRUE,tst.ind=5L,signif.stars=FALSE, cs.ind=1L:2)
  }
}

summary.qgcompfit <- function(object, ...){
  #' @export
  print(object)
}

family.qgcompfit <- function(object, ...){
  #' @export
  object$fit$family
}


plot.qgcompfit <- function(x, 
                           suppressprint=FALSE, 
                           pointwisebars=TRUE, 
                           modelfitline=TRUE, 
                           modelband=TRUE, 
                           flexfit=TRUE, 
                           pointwiseref = 1,
                           ...){
  #' @title plot.qgcompfit: default plotting method for a qgcompfit object
  #'
  #' @description Plot a quantile g-computation object. For qgcomp.noboot, this function will
  #' create a butterfly plot of weights. For qgcomp.boot, this function will create
  #' a box plot with smoothed line overlaying that represents a non-parametric
  #' fit of a model to the expected outcomes in the population at each quantile
  #' of the joint exposures (e.g. '1' represents 'at the first quantile for
  #' every exposure')
  #' 
  #' @param x "qgcompfit" object from `qgcomp.noboot`,  `qgcomp.boot`, 
  #'   `qgcomp.cox.noboot`,  `qgcomp.cox.boot`, `qgcomp.zi.noboot` or `qgcomp.zi.boot` functions
  #' @param suppressprint If TRUE, suppresses the plot, rather than printing it 
  #'   by default (it can be saved as a ggplot2 object (or list of ggplot2 objects if x is from a zero-
  #'   inflated model) and used programmatically)
  #'   (default = FALSE)
  #' @param pointwisebars (boot.gcomp only) If TRUE, adds 95\%  error bars for pointwise comparisons
  #' of E(Y|joint exposure) to the smooth regression line plot
  #' @param modelfitline (boot.gcomp only) If TRUE, adds fitted (MSM) regression line
  #' of E(Y|joint exposure) to the smooth regression line plot
  #' @param modelband If TRUE, adds 95\% prediction bands for E(Y|joint exposure) (the MSM fit)
  #' @param flexfit (boot.gcomp only) if TRUE, adds flexible interpolation of predictions from 
  #' underlying (conditional) model
  #' @param pointwiseref (boot.gcomp only) integer: which category of exposure (from 1 to q) 
  #' should serve as the referent category for pointwise comparisons? (default=1)
  #' @param ... unused
  #' @seealso \code{\link[qgcomp]{qgcomp.noboot}}, \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @import ggplot2 grid gridExtra
  #' @export
  #' @examples
  #' set.seed(12)
  #' dat <- data.frame(x1=(x1 <- runif(100)), x2=runif(100), x3=runif(100), z=runif(100),
  #'                   y=runif(100)+x1+x1^2)
  #' ft <- qgcomp.noboot(y ~ z + x1 + x2 + x3, expnms=c('x1','x2','x3'), data=dat, q=4)
  #' ft
  #' plot(ft)
  #' # examining fit
  #' plot(ft$fit, which=1) # residual vs. fitted is not straight line!
  #' 
  #' # using non-linear outcome model
  #' ft2 <- qgcomp.boot(y ~ z + x1 + x2 + x3 + I(x1*x1), expnms=c('x1','x2','x3'), 
  #' data=dat, q=4, B=10)
  #' ft2
  #' plot(ft2$fit, which=1) # much better looking fit diagnostics suggests
  #' # it is better to include interaction term for x
  #' plot(ft2) # the msm predictions don't match up with a smooth estimate
  #' # of the expected outcome, so we should consider a non-linear MSM

  #' # using non-linear marginal structural model
  #' ft3 <- qgcomp.boot(y ~ z + x1 + x2 + x3 + I(x1*x1), expnms=c('x1','x2','x3'), 
  #' data=dat, q=4, B=10, degree=2)
  #' # plot(ft3$fit, which=1) - not run - this is identical to ft2 fit
  #' plot(ft3) # the MSM estimates look much closer to the smoothed estimates
  #' # suggesting the non-linear MSM fits the data better and should be used
  #' # for inference about the effect of the exposure
  requireNamespace("ggplot2")
  requireNamespace("grid")
  requireNamespace("gridExtra")
  ymin <- ymax <- w <- v <- NULL # appease R CMD check
  theme_butterfly_l <- list(theme(
    legend.position = c(0,0), 
    legend.justification = c(0,0),
    legend.background = element_blank(), 
    panel.background = element_blank(), 
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(colour="black", face="bold", size=14), 
    axis.title = element_text(size=16, face="bold"), 
    legend.key = element_blank(),
    plot.margin = unit(c(t=1, r=0, b=.75, l=0.5), "cm"),
    panel.border = element_blank()))

  theme_butterfly_r <- list(theme(
    panel.background = element_blank(), 
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "black"), 
    axis.text.x = element_text(colour="black", face="bold", size=14), 
    axis.title.x = element_text(size=16, face="bold"), 
    axis.ticks.y = element_blank(), 
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(), 
    legend.key = element_blank(),
    plot.margin = unit(c(t=1, r=0.5, b=.75, l=0.0), "cm"),
    panel.border = element_blank()))

  if(!is.null(x$fit$family)) nms = unique(names(sort(c(x$pos.weights, x$neg.weights), decreasing = FALSE)))
  
  #vpl <- grid::viewport(width=0.525, height=1, x=0, y=0, just=c("left", "bottom"))
  #vpr <- grid::viewport(width=0.475, height=1, x=0.525, y=0, just=c("left", "bottom"))
  if(!x$bootstrap){
    if(!is.null(x$fit$family)){
      # glm
      poscolwt = 1-x$pos.psi/(x$pos.psi - x$neg.psi)
      if(length(x$pos.weights)==0) x$pos.weights = x$neg.weights*0
      if(length(x$neg.weights)==0) x$neg.weights = x$pos.weights*0
      pright <- ggplot() + 
        stat_identity(aes(x=v, y=w), position = "identity", geom="bar", 
                      data=data.frame(w=x$pos.weights, v=names(x$pos.weights)),
                      fill=gray(poscolwt)) + 
        scale_y_continuous(name="Positive weights", expand=c(0.000,0.000), breaks=c(0.25, 0.5, 0.75)) +
        scale_x_discrete(limits=nms, breaks=nms, labels=nms, drop=FALSE, position="top") +
        geom_hline(aes(yintercept=0)) + 
        coord_flip(ylim=c(0,1)) + 
        theme_butterfly_r
      pleft <- ggplot() + 
        stat_identity(aes(x=v, y=w), position = "identity", geom="bar", 
                      data=data.frame(w=x$neg.weights, v=names(x$neg.weights)),
                      fill=gray(1-poscolwt)) + 
        scale_y_reverse(name="Negative weights", expand=c(0.000,0.000), breaks=c(0.25, 0.5, 0.75)) +
        scale_x_discrete(name="Variable", limits=nms, breaks=nms, labels=nms, drop=FALSE) +
        geom_hline(aes(yintercept=0)) + 
        coord_flip(ylim=c(0,1)) + 
        theme_butterfly_l
      if((length(x$neg.weights)>0 & length(x$pos.weights)>0)){
        maxstr = max(mapply(nchar, c(names(x$neg.weights), names(x$pos.weights))))
        lw = 1+maxstr/20
        p1 <- gridExtra::arrangeGrob(grobs=list(pleft, pright), ncol=2, padding=0.0, widths=c(lw,1))
        if(!suppressprint) {
          grid::grid.newpage()
          grid::grid.draw(p1)
        }
        if(suppressprint) return(p1)
      }
    }
    if(is.null(x$fit$family)){
      # zero inflated
      p1 = list()
      for(modtype in names(x$psi)){
        nms = unique(names(sort(c(x$pos.weights[[modtype]], x$neg.weights[[modtype]]), decreasing = FALSE)))
        poscolwt = 1-x$pos.psi[[modtype]]/(x$pos.psi[[modtype]] - x$neg.psi[[modtype]])
        if(length(x$pos.weights[[modtype]])==0) x$pos.weights[[modtype]] = x$neg.weights[[modtype]]*0
        if(length(x$neg.weights[[modtype]])==0) x$neg.weights[[modtype]] = x$pos.weights[[modtype]]*0
        pright <- ggplot() + 
          stat_identity(aes(x=v, y=w), position = "identity", geom="bar", 
                        data=data.frame(w=x$pos.weights[[modtype]], v=names(x$pos.weights[[modtype]])),
                        fill=gray(poscolwt)) + 
          scale_y_continuous(name="Positive weights", expand=c(0.000,0.000), breaks=c(0.25, 0.5, 0.75)) +
          scale_x_discrete(limits=nms, breaks=nms, labels=nms, drop=FALSE, position="top") +
          geom_hline(aes(yintercept=0)) + 
          coord_flip(ylim=c(0,1)) + 
          theme_butterfly_r
        pleft <- ggplot() + 
          stat_identity(aes(x=v, y=w), position = "identity", geom="bar", 
                        data=data.frame(w=x$neg.weights[[modtype]], v=names(x$neg.weights[[modtype]])),
                        fill=gray(1-poscolwt)) + 
          scale_y_reverse(name="Negative weights", expand=c(0.000,0.000), breaks=c(0.25, 0.5, 0.75)) +
          scale_x_discrete(name=paste0("Variable (", modtype, " model)"), limits=nms, breaks=nms, labels=nms, drop=FALSE) +
          geom_hline(aes(yintercept=0)) + 
          coord_flip(ylim=c(0,1)) + 
          theme_butterfly_l
        if((length(x$neg.weights[[modtype]])>0 & length(x$pos.weights[[modtype]])>0)){
          maxstr = max(mapply(nchar, c(names(x$neg.weights[[modtype]]), names(x$pos.weights[[modtype]]))))
          lw = 1+maxstr/20
          p1[[modtype]] <- gridExtra::arrangeGrob(grobs=list(pleft, pright), ncol=2, padding=0.0, widths=c(lw,1))
        }
      }
      if(!suppressprint) {
        plfun <- function(plt){ 
          grid::grid.newpage()
          grid::grid.draw(plt)
        }
        lapply(p1, plfun)
      }
      if(suppressprint) return(p1)
    }
  }
  if(x$bootstrap){
       # variance based on delta method and knowledge that non-linear
       #functions will always be polynomials in qgcomp
       # default plot for bootstrap results (no weights obtained)
    surv <- NULL # appease R CMD check
    p <- ggplot() 
    if(x$msmfit$family$family=='cox'){
      requireNamespace("survival")
      #construction("warning", "Plot type may change in future releases.")
      rootdat <- as.data.frame(x$fit$x)
      psidat <- data.frame(psi=0)
      rootfun <- function(idx, df){
        df[,x$expnms] <- idx
        df
      }
      rootfun2 <- function(idx, df){
        df[,"psi"] <- idx
        df[,"psi1"] <- idx
        df[,"psi2"] <- idx^2
        df[,"psi3"] <- idx^3
        df[,"psi4"] <- idx^4
        df
      }
      newmarg = lapply(0:(x$q-1), rootfun2, df=psidat)
      margdf = data.frame(do.call("rbind", newmarg))
      newcond = lapply(0:(x$q-1), rootfun, df=rootdat)
      conddf = data.frame(do.call("rbind", newcond))
      msmobj = survfit(x$msmfit, newdata=margdf)
      gcompobj = survfit(x$fit, newdata=conddf)
      mdf = with(msmobj, data.frame(time=time, surv=apply(surv, 1, mean)))
      gdf = with(gcompobj, data.frame(time=time, surv=apply(surv, 1, mean)))
      mdf0 = with(survfit(x$msmfit, newdata=newmarg[[1]]), 
                  data.frame(time=time, surv=surv))
      gdf0 = with(survfit(x$fit, newdata=newcond[[1]]), 
                  data.frame(time=time, surv=apply(surv, 1, mean)))
      mdfx = with(survfit(x$msmfit, newdata=newmarg[[x$q]]), 
                  data.frame(time=time, surv=surv))
      gdfx = with(survfit(x$fit, newdata=newcond[[x$q]]), 
                  data.frame(time=time, surv=apply(surv, 1, mean)))
      p <- p +
        geom_step(aes(x=time, y=surv, color="MSM", linetype="Average (all quantiles)"), data=mdf)+
        geom_step(aes(x=time, y=surv, color="Conditional", linetype="Average (all quantiles)"), data=gdf) + 
        geom_step(aes(x=time, y=surv, color="MSM", linetype="Lowest quantile"), data=mdf0)+
        geom_step(aes(x=time, y=surv, color="Conditional", linetype="Lowest quantile"), data=gdf0) + 
        geom_step(aes(x=time, y=surv, color="MSM", linetype="Highest quantile"), data=mdfx)+
        geom_step(aes(x=time, y=surv, color="Conditional", linetype="Highest quantile"), data=gdfx) + 
        scale_y_continuous(name="Survival", limits=c(0,1)) + 
        scale_x_continuous(name="Time") +
        scale_linetype_discrete(name="")+
        theme(legend.position = c(0.01, 0.01), legend.justification = c(0,0))
    }
    if(x$msmfit$family$family=='gaussian'){
      p <- p + labs(x = "Joint exposure quantile", y = "Y")
       #confidence band
       y = x$y.expectedmsm
       COV = x$covmat.psi
       intvals = as.numeric(names(table(x$index)))
       grad = grad.poly(intvals, x$degree)
       py = tapply(x$y.expectedmsm, x$index, mean)
       varpy = 0*py
       for(i in 1:length(intvals)){
         varpy[i] = se_comb(expnms = paste0('psi', 1:x$degree), 
                        covmat=COV, grad = grad[i,])
       }
       pyup = py + qnorm(.975)*sqrt(varpy)
       pydo = py + qnorm(.025)*sqrt(varpy)
       if(modelband){
         p <- p + geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax, 
                                  fill="Model confidence band"),
                              data=data.frame(ymin=pydo, ymax=pyup, x=intvals/max(intvals)))
       }
       if(flexfit){
         p <- p + geom_smooth(aes(x=x,y=y, color="Smooth conditional fit"),se = FALSE,
                              method = "gam", formula = y ~ s(x, k=length(table(x$index))-1, bs = "cs"),
                              data=data.frame(y=x$y.expected, x=x$index/max(x$index)))
       }
       if(modelfitline){
         p <- p + geom_line(aes(x=x,y=y, color="MSM fit"),
                            data=data.frame(y=y, x=x$index/max(x$index)))
       }
      if(pointwisebars){
        # pairwise comparisons with referent
        ycovmat = x$cov.yhat # bootstrap covariance matrix of E(y|x) from MSM
        pw.diff = c(0,diff(py))
        pw.vars = numeric(length(pw.diff))
        pw.vars[pointwiseref] = 0
        pw.idx = (1:length(py))[-pointwiseref]
        for(j in pw.idx){
          pw.vars[j] = sum(ycovmat[c(pointwiseref,j),c(pointwiseref,j)])  
        }
        pw.up = py + qnorm(.975)*sqrt(pw.vars)
        pw.lo = py + qnorm(.025)*sqrt(pw.vars)
        p <- p + geom_errorbar(aes(x=x,ymin=ymin,ymax=ymax, color="Pointwise 95% CI"), width = 0.05,
                                  data=data.frame(ymin=pw.lo, ymax=pw.up, x=intvals/max(intvals)))
      }
     }
    if(x$msmfit$family$family=='binomial'){
       y = x$y.expectedmsm # probabilities (not on model scale)
       #variance/gradient on model scales
       COV = x$covmat.psi
       intvals = as.numeric(names(table(x$index)))
       grad = grad.poly(intvals, x$degree)
       py = tapply(x$y.expectedmsm, x$index, mean)
       varpy = 0*py   
       for(i in 1:length(intvals)){
         varpy[i] = se_comb(expnms = paste0('psi', 1:x$degree), 
                        covmat=COV, grad = grad[i,])
       }

       if(x$msmfit$family$link=='log'){
         p <- p + labs(x = "Joint exposure quantile", y = "Pr(Y=1)")
         pyup = pmin(exp(log(py) + qnorm(.975)*sqrt(varpy)), 1)
         pydo = pmax(exp(log(py) + qnorm(.025)*sqrt(varpy)), 0)       
         #pyup = exp(log(py) + qnorm(.975)*sqrt(varpy))
         #pydo = exp(log(py) + qnorm(.025)*sqrt(varpy))
       }
       if(x$msmfit$family$link=='logit'){
         p <- p + labs(x = "Joint exposure quantile", y = "Odds(Y=1)")
         pyup = 1/(1+exp(-(log(py/(1-py)) + qnorm(.975)*sqrt(varpy))))
         pydo = 1/(1+exp(-(log(py/(1-py)) + qnorm(.025)*sqrt(varpy))))      
       }

       if(modelband){
         p <- p + geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax, 
                                  fill="Model confidence band"),
                              data=data.frame(ymin=pydo, ymax=pyup, x=intvals/max(intvals)))
       }
       if(flexfit){
         p <- p + geom_smooth(aes(x=x,y=y, color="Smooth conditional fit"),se = FALSE,
                              method = "gam", formula = y ~ s(x, k=length(table(x$index))-1, bs = "cs"),
                              data=data.frame(y=x$y.expected, x=x$index/max(x$index)))
       }
       if(modelfitline){
         p <- p + geom_line(aes(x=x,y=y, color="Model fit"),
                            data=data.frame(y=y, x=x$index/max(x$index)))
       }
       if(pointwisebars){
         # pairwise comparisons with referent
         ycovmat = x$cov.yhat # bootstrap covariance matrix of E(mu|x) from MSM
         pw.diff = c(0,diff(py))
         pw.vars = numeric(length(pw.diff))
         pw.vars[pointwiseref] = 0
         pw.idx = (1:length(py))[-pointwiseref]
         for(j in pw.idx){
           pw.vars[j] = sum(ycovmat[c(pointwiseref,j),c(pointwiseref,j)])  
         }
         if(x$msmfit$family$link=='log'){
           pw.up = pmin(exp(log(py) + qnorm(.975)*sqrt(pw.vars)), 1)
           pw.lo = pmax(exp(log(py) + qnorm(.025)*sqrt(pw.vars)), 0)       
         }
         if(x$msmfit$family$link=='logit'){
           pw.up = 1/(1+exp(-(log(py/(1-py)) + qnorm(.975)*sqrt(pw.vars))))
           pw.lo = 1/(1+exp(-(log(py/(1-py)) + qnorm(.025)*sqrt(pw.vars))))      
         }
         p <- p + geom_errorbar(aes(x=x,ymin=ymin,ymax=ymax, color="Pointwise 95% CI"), width = 0.05,
                                data=data.frame(ymin=pw.lo, ymax=pw.up, x=intvals/max(intvals)))
       }
       
     }
    if(x$msmfit$family$family!='cox'){
     # p <- p + geom_smooth(aes(x=x,y=y, color="Smooth fit"),
     #                     data=data.frame(y=x$y.expected, x=x$index/max(x$index)), 
     #                     method = 'gam', 
     #                     formula=y~s(x, k=4,fx=TRUE), se = FALSE) + 
     #scale_x_continuous(name=("Joint exposure quantile")) + 
     #scale_y_continuous(name="E(outcome)") 
    }
    p <- p + scale_fill_grey(name="", start=.9) + 
      scale_colour_grey(name="", start=0.0, end=0.6) + 
      theme_classic()
    if(!suppressprint) print(p)
  }
  if(suppressprint) return(p)
  #grid.text("Density", x=0.55, y=0.1, gp=gpar(fontsize=14, fontface="bold", fontfamily="Helvetica"))
}

predict.qgcompfit <- function(object, expnms=NULL, newdata=NULL, type="response", ...){
  #' @title predict.qgcompfit: default prediction method for a qgcompfit object (non-survival 
  #' outcomes only)
  #'
  #' @description get predicted values from a qgcompfit object, or make predictions
  #' in a new set of data based on the qgcomfit object. Note that when making predictions
  #' from an object from qgcomp.boot, the predictions are made from the g-computation
  #' model rather than the marginal structural model. Predictions from the marginal
  #' structural model can be obtained via \code{\link[qgcomp]{msm.predict}}
  #' 
  #' @param object "qgcompfit" object from `qgcomp.noboot` or  `qgcomp.boot` functions
  #' @param expnms character vector of exposures of interest
  #' @param newdata (optional) new set of data with all predictors from "qgcompfit" object
  #' @param type  (from predict.glm) the type of prediction required. The default 
  #' is on the scale of the linear predictors; the alternative "response" is on 
  #' the scale of the response  variable. Thus for a default binomial model the 
  #' default predictions are of log-odds (probabilities on logit scale) and 
  #' type = "response" gives the predicted probabilities. The "terms" option 
  #' returns a matrix giving the fitted values of each term in the model formula 
  #' on the linear predictor scale.
  #' @param ... arguments to predict.glm
  #' @import grDevices
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' obj1 <- qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' obj2 <- qgcomp.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=10, seed=125)
  #' set.seed(52)
  #' dat2 <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' summary(predict(obj1, expnms = c('x1', 'x2'), newdata=dat2))
  #' summary(predict(obj2, expnms = c('x1', 'x2'), newdata=dat2))
  if(is.null(newdata)){
   pred <- predict(object$fit, type=type, ...) 
  }
 if(!is.null(newdata)){
   newqdata <- quantize(newdata, expnms, q=NULL, object$breaks)$data
   pred <- predict(object$fit, newdata=newqdata, type=type, ...) 
 }
  return(pred)
}

msm.predict <- function(object, newdata=NULL){
  #' @title msm.predict: secondary prediction method for the MSM within a qgcompfit 
  #' object (non-survival outcomes only). 
  #' 
  #' @description Makes predictions from the MSM (rather than the g-computation 
  #' model) from a "qgcompfit" object. Generally, this should not be used in 
  #' favor of the default \code{\link[qgcomp]{predict.qgcompfit}} function. This
  #' function can only be used following the `qgcomp.boot` function. For the 
  #' `qgcomp.noboot` function, \code{\link[qgcomp]{predict.qgcompfit}} gives 
  #' identical inference to predicting from an MSM.
  #'
  #' @description get predicted values from a qgcompfit object from
  #' \code{\link[qgcomp]{qgcomp.boot}}
  #' 
  #' @param object "qgcompfit" object from `qgcomp.boot` function
  #' @param newdata (optional) new set of data (data frame) with a variable 
  #' called `psi` representing the joint exposure level of all exposures
  #' under consideration
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' obj <- qgcomp.boot(y ~ z + x1 + x2 + I(z*x1), expnms = c('x1', 'x2'), data=dat, q=4, B=10, seed=125)
  #' dat2 <- data.frame(psi=seq(1,4, by=0.1))
  #' summary(msm.predict(obj))
  #' summary(msm.predict(obj, newdata=dat2))
 if(!object$bootstrap) stop("only valid for results from qgcomp.boot function")
 if(is.null(newdata)){
   pred <- predict(object$msmfit, type='response') 
  }
 if(!is.null(newdata)){
   pred <- predict(object$msmfit, newdata=newdata, type='response') 
 }
  return(pred)
}

