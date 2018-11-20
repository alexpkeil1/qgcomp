
se_comb <- function(expnms, covmat){
  #calculate standard error of weighted linear combination of random variables
  # given a vector of weights and a covariance matrix (not exported)
  # usage: qgcomp:::se_comb(expnms='x', covmat=summary(lmfit)$cov.scaled)
  #calculate standard error of weighted linear combination of random variables
  weightvec = rep(0, dim(covmat)[1])
  weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] = 1
  wcovmat <- weightvec %*% t(weightvec) * covmat
  var <- sum(wcovmat)
  sqrt(var)
}

quantize <- function (data, expnms, q) {
  #' @title create variables representing indicator functions with cutpoints defined
  #' by quantiles
  #' @details This function vectorized version of `quantile_f` from the `gWQS` package
  #' @param data a data frame
  #' @param expnms a character vector with the names of  the columns to be
  #' quantized
  #' @param q integer, number of quantiles used in creating quantized variables
  #' @keywords variance, mixtures
  #' @import stats
  #' @export
  #' @examples
  #' set.seed(1232)
  #' dat = data.frame(y=runif(100), x1=runif(100), x2=runif(100), z=runif(100))
  #' qdata = quantize(data=dat, expnms=c("x1", "x2"), q=4)
  #' table(qdata$x1)
  #' table(qdata$x2)
  #' summary(dat[c('y', 'z')]);summary(qdata[c('y', 'z')]) # not touched
    qt <- function(i){
      # not exported
        datmat <- as.numeric(unlist(data[, expnms[i]]))
        cut(datmat, breaks = unique(quantile(datmat,
             probs = seq(0, 1, by = 1 / q), na.rm = TRUE)), labels = FALSE,
             include.lowest = TRUE) - 1
    }
    data[, expnms] = sapply(1:length(expnms), qt)
    return(data)
}


msm.fit <- function(f, qdata, q, expnms, rr=TRUE, main=TRUE, ...){
  #' @title fitting marginal structural model (MSM) based on g-computation with
  #' quantized exposures
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
  #' @param f an r function representing the conditional model for the outcome, given all
  #' exposures and covariates
  #' @param data a data frame
  #' @param expnms a character vector with the names of  the columns to be
  #' quantized
  #' @param q integer, number of quantiles used in creating quantized variables
  #' @param rr logical, estimate log(risk ratio) (family='binary' only)
  #' @param main logical, internal use: produce estimates of exposure effect (gamma)
  #'  and expected outcomes under g-computation and the MSM
  #' @keywords variance, mixtures
  #' @import stats
    # conditional outcome regression fit
    fit <- glm(f, data = qdata, ...)
    if(fit$family$family=="gaussian") rr=FALSE
    ### 
    # get predictions (set exposure to 1,2,...,q)
    if(is.null(q)){
      q = length(table(qdata[expnms[1]]))
    }
    predit <- function(idx){
      newdata <- qdata
      newdata[,expnms] <- idx-1
      suppressWarnings(predict(fit, newdata=newdata, type='response'))
    }
    predmat = lapply(1:(q+1), predit)
    # fit MSM using g-computation estimates of expected outcomes under joint 
    #  intervention
    nobs <- dim(qdata)[1]
    msmdat <- data.frame(
      Ya = unlist(predmat),
      gamma = rep(0:q, each=nobs))
    if(!rr) suppressWarnings(msmfit <- glm(Ya ~ gamma, data=msmdat,...))
    if(rr)  suppressWarnings(msmfit <- glm(Ya ~ gamma, data=msmdat, family=binomial(link='log')))
    res = list(fit=fit, msmfit=msmfit)
    if(main) {
      res$Ya = msmdat$Ya   # expected outcome under joint exposure, by gcomp
      res$Yamsm = predict(msmfit, type='response')
      res$A =  msmdat$gamma # joint exposure (1 = all exposures set to first quantile)
    }
    res
}


qgcomp.noboot <- function(f, data, expnms=NULL, q=4, alpha=0.05, ...){
  #' @title estimation of quantile g-computation fit (continuous outcome)
  #'  or conditional quantile odds ratio (binary outcome)
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
  #' @param q number of quantiles used to create quantile indicator variables
  #' representing the exposure variables
  #' @param alpha alpha level for confidence limit calculation
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[gqcomp]{qgcomp.boot}}, and \code{\link[gqcomp]{qgcomp}}
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (gamma) and associated variance (var.gamma), as well
  #'  as information on the model fit (fit) and information on the 
  #'  weights/standardized coefficients in the positive (pweights) and 
  #'  negative (nweight) directions.
  #' @keywords variance, mixtures
  #' @import stats
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat = data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
    if (is.null(expnms)) {
      cat("Including all model terms as exposures of interest")
      expnms <- attr(terms(f, data = data), "term.labels")
    }
    if (!is.null(q)){
      qdata <- quantize(data, expnms, q)
    } else qdata <- data
    fit <- glm(f, data = qdata, ...)
    mod <- summary(fit)
    estb <- sum(mod$coefficients[expnms,1])
    seb <- se_comb(expnms, covmat = mod$cov.scaled)
    tstat <- estb / seb
    df <- mod$df.null - length(expnms)
    pval <- 2 - 2 * pt(abs(tstat), df = df)
    pvalz <- 2 - 2 * pnorm(abs(tstat))
    ci <- c(estb + seb * qnorm(alpha / 2), estb + seb * qnorm(1 - alpha / 2))
    # 'weights'
    wcoef <- fit$coefficients[expnms]
    names(wcoef) <- gsub("_q", "", names(wcoef))
    poscoef <- which(wcoef > 0)
    negcoef <- which(wcoef <= 0)
    pweights <- abs(wcoef[poscoef]) / sum(abs(wcoef[poscoef]))
    nweights <- abs(wcoef[negcoef]) / sum(abs(wcoef[negcoef]))
    # 'post-hoc' positive and negative estimators 
    # similar to constrained gWQS
    pos.gamma <- sum(wcoef[poscoef])
    neg.gamma <- sum(wcoef[negcoef])
    nmpos = names(pweights)
    nmneg = names(nweights)
    se.pos.gamma <- se_comb(nmpos, covmat = mod$cov.scaled)
    se.neg.gamma <- se_comb(nmneg, covmat = mod$cov.scaled)
    qx <- qdata[, expnms]
    names(qx) <- paste0(names(qx), "_q")
    res <- list(
      qx = qx, fit = fit, gamma = estb, var.gamma = seb ^ 2, ci = ci,
      expnms=expnms, q=q,
      pos.gamma = pos.gamma, var.pos.gamma = se.pos.gamma^2,
      neg.gamma = neg.gamma, var.neg.gamma = se.neg.gamma^2,
                pweights = sort(pweights, decreasing = TRUE),
                nweights = sort(nweights, decreasing = TRUE), 
                psize = sum(abs(wcoef[poscoef])),
                nsize = sum(abs(wcoef[negcoef])),
                bootstrap=FALSE
                )
      if(fit$family$family=='gaussian'){
        res$tstat = tstat
        res$df = df
        res$pval = pval
      }
      if(fit$family$family=='binomial'){
        res$zstat = tstat
        res$pval = pvalz
      }
    attr(res, "class") <- "qgcompfit"
    res
}


qgcomp.boot <- function(f, data, expnms=NULL, q=4, alpha=0.05, B=200, rr=TRUE, ...){
  #' @title estimation of quantile g-computation fit, using bootstrap confidence intervals
  #'  
  #' @description This function yields population average effect estimates for both continuous and binary outcomes
  #'  
  #' @details Estimates correspond to the average expected change in the
  #'  (log) outcome per quantile increase in the joint exposure to all exposures 
  #'  in `expnms'. Test statistics and confidence intervals are based on 
  #'  a non-parametric bootstrap, using the standard deviation of the bootstrap
  #'  estimates to estimate the standard error. The bootstrap standard error is 
  #'  then used to estimate Wald-type confidence intervals.
  #'
  #' @param f R style formula
  #' @param data data frame
  #' @param expnms character vector of exposures of interest
  #' @param q number of quantiles used to create quantile indicator variables
  #' representing the exposure variables
  #' @param alpha alpha level for confidence limit calculation
  #' @param B integer: number of bootstrap iterations
  #' @param rr logical: if using binary outcome and rr=TRUE, qgcomp.boot will estimate risk ratio rather than odds ratio
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[gqcomp]{qgcomp.noboot}}, and \code{\link[gqcomp]{qgcomp}}
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (gamma) and associated variance (var.gamma), as well
  #'  as information on the model fit (fit) and information on the 
  #'  marginal structural model (msmfit) used to estimate the final effect
  #'  estimates.
  #' @keywords variance, mixtures
  #' @import stats
  #' @export
  #' @examples
  #' set.seed(30)
  #' # continuous outcome
  #' dat = data.frame(y=rnorm(100), x1=runif(100), x2=runif(100), z=runif(100))
  #' #Conditional linear slope
  #' qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4, family=gaussian())
  #' #Marginal linear slope (population average slope, for a purely linear, additive model this will equal the conditional)
  #' qgcomp.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4, family=gaussian(), B=10) #increase B to at least 200 in actual examples
  #' #Population average mixture slope which accounts for non-linearity and interactions
  #' qgcomp.boot(y ~ z + x1 + x2 + I(x1^2) + I(x2*x1), family="gaussian", expnms = c('x1', 'x2'), data=dat, q=4, rr=TRUE, B=10)
  #' # binary outcome
  #' dat = data.frame(y=rbinom(50,1,0.5), x1=runif(50), x2=runif(50), z=runif(50))
  #' #Conditional mixture OR
  #' qgcomp.noboot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), data=dat, q=2)
  #' #Marginal mixture OR (population average OR - in general, this will not equal the conditional mixture OR due to non-collapsibility of the OR)
  #' qgcomp.boot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), data=dat, q=2, B=10)
  #' #Population average mixture RR
  #' qgcomp.boot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), data=dat, q=2, rr=TRUE, B=10)
    # character names of exposure mixture components
    if (is.null(expnms)) {
      cat("Including all model terms as exposures of interest")
      expnms <- attr(terms(f, data = data), "term.labels")
    }
    if (!is.null(q)){
      qdata <- quantize(data, expnms, q)
    } else qdata <- data
    ###
    msmfit <- msm.fit(f, qdata, q, expnms, rr, main=TRUE, ...)
    # main estimate  
    estb <- msmfit$msmfit$coefficients['gamma']
    #bootstrap to get std. error
    nobs = dim(qdata)[1]
    gamma.only <- function(i=1, f=f, qdata=qdata, q=q, expnms=expnms, rr=rr, nobs=nobs, ...){
      set.seed(i)
      msm.fit(f, qdata[sample(1:nobs, nobs, replace = TRUE),], q, expnms, rr, ...)$msmfit$coefficients['gamma']
    }
    bootsamps = sapply(X=1:B, FUN=gamma.only,f=f, qdata=qdata, q=q, expnms=expnms, nobs=nobs)
    seb <- sd(bootsamps)
    tstat <- estb / seb
    df <- nobs - length(attr(terms(f, data = data), "term.labels")) - 2 # df based on obs - gcomp terms - msm terms
    pval <- 2 - 2 * pt(abs(tstat), df = df)
    pvalz <- 2 - 2 * pnorm(abs(tstat))
    ci <- c(estb + seb * qnorm(alpha / 2), estb + seb * qnorm(1 - alpha / 2))
    # 'weights' not applicable in this setting, generally (i.e. if using this function for non-linearity, 
    #   then weights will vary with level of exposure)
    qx <- qdata[, expnms]
    res <- list(
      qx = qx, fit = msmfit$fit, msmfit = msmfit$msmfit, gamma = estb, var.gamma = seb ^ 2, ci = ci,
      expnms=expnms, q=q,
      pos.gamma = NULL, var.pos.gamma = NULL,neg.gamma = NULL, var.neg.gamma = NULL,
      pweights = NULL,nweights = NULL, psize = NULL,nsize = NULL, bootstrap=TRUE,
      y.expected=msmfit$Ya, index=msmfit$A
    )
      if(msmfit$fit$family$family=='gaussian'){
        res$tstat = tstat
        res$df = df
        res$pval = pval
      }
      if(msmfit$fit$family$family=='binomial'){
        res$zstat = tstat
        res$pval = pvalz
      }
    attr(res, "class") <- "qgcompfit"
    res
}


qgcomp <- function(f,data=data,family=gaussian(),rr=TRUE,...){
  #' @title estimation of quantile g-computation fit
  #' 
  #' @description This function automatically selects between qgcomp.noboot and qgcomp.boot
  #'  to select the most efficient approach to estimate the average expected 
  #'  change in the (log) outcome per quantile increase in the joint 
  #'  exposure to all exposures in `expnms'
  #'
  #' @param f R style formula
  #' @param data data frame
  #' @param family `gaussian()` or `binomial()`
  #' @param rr logical: if using binary outcome and rr=TRUE, qgcomp.boot will 
  #' estimate risk ratio rather than odds ratio. Note, to get population average 
  #' effect estimates for a binary outcome, set rr=TRUE (default: ORs are generally not
  #' of interest as population average effects, so if rr=FALSE then a conditional
  #' OR will be estimated, which cannot be interpreted as a population average
  #' effect
  #' @param ... arguments to qgcomp.noboot or qgcomp.boot (e.g. q)
  #' @seealso \code{\link[gqcomp]{qgcomp.noboot}} and \code{\link[gqcomp]{qgcomp.boot}}
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (gamma) and associated variance (var.gamma), as well
  #'  as information on the model fit (fit) and possibly information on the 
  #'  marginal structural model (msmfit) used to estimate the final effect
  #'  estimates (qgcomp.boot only). If appropriate, weights are also reported.
  #' @keywords variance, mixtures
  #' @import stats
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat = data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' qgcomp.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=100)
  #' # automatically selects appropriate method
  #' qgcomp(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' # note for binary outcome this will 
  #' dat = data.frame(y=rbinom(100, 1, 0.5), x1=runif(50), x2=runif(50), z=runif(50))
  #' qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=binomial())
  #' qgcomp.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=100, family=binomial())
  #' # automatically selects appropriate method
  #' qgcomp(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=binomial())
  #' qgcomp(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=binomial(), rr=TRUE)
  terms = attr(terms(f,data=data), 'term.labels')
  doboot = ifelse(isTRUE(grep("I\\(", terms)>0), TRUE, FALSE)
  if(rr | doboot){
    res = qgcomp.boot(f=f,data=data,family=family,rr=rr,...)
  }else res = qgcomp.noboot(f=f,data=data,family=family,...)
  res
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
  #' @seealso \code{\link[gqcomp]{qgcomp.noboot}}, \code{\link[gqcomp]{qgcomp.boot}}, and \code{\link[gqcomp]{qgcomp}}
  #' @keywords variance, mixtures
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat = data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' obj1 = qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' obj2 = qgcomp.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=100)
  #' # does not need to be explicitly called, but included here for clarity
  #' print(obj1)
  #' print(obj2)
  fam <- x$fit$family$family
  if(!is.null(x$psize)) {
    cat(paste0("Scaled effect size (positive direction, sum of positive coefficients = ", signif(x$psize, 3) , ")\n"))
    if (length(x$pweights) > 0) {
      print(x$pweights, digits = 3)
    } else cat("None\n")
    cat("\n")
  }
  if(!is.null(x$nsize)) {
    cat(paste0("Scaled effect size (negative direction, sum of negative coefficients = ", signif(-x$nsize, 3) , ")\n"))
    if (length(x$nweights) > 0) {
      print(x$nweights, digits = 3)
    } else cat("None\n")
    cat("\n")
  }
  if (fam == "binomial"){
    estimand = 'OR'
    if(x$bootstrap && x$msmfit$family$link=='log') estimand = 'RR'
    cat(paste0("Mixture log(",estimand,")", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n"))
    cat(paste0("gamma (CI): ", signif(x$gamma, 3), " (",
             signif(x$ci[1], 3), ",", signif(x$ci[2], 3), "), z=",
             signif(x$zstat, 3), ", p=",
             signif(x$pval, 3), "\n"))
  }
  if (fam == "gaussian"){
    cat(paste0("Mixture slope", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n"))
    cat(paste0("gamma (CI): ", signif(x$gamma, 3), " (",
             signif(x$ci[1], 3), ",", signif(x$ci[2], 3), "), t=",
             signif(x$tstat, 3), ", df=", x$df, ", p=",
             signif(x$pval, 3), "\n"))
  }
}


plot.qgcompfit <- function(x, ...){
  #' @title plot.qgcompfit: default plotting method for a qgcompfit object
  #'
  #' @description Plot a quantile g-computation object. For qgcomp.noboot, this function will
  #' createa butterfly plot of weights. For qgcomp.boot, this function will create
  #' a box plot with smoothed line overlaying that represents a non-parametric
  #' fit of a model to the expected outcomes in the population at each quantile
  #' of the joint exposures (e.g. '1' represents 'at the first quantile for
  #' every exposure')
  #' 
  #' @param x "qgcompfit" object from `qgcomp.noboot` or  `qgcomp.boot` functions
  #' @param ... unused
  #' @seealso \code{\link[gqcomp]{qgcomp.noboot}}, \code{\link[gqcomp]{qgcomp.boot}}, and \code{\link[gqcomp]{qgcomp}}
  #' @import ggplot2 grid gridExtra
  #' @export
  #' @examples
  #' set.seed(12)
  #' dat = data.frame(y=runif(100), x1=runif(100), x2=runif(100), z=runif(100))
  #' ft = qgcomp.noboot(y ~ z + x1 + x2, expnms=c('x1','x2'), data=dat, q=8)
  #' plot(ft)
  #' ft2 = qgcomp.boot(y ~ z + x1 + x2, expnms=c('x1','x2'), data=dat, q=8)
  #' plot(ft2)

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
    axis.text = element_text(colour="black", face="bold", size=14, family="Helvetica"), 
    axis.title = element_text(size=16, face="bold", family="Helvetica"), 
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
    axis.text.x = element_text(colour="black", face="bold", size=14, family="Helvetica"), 
    axis.title.x = element_text(size=16, face="bold", family="Helvetica"), 
    axis.ticks.y = element_blank(), 
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(), 
    legend.key = element_blank(),
    plot.margin = unit(c(t=1, r=0.5, b=.75, l=0.0), "cm"),
    panel.border = element_blank()))

  nms = names(sort(c(x$pweights, x$nweights), decreasing = FALSE))
  
  #vpl <- grid::viewport(width=0.525, height=1, x=0, y=0, just=c("left", "bottom"))
  #vpr <- grid::viewport(width=0.475, height=1, x=0.525, y=0, just=c("left", "bottom"))
  if(!x$bootstrap){
    pright <- ggplot() + 
    stat_identity(aes(x=v, y=w), position = "identity", geom="bar", data=data.frame(w=x$pweights, v=names(x$pweights))) + 
    scale_y_continuous(name="Positive weights", expand=c(0.000,0.000), breaks=c(0.25, 0.5, 0.75)) +
    scale_x_discrete(limits=nms, breaks=nms, labels=nms, drop=FALSE, position="top") +
    geom_hline(aes(yintercept=0)) + 
    coord_flip(ylim=c(0,1)) + 
    theme_butterfly_r
    pleft <- ggplot() + 
    stat_identity(aes(x=v, y=w), position = "identity", geom="bar", data=data.frame(w=x$nweights, v=names(x$nweights))) + 
    scale_y_reverse(name="Negative weights", expand=c(0.000,0.000), breaks=c(0.25, 0.5, 0.75)) +
    scale_x_discrete(name="Variable", limits=nms, breaks=nms, labels=nms, drop=FALSE) +
    geom_hline(aes(yintercept=0)) + 
    coord_flip(ylim=c(0,1)) + 
    theme_butterfly_l
    if(length(x$pweights)==0) print(pleft)
    if(length(x$nweights)==0) print(pright)
    if((length(x$nweights)>0 & length(x$pweights)>0)){
      maxstr = max(mapply(nchar, c(names(x$nweights), names(x$pweights))))
      lw = 1+maxstr/20
      p1 <- gridExtra::arrangeGrob(grobs=list(pleft, pright), ncol=2, padding=0.0, widths=c(lw,1))
      grid::grid.newpage()
      grid::grid.draw(p1)
    }
  }
  if(x$bootstrap){
   # default plot for bootstrap results (no weights obtained)
   p <- ggplot() 
     if(x$msmfit$family$family=='gaussian'){
       #prediction interval (large sample estimator under normal assumption)
       resvar = summary(x$fit)$dispersion
       gammavar = x$var.gamma
       yup = x$y.expected + qnorm(.975)*sqrt(resvar+gammavar)
       ydo = x$y.expected + qnorm(.025)*sqrt(resvar+gammavar)
       p <- p + geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax), data=data.frame(ymin=ydo, ymax=yup, x=x$index)) 
     }
     if(x$msmfit$family$family=='binomial'){
       resvar = x$y.expected*(1-x$y.expected)
       gammavar = x$var.gamma
       yup = x$y.expected + qnorm(.975)*sqrt(resvar+gammavar)
       ydo = x$y.expected + qnorm(.025)*sqrt(resvar+gammavar)
       p <- p + geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax), data=data.frame(ymin=ydo, ymax=yup, x=x$index)) 
     }
     p <- p + geom_smooth(aes(x=x,y=y),data=data.frame(y=x$y.expected, x=x$index), method = 'gam', formula=y~s(x, k=4,fx=TRUE)) + 
     scale_x_continuous(name=("Joint exposure quantile")) + 
     scale_y_continuous(name="E(outcome)") + 
     theme_classic()
   print(p)
  }
  #grid.text("Density", x=0.55, y=0.1, gp=gpar(fontsize=14, fontface="bold", fontfamily="Helvetica"))
}

# todo: predict methods
