#----------------------------------------------------------------#
# Helper functions for qgcomp with multinomial outcome ####
#----------------------------------------------------------------#


.get_baseline_prob_multinom <- function(x){
  cf = coef(x)
  ints = grep("[Ii]ntercept", names(cf))
  baseprob = 1/(sum(exp(cf[ints]))-1)
  baseprob
}


.flatten <- function(xm, nm="value"){
  xx = dimnames(xm)
  stodf = expand.grid(xx[[1]],xx[[2]])
  rownames(stodf) = paste(stodf$Var1, stodf$Var2, sep =".")
  stodf[,nm] = 0
  #for(i in seq_along(nrow(stodf))){
  for(i in seq_len(nrow(stodf))){
      stodf[i,nm] = xm[stodf$Var1[i], stodf$Var2[i]]
  }
  t(as.matrix(stodf[,nm,drop=FALSE]))[1,]
}


.psiest_qgcomp_multi <- function(
    ufit, 
    expnms
){
  ucoef = coef(ufit)
  psi = apply(ucoef[,expnms],1,sum)
  psi
}

.partialpsiest_qgcomp_multi <- function(
    ufit, 
    expnms
){
  ucoef = coef(ufit)
  possum = function(x){sum(x[x>0])}
  negsum = function(x){sum(x[x<0])}
  pospsi = apply(ucoef[,expnms],1,possum)
  negpsi = apply(ucoef[,expnms],1,negsum)
  list(positive_psi=pospsi,negative_psi=negpsi)
}


.vcov_qgcomp_multi <- function(
    ufit, 
    expnms
){

  labs = ufit$lab
  reflabs = labs[2:length(labs)]
  nlevels = length(reflabs)
  uvcov = vcov(ufit)
  psi_vcov = matrix(NA, nlevels, nlevels, dimnames=list(reflabs, reflabs))
  #psi_vcov = 
  for(i in 1:(nlevels)){
    inames = paste0(reflabs[i],":", expnms)
    weightvec <- rep(0, dim(uvcov)[1])
    weightvec[which(colnames(as.matrix(uvcov)) %in% inames)] <- 1
    psi_vcov[i, i] <- weightvec %*% uvcov %*% weightvec
    for(j in 1:(i-1)){
      jnames = paste0(reflabs[j],":", expnms)
      acol = which(colnames(as.matrix(uvcov)) %in% inames)
      bcol = which(colnames(as.matrix(uvcov)) %in% jnames)
      psi_vcov[i, j] <- psi_vcov[j, i] <- sum(uvcov[acol, bcol])
    }
  }
  psi_vcov
}

.vcov_qgcomp_multi_coef <- function(
    ufit, 
    expnms
){
  
  labs = ufit$lab
  reflabs = labs[2:length(labs)]
  intlabs = paste0(reflabs,":(Intercept)")
  psilabs = paste0(reflabs,":psi")
  matlabs = c(intlabs, psilabs)
  nlevels = length(reflabs)
  nentries = length(matlabs)
  uvcov = vcov(ufit)
  coef_vcov = matrix(NA, nentries, nentries, dimnames=list(matlabs, matlabs))
  psi_vcov = matrix(NA, nlevels, nlevels, dimnames=list(reflabs, reflabs))
  #psi_vcov = 
  inames = paste0(reflabs,":(Intercept)")
  psinames = paste0(reflabs,":psi")
  for(i in 1:(nlevels)){
    xnames = paste0(reflabs[i],":", expnms)
    iname = inames[i]
    psiname = psinames[i]
    weightvec <- rep(0, dim(uvcov)[1])
    weightvec[which(colnames(as.matrix(uvcov)) %in% xnames)] <- 1
    
    coef_vcov[c(iname, psiname),c(iname, psiname)] = vc_comb(aname=iname, xnames, uvcov, grad=weightvec) 
    #psi_vcov[i, i] <- weightvec %*% uvcov %*% weightvec
    j = 1
    while(j<i){
      jname = inames[j]
      # intercepts
      coef_vcov[iname,jname] <- coef_vcov[jname,iname] <- uvcov[iname,jname]
      xjnames = paste0(reflabs[j],":", expnms)
      psinamej = psinames[j]
      weightvec <- weightvecj <- rep(0, dim(uvcov)[1])
      weightvecj[which(colnames(as.matrix(uvcov)) %in% xjnames)] <- 1
      weightvec[which(colnames(as.matrix(uvcov)) %in% xnames)] <- 1
      # psi coefficient covariance with intercepts
      coef_vcov[c(iname, psinamej),c(iname, psinamej)] = vc_comb(aname=iname, xjnames, uvcov, grad=weightvecj) 
      coef_vcov[c(jname, psiname),c(jname, psiname)] = vc_comb(aname=jname, xnames, uvcov, grad=weightvec) 
      # psi coefficient covariance
      acol = which(colnames(as.matrix(uvcov)) %in% xnames)
      bcol = which(colnames(as.matrix(uvcov)) %in% xjnames)
      coef_vcov[psiname, psinamej] <- coef_vcov[psinamej, psiname] <- sum(uvcov[acol, bcol])
      j = j+1
    }
  }
  coef_vcov
}


.calc_qgcomp_weights <- function(
    ufit, 
    expnms
){
  labs = ufit$lab
  nlevels = length(ufit$lab)-1
  pos.weights = list()
  neg.weights = list()
  weights = matrix(NA, ncol = length(expnms), nrow=nlevels)
  #dimnames(weights)[2] = labs[2:length(labs)]
  
  for(i in 1:nlevels){
    #inames = paste0(i,":", expnms)
    icoef <- coef(ufit)[i,expnms]
    pos.coef <- which(icoef > 0)
    neg.coef <- which(icoef <= 0)
    pos.weights <- abs(icoef[pos.coef])/sum(abs(icoef[pos.coef]))
    neg.weights <- abs(icoef[neg.coef])/sum(abs(icoef[neg.coef]))
    weights[i,] <- c(pos.weights, -neg.weights)[expnms]
  }
  rownames(weights) = labs[2:length(labs)]
  colnames(weights) = expnms
  weights
}

#----------------------------------------------------------------#
# functions for testing hypotheses ####
#----------------------------------------------------------------#

#' Hypothesis testing about a joint effect of exposures on a multinomial outcome
#'
#' @description Tests the null hypothesis that the joint effect of the mixture components is homogenous across all referent outcome types
#' 
#' @param x Result from qgcomp multinomial fit (qgcompmultfit object).
#' @param ... Unused
#' @export
homogeneity_test <- function(x,...){
  UseMethod("homogeneity_test")
}

#' Hypothesis testing about a joint effect of exposures on a multinomial outcome
#' 
#' @description Tests the null hypothesis that the joint effect of the mixture components is null across all referent outcome types (Test of global null effect of the mixture on a quantized basis)
#'
#' @param x Result from qgcomp multinomial fit (qgcompmultfit object).
#' @param ... Unused
#' @returns qgcompmulttest object (list) with results of a chi-squared test
#' @export
joint_test <- function(x,...){
  UseMethod("joint_test")
}

#' @importFrom utils combn
#' @rdname homogeneity_test
#' @export
homogeneity_test.qgcompmultfit <- function(x,...){
  nullh = paste0("H_0: ", paste0("psi_", x$labs[2:length(x$labs)], collapse=" = "))
  df = x$nlevels-1
  whichcont = combn(x$nlevels,2)
  L = matrix(0, nrow=df, ncol=x$nlevels)
  for(j in 1:df){
    L[j,whichcont[,j]] = c(-1,1)
  }
  B = x$psi
  V = x$covmat.psi
  lbc = L %*% B# - C
  lvl = L %*% V %*% t(L)
  X2 = as.numeric(t(lbc) %*% solve(lvl) %*% lbc)
  p = pchisq(X2, df, lower.tail = FALSE)
  res = list(nullh=nullh,labs=x$labs,B=B,V=V,L=L,chisq=X2,df=df,pvalue=p)
  attr(res, "class") <- c("qgcompmulttest", "list")
  res
}



#' @importFrom utils combn
#' @rdname joint_test
#' @export
joint_test.qgcompmultfit <- function(x,...){
  nullh = paste0("H_0:", paste0(" psi_", x$labs[2:length(x$labs)], " = 0", collapse=" &"))
  df = x$nlevels
  whichcont = combn(df,1)
  L = matrix(0, nrow=df, ncol=df)
  for(j in 1:(df)){
    L[j,whichcont[,j]] = 1
  }
  B = x$psi
  V = x$covmat.psi
  lbc = L %*% B# - C
  lvl = L %*% V %*% t(L)
  X2 = as.numeric(t(lbc) %*% solve(lvl) %*% lbc)
  p = pchisq(X2, df, lower.tail = FALSE)
  res = list(nullh=nullh,labs=x$labs,B=B,V=V,L=L,chisq=X2,df=df,pvalue=p)
  attr(res, "class") <- c("qgcompmulttest", "list")
  res
}


#----------------------------------------------------------------#
# generic printing/summary functions ####
#----------------------------------------------------------------#



#' @exportS3Method base::print
print.qgcompmultfit <- function(x, ...){
  if(!x$bootstrap){
    cat("Weights\n")
    print(x$weights)
    cat("Partial effects (positive)\n")
    print(x$partial_psi$positive_psi)
    cat("Partial effects (negative)\n")
    print(x$partial_psi$negative_psi)
  }
  if(!x$bootstrap)
    cat("\nMixture slope parameters (Standard CI):\n")
  if(x$bootstrap)
    cat("\nMixture slope parameters (Bootstrap CI):\n")
  #print(qgcompobj$coeftable)
  coeftable = cbind(Estimate=x$coef, 
                    `Std. Error`=sqrt(x$var.coef), 
                    `Lower CI`=x$ci.coef[,1], 
                    `Upper CI`=x$ci.coef[,2], 
                    `Z value`=x$Z, 
                    `Pr(>|Z|)` = x$pvalues) 
  printCoefmat(coeftable)
}

#' Summarize gcompmultfit object
#' @description Summary printing to include coefficients, standard errors, hypothesis tests, weights
#'
#' @param object Result from qgcomp multinomial fit (qgcompmultfit object).
#' @param ... Unused
#' @param tests Character vector (e.g. c("global", "homogeneity")) that determine the types of hypothesis tests that are printed
#' @returns qgcompmulttest object (list) with results of a chi-squared test
#' @export
#' @exportS3Method base::summary
#' @rdname summary
summary.qgcompmultfit <- function(object, ..., tests=NULL){
  # to do
  cat("Reference outcome levels:\n")
  cat(object$labs)
  if(!object$bootstrap){
    cat("\nWeights\n")
    print(object$weights)
    cat("\nSum of positive coefficients \n")
    print(object$partial_psi$positive_psi)
    cat("Sum of negative coefficients \n")
    print(object$partial_psi$negative_psi)
  }
  
  if(!object$bootstrap)
    cat("\nMixture slope parameters (Standard CI):\n")
  if(object$bootstrap)
    cat("\nMixture slope parameters (Bootstrap CI):\n")
  coeftable = cbind(Estimate=object$coef, 
                    `Std. Error`=sqrt(object$var.coef), 
                    `Lower CI`=(object$ci[,1]), 
                    `Upper CI`=(object$ci[,2]), 
                    `Z value`=object$Z, 
                    `Pr(>|Z|)` = object$pvalues) 
  printCoefmat(coeftable, P.values=FALSE, has.Pvalue=FALSE, cs.ind=c(1,2), tst.ind=NULL, na.print="")
  #cat("\nStandard errors\n")
  #printCoefmat(object$stderrs, P.values=FALSE, has.Pvalue=FALSE, cs.ind=c(1,2), tst.ind=NULL, na.print="")
  #cat("\nWald Z\n")
  #printCoefmat(object$Z, P.values=FALSE, has.Pvalue=FALSE, tst.ind=c(1,2), na.print="")
  #cat("\n p(>|Z|)\n")
  #printCoefmat(object$pvalues, P.values=TRUE, has.Pvalue=TRUE, na.print="")
  #print(x$pvalues)
  testvals = sapply(tests, function(x) tolower(substr(x,1,1)))
  if("g" %in% testvals){
    cat("\nWald global test, ")
    print(joint_test(object))
    
  }
  if("h" %in% testvals){
    cat("\nWald homogeneity test, ")
    print(homogeneity_test(object))
  }
}

#' @exportS3Method base::print
print.qgcompmulttest <- function(x,...){
  nullh = x$nullh
  cat(paste0(nullh, "\n"))
  cat(paste0("Chi^2 (df=", x$df, ") = ", signif(x$chisq), ", p = ", format.pval(x$pvalue),"\n"))
}

#----------------------------------------------------------------#
# modeling functions ####
#----------------------------------------------------------------#


#' Quantile g-computation for multinomial outcomes
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
#' id/cluster). Note that qgcomp.multinomial.noboot will not produce cluster-appropriate
#' standard errors (this parameter is essentially ignored in qgcomp.multinomial.noboot).
#' qgcomp.qgcomp.multinomial.boot can be used for this, which will use bootstrap
#' sampling of clusters/individuals to estimate cluster-appropriate standard
#' errors via bootstrapping.
#' @param weights "case weights" - passed to the "weight" argument of
#' \code{\link[nnet]{multinom}}
#' @param alpha alpha level for confidence limit calculation
#' @param bayes Logical, Not yet implemented (gives and error if set to TRUE)
#' @param ... arguments to nnet::multinom
#' @seealso \code{\link[qgcomp]{qgcomp.glm.noboot}}, \code{\link[nnet]{multinom}}
#' @family qgcomp_methods
#' @return a qgcompmultfit object, which contains information about the effect
#'  measure of interest (psi) and associated variance (var.psi), as well
#'  as information on the model fit (fit) and information on the
#'  weights/standardized coefficients in the positive and
#'  negative directions (weights).
#' @concept variance mixtures
#' @import stats nnet
#' @export
#' @examples
#' data("metals") # from qgcomp package
#' # create categorical outcome from the existing continuous outcome (usually, one will already exist)
#' metals$ycat = factor(quantize(metals, "y",q=4)$data$y, levels=c("0", "1", "2", "3"), 
#'                      labels=c("cct", "ccg", "aat", "aag")) 
#' # restrict to smaller dataset for simplicity
#' smallmetals = metals[,c("ycat", "arsenic", "lead", "cadmium", "mage35")]
#' 
#' ### 1: Define mixture and underlying model ####
#' mixture = c("arsenic", "lead", "cadmium")
#' f0 = ycat ~ arsenic + lead + cadmium # the multinomial model 
#' # (be sure that factor variables are properly coded ahead of time in the dataset)
#' 
#' rr = qgcomp.multinomial.noboot(
#'  f0, 
#'  expnms = mixture,
#'  q=4, 
#'  data = smallmetals, 
#'  )
#'  
#'  ### 5: Create summary qgcomp object for nice printing ####
#'  
#'  summary(rr, tests=c("H")) # include homogeneity test
#'  
#'  # 95% confidence intervals
#'  confint(rr, level=0.95)
#'  rr$breaks # quantile cutpoints for exposures
#'  # homogeneity_test(rr)
#'  joint_test(rr)
#'
qgcomp.multinomial.noboot <- function(f,
                          data,
                          expnms=NULL,
                          q=4,
                          breaks=NULL,
                          id=NULL,
                          weights,
                          alpha=0.05,
                          bayes=FALSE,
                          ...){
  
  newform <- terms(f, data = data)
  hasintercept = as.logical(attr(newform, "intercept"))
  
  nobs = nrow(data)
  origcall <- thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("formula", "data", "weights", "offset"), names(thecall), 0L)
  #m <- match(c("f", "data", "weights", "offset"), names(thecall), 0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE
  
  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame()) # a model frame pulled in from the environment in which the function was called
  if(hasweights){
    data$weights <- as.vector(model.weights(thecalle))
  } else data$weights = rep(1, nobs)
  ####
  if (is.null(expnms)) {
    
    #expnms <- attr(terms(f, data = data), "term.labels")
    expnms <- attr(newform, "term.labels")
    
    message("Including all model terms as exposures of interest\n")
  }
  lin = checknames(expnms)
  if(!lin) stop("Model appears to be non-linear: use qgcomp.multinomial.boot instead")
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
    qdata$id__ = seq_len(dim(qdata)[1])
  }
  
  if(!bayes) {
    fit = nnet::multinom(newform, data = qdata, 
                                         weights=weights, Hess=TRUE, trace=FALSE,
                                         ...)
  }
  if(bayes){
    stop("bayes=TRUE is not implemented for this model")
    #requireNamespace("arm")
    #fit <- bayesglm(newform, data = qdata,
    #                weights=weights,
    #                ...)
  }
  
  #if(length(setdiff(expnms, rownames(mod$coefficients)))>0){
  #  stop("Model aliasing occurred, likely due to perfectly correlated quantized exposures.
  #         Try one of the following:
  #           1) set 'bayes' to TRUE in the qgcomp function (recommended)
  #           2) set 'q' to a higher value in the qgcomp function (recommended)
  #           3) check correlation matrix of exposures, and drop all but one variable in each highly correlated set  (not recommended)
  #         ")
  #}
  
  psi = .psiest_qgcomp_multi(fit, expnms)
  partpsi = .partialpsiest_qgcomp_multi(fit, expnms)
  #
  psi_vcov = .vcov_qgcomp_multi(fit, expnms)
  coef_vcov = .vcov_qgcomp_multi_coef(fit, expnms)
  #
  qgcweights <- .calc_qgcomp_weights(fit,expnms)
  coeftable = cbind(`(Intercept)`=coef(fit)[,1], psi=psi)
  labs = fit$lab
  reflabs = labs[2:length(fit$lab)]
  stderrs = sqrt(cbind(`(Intercept)`=as.numeric(diag(vcov(fit)[paste0(reflabs, ":(Intercept)"),paste0(reflabs, ":(Intercept)")])), psi=diag(psi_vcov)))
  estb = .flatten(coeftable)
  seb = .flatten(stderrs)
  Z = estb / seb
  psiidx = which(!grepl("[Ii]ntercept", names(estb)))
  ci <- cbind(estb + seb * qnorm(alpha / 2), estb + seb * qnorm(1 - alpha / 2))
  #.qgcomp_object
  #qgcompobj = list(
  qx <- qdata[, expnms]
  names(qx) <- paste0(names(qx), "_q")
  fit$family = multinom_family()
  qgcompobj = .qgcompmult_object(
    qx = qx,
    fit = fit,
    labs = labs,
    nlevels = length(labs)-1,
    psi = psi,
    var.psi = diag(psi_vcov),
    covmat.psi = psi_vcov,
    ci.psi = ci[psiidx,],
    #
    coef = estb,
    var.coef=seb^2,
    covmat.coef = coef_vcov,
    ci.coef = ci,
    zstat = Z,
    pval = pnorm(abs(Z), lower.tail=FALSE)*2,
    #
    partial_psi = partpsi,
    expnms=expnms, q=q, breaks=br, degree=NULL,
    weights=qgcweights,
    alpha=alpha,
    call=origcall,
    hasintercept=hasintercept,
    bootstrap=FALSE
  )
  qgcompobj
}


msm_multinomial_fit <- function(f,
                                qdata,
                                intvals,
                                expnms,
                                main=TRUE,
                                degree=1,
                                id=NULL,
                                weights,
                                bayes=FALSE,
                                MCsize=nrow(qdata), 
                                hasintercept=TRUE, 
                                ...){
  #' @title Fitting marginal structural model (MSM) within quantile g-computation
  #' @description This is an internal function called by \code{\link[qgcomp]{qgcomp.multinomial.boot}},
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
  #' @param main logical, internal use: produce estimates of exposure effect (psi)
  #'  and expected outcomes under g-computation and the MSM
  #' @param degree polynomial bases for marginal model (e.g. degree = 2
  #'  allows that the relationship between the whole exposure mixture and the outcome
  #'  is quadratic. Default=1)
  #' @param id (optional) NULL, or variable name indexing individual units of
  #' observation (only needed if analyzing data with multiple observations per
  #' id/cluster)
  #' @param weights "case weights" - passed to the "weight" argument of
  #' \code{\link[nnet]{multinom}}
  #' @param bayes use underlying Bayesian model (`arm` package defaults). Results
  #' in penalized parameter estimation that can help with very highly correlated
  #' exposures. Note: this does not lead to fully Bayesian inference in general,
  #' so results should be interpreted as frequentist.
  #' @param MCsize integer: sample size for simulation to approximate marginal
  #'  zero inflated model parameters. This can be left small for testing, but should be as large
  #'  as needed to reduce simulation error to an acceptable magnitude (can compare psi coefficients for
  #'  linear fits with qgcomp.zi.noboot to gain some intuition for the level of expected simulation
  #'  error at a given value of MCsize)
  #' @param hasintercept (logical) does the model have an intercept?
  #' @param ... arguments to nnet::multinom
  #' @seealso \code{\link[qgcomp]{qgcomp.glm.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @concept variance mixtures
  #' @import stats arm
  #' @export
  #' @examples
  #' data("metals") # from qgcomp package
  #' # create categorical outcome from the existing continuous outcome (usually, one will already exist)
  #' metals$ycat = factor(quantize(metals, "y",q=4)$data$y, levels=c("0", "1", "2", "3"), 
  #'                      labels=c("cct", "ccg", "aat", "aag")) 
  #' # restrict to smaller dataset for simplicity
  #' smallmetals = metals[,c("ycat", "arsenic", "lead", "cadmium", "mage35")]
  #' 
  #' ### 1: Define mixture and underlying model ####
  #' mixture = c("arsenic", "lead", "cadmium")
  #' f0 = ycat ~ arsenic + lead + cadmium # the multinomial model 
  #' # (be sure that factor variables are properly coded ahead of time in the dataset)
  #' qdat <- quantize(smallmetals, mixture, q=4)$data
  #' mod <- msm_multinomial_fit(f0,
  #'         expnms = mixture, qdata=qdat, intvals=1:4, bayes=FALSE)
  #' summary(mod$fit) # outcome regression model
  #' summary(mod$msmfit) # msm fit (variance not valid - must be obtained via bootstrap)
  
  newform <- terms(f, data = qdata)
  origY <- eval(attr(newform, "variables")[[2]], qdata)
  ytype = "factor"
  if(is.numeric(origY))
    ytype = "numeric"
  nobs = nrow(qdata)
  thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("qdata", "data", names(thecall))
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("formula", "data", "weights", "offset"), names(thecall), 0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE
  
  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame())
  if(hasweights){
    qdata$weights <- as.vector(model.weights(thecalle))
  } else qdata$weights = rep(1, nobs)
  
  if(is.null(id)) {
    id <- "id__"
    qdata$id__ <- seq_len(dim(qdata)[1])
  }
  # conditional outcome regression fit
  nidx = which(!(names(qdata) %in% id))
  if(!bayes) {
    fit = nnet::multinom(newform, data = qdata, 
                         weights=weights, Hess=TRUE, trace=FALSE,model=TRUE,
                         ...)
    fit$model = qdata
  }
  if(bayes){
    stop("bayes=TRUE is not implemented for this model")
    #requireNamespace("arm")
    #fit <- bayesglm(newform, data = qdata,
    #                weights=weights,
    #                ...)
  }
  #if(fit$family$family %in% c("gaussian", "poisson")) rr=FALSE
  ###
  # get predictions (set exposure to 0,1,...,q-1)
  if(is.null(intvals)){
    intvals <- (seq_len(length(table(qdata[expnms[1]])))) - 1
  }
  predit <- function(idx, newdata, ytype=ytype){
    #newdata <- qdata
    newdata[,expnms] <- idx
    classpreds = suppressWarnings(predict(fit, newdata=newdata, type='probs'))
    classlabs = fit$lev
    p1 = classlabs[which(rmultinom(1, 1, classpreds[1,])==1)]
    pvec = rep(p1, nrow(newdata))
    for(j in 2:length(pvec))
      pvec[j] = classlabs[which(rmultinom(1, 1, classpreds[j,])==1)]
    pvec
  }
  if(MCsize==nrow(qdata)){
    newdata <- qdata
  }else{
    newids <- data.frame(temp=sort(sample(unique(qdata[,id, drop=TRUE]), MCsize,
                                          #probs=weights, #bootstrap sampling with weights works with fixed weights, but not time-varying weights
                                          replace = TRUE
    )))
    names(newids) <- id
    newdata <- merge(qdata,newids, by=id, all.x=FALSE, all.y=TRUE)[seq_len(MCsize),]
  }
  predmat = lapply(intvals, predit, newdata=newdata)
  if(ytype=="factor"){
    factfun = function(x, levs=levels(origY)){
      factor(x,levs)
    }
    predmat = lapply(predmat, factfun)
  }
  if(ytype=="numeric"){
    predmat = lapply(predmat, as.numeric)
  }
  # fit MSM using g-computation estimates of expected outcomes under joint
  #  intervention
  #nobs <- dim(qdata)[1]
  
  msmdat <- data.frame(
    #cbind(
    Ya = unlist(predmat),
    psi = rep(intvals, each=MCsize),
    weights = rep(newdata$weights, times=length(intvals))
    #times=length(table(qdata[expnms[1]])))
    #)
  )
  # to do: allow functional form variations for the MSM via specifying the model formula
  hasintercept=TRUE
  msmforms = paste0("Ya ~ ", 
                    ifelse(hasintercept, "1 +", "-1 +"), 
                    "poly(psi, degree=",degree,", raw=TRUE)"
  )
  msmform = as.formula(msmforms)
  
  if(bayes){
    stop("bayes=TRUE is not implemented for this model")
    #
    #suppressWarnings(msmfit <- bayesglm(msmform, data=msmdat,
    #                                            weights=weights, x=TRUE,
    #                                            ...))
  }
  if(!bayes){
    suppressWarnings(msmfit <- nnet::multinom(msmform, data=msmdat,
                                              weights=weights, Hess=TRUE, trace=FALSE, model=TRUE,
                                              ...))
    msmfit$model = msmdat
    idx = ifelse(hasintercept, 2:length(msmfit$coefnames), seq_along(msmfit$coefnames))
    nm = paste0("psi", gsub("[a-zA-Z= (),]+", "", msmfit$coefnames[idx]))
    msmfit$coefnames[idx] <- msmfit$vcoefnames[idx] <- nm
    
  }
  res <- list(fit=fit, msmfit=msmfit)
  if(main) {
    res$Ya <- msmdat$Ya   # expected outcome under joint exposure, by gcomp
    res$Yamsm <- predict(msmfit, type='probs')
    res$Yamsml <- NULL
    res$A <- msmdat$psi # joint exposure (0 = all exposures set category with
    # upper cut-point as first quantile)
  }
  res
}



#' @title Quantile g-computation for multinomial outcomes
#'
#' @description This function estimates a dose-response parameter representing a one quantile
#' increase in a set of exposures of interest. This model estimates the parameters of a marginal
#' structural model (MSM) based on g-computation with quantized exposures. Note: this function
#' allows linear and non-additive effects of individual components of the exposure, as well as
#' non-linear joint effects of the mixture via polynomial basis functions, which increase the
#' computational computational burden due to the need for non-parametric bootstrapping.
#'
#' @details Estimates correspond to the average expected change in the
#'  probability of an outcome type per quantile increase in the joint exposure to all exposures
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
#' id/cluster). Note that qgcomp.multinomial.noboot will not produce cluster-appropriate
#' standard errors. qgcomp.multinomial.boot can be used for this, which will use bootstrap
#' sampling of clusters/individuals to estimate cluster-appropriate standard
#' errors via bootstrapping.
#' @param weights "case weights" - passed to the "weight" argument of
#' \code{\link[nnet]{multinom}}
#' @param alpha alpha level for confidence limit calculation
#' @param B integer: number of bootstrap iterations (this should typically be >=200,
#'  though it is set lower in examples to improve run-time).
#' @param rr (not used)
#' @param degree polynomial bases for marginal model (e.g. degree = 2
#'  allows that the relationship between the whole exposure mixture and the outcome
#'  is quadratic (default = 1).
#' @param seed integer or NULL: random number seed for replicable bootstrap results
#' @param bayes use underlying Bayesian model (`arm` package defaults). Results
#' in penalized parameter estimation that can help with very highly correlated
#' exposures. Note: this does not lead to fully Bayesian inference in general,
#' so results should be interpreted as frequentist.
#' @param MCsize integer: sample size for simulation to approximate marginal
#'  zero inflated model parameters. This can be left small for testing, but should be as large
#'  as needed to reduce simulation error to an acceptable magnitude (can compare psi coefficients for
#'  linear fits with qgcomp.multinomial.noboot to gain some intuition for the level of expected simulation
#'  error at a given value of MCsize). This likely won't matter much in linear models, but may
#'  be important with binary or count outcomes.
#' @param parallel use (safe) parallel processing from the future and future.apply packages
#' @param parplan (logical, default=FALSE) automatically set future::plan to plan(multisession) (and set to existing plan, if any, after bootstrapping)
#' @param ... arguments to nnet::multinom
#' @return a qgcompfit object, which contains information about the effect
#'  measure of interest (psi) and associated variance (var.psi), as well
#'  as information on the model fit (fit) and information on the
#'  marginal structural model (msmfit) used to estimate the final effect
#'  estimates.
#' @concept variance mixtures
#' @import stats
#' @family qgcomp_methods
#' @export
#' @examples
#' data("metals") # from qgcomp package
#' # create categorical outcome from the existing continuous outcome (usually, one will already exist)
#' metals$ycat = factor(quantize(metals, "y",q=4)$data$y, levels=c("0", "1", "2", "3"), 
#'                      labels=c("cct", "ccg", "aat", "aag")) 
#' # restrict to smaller dataset for simplicity
#' smallmetals = metals[,c("ycat", "arsenic", "lead", "cadmium", "mage35")]
#' 
#' ### 1: Define mixture and underlying model ####
#' mixture = c("arsenic", "lead", "cadmium")
#' f0 = ycat ~ arsenic + lead + cadmium # the multinomial model 
#' # (be sure that factor variables are properly coded ahead of time in the dataset)
#' rr = qgcomp.multinomial.boot(
#'  f0, 
#'  expnms = mixture,
#'  q=4, 
#'  data = smallmetals, 
#'  B = 5, # set to higher values in real examples
#'  MCsize = 100,  # set to higher values in small samples
#'  )
#'
#' rr2 = qgcomp.multinomial.noboot(
#'  f0, 
#'  expnms = mixture,
#'  q=4, 
#'  data = smallmetals
#'  )
#'  
#'  ### 5: Create summary qgcomp object for nice printing ####
#'  
#'  summary(rr, tests=c("H")) # include homogeneity test
#'  
#'  # 95% confidence intervals
#'  #confint(rr, level=0.95)
#'  #rr$breaks # quantile cutpoints for exposures
#'  # homogeneity_test(rr)
#'  #joint_test(rr)
#'
#' qdat = simdata_quantized(
#'   outcometype="multinomial", 
#'   n=10000, corr=c(-.9), coef=cbind(c(.2,-.2,0,0), c(.1,.1,.1,.1)), 
#'   q = 4
#' )
#' 
#'  rr_sim = qgcomp.multinomial.noboot(
#'   y~x1+x2+x3+x4, 
#'   expnms = c("x1", "x2", "x3", "x4"),
#'   q=4, 
#'   data = qdat
#'  )
#'  
#'  rr_sim2 = qgcomp.multinomial.boot(
#'   y~x1+x2+x3+x4, 
#'   expnms = c("x1", "x2", "x3", "x4"),
#'   q=4, 
#'   data = qdat,
#'   B=1
#'  )

#'
qgcomp.multinomial.boot <- function(
    f,
    data,
    expnms=NULL,
    q=4,
    breaks=NULL,
    id=NULL,
    weights,
    alpha=0.05,
    B=200,
    rr=TRUE,
    degree=1,
    seed=NULL,
    bayes=FALSE,
    MCsize=nrow(data),
    parallel=FALSE, 
    parplan = FALSE,
    ...
){
  oldq = NULL
  if(is.null(seed)) seed = round(runif(1, min=0, max=1e8))
  
  newform <- terms(f, data = data)
  origY <- eval(attr(newform, "variables")[[2]], data)
  class(origY)
  hasintercept = as.logical(attr(newform, "intercept"))
  class(newform) <- "formula"
  
  nobs = nrow(data)
  origcall <- thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("formula", "data", "weights", "offset"), names(thecall), 0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE
  
  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame())
  if(hasweights){
    data$weights <- as.vector(model.weights(thecalle))
  } else data$weights = rep(1, nobs)
  
  
  if (is.null(expnms)) {
    
    #expnms <- attr(terms(f, data = data), "term.labels")
    expnms <- attr(newform, "term.labels")
    
    message("Including all model terms as exposures of interest\n")
  }
  lin = checknames(expnms)
  if(!lin) stop("Model appears to be non-linear and I'm having trouble parsing it:
                  please use `expnms` parameter to define the variables making up the exposure")
  if (!is.null(q) & !is.null(breaks)){
    # if user specifies breaks, prioritize those
    oldq = q
    q <- NULL
  }
  if (!is.null(q) | !is.null(breaks)){
    ql <- qgcomp::quantize(data, expnms, q, breaks)
    qdata <- ql$data
    br <- ql$breaks
    if(is.null(q)){
      # rare scenario with user specified breaks and q is left at NULL
      nvals <- length(br[[1]])-1
    } else{
      nvals <- q
    }
    intvals <- (seq_len(nvals))-1
  } else {
    # if( is.null(breaks) & is.null(q)) # also includes NA
    qdata <- data
    # if no transformation is made (no quantiles, no breaks given)
    # then draw distribution values from quantiles of all the exposures
    # pooled together
    # : allow user specification of this
    nvals = length(table(unlist(data[,expnms])))
    if(nvals < 10){
      message("\nNote: using all possible values of exposure as the
              intervention values\n")
      p = length(expnms)
      intvals <- as.numeric(names(table(unlist(data[,expnms]))))
      
      br <- lapply(seq_len(p), function(x) c(-1e16, intvals[2:nvals]-1e-16, 1e16))
    }else{
      message("\nNote: using quantiles of all exposures combined in order to set
          proposed intervention values for overall effect (25th, 50th, 75th %ile)
        You can ensure this is valid by scaling all variables in expnms to have similar ranges.")
      intvals = as.numeric(quantile(unlist(data[,expnms]), c(.25, .5, .75)))
      br <- NULL
    }
  }
  if(is.null(id)) {
    id <- "id__"
    qdata$id__ <- seq_len(dim(qdata)[1])
  }
  ###
  msmfit <- msm_multinomial_fit(newform, qdata, intvals, expnms, main=TRUE,degree=degree, id=id,
                                weights,
                                bayes,
                                MCsize=MCsize, hasintercept = hasintercept,
                                ...)
  # main estimate
  #estb <- as.numeric(msmfit$msmfit$coefficients[-1])
  estb <- as.matrix(coef(msmfit$msmfit))
  #bootstrap to get std. error
  nobs <- dim(qdata)[1]
  nids <- length(unique(qdata[,id, drop=TRUE]))
  starttime = Sys.time()
  psi.only <- function(i=1, f=f, qdata=qdata, intvals=intvals, expnms=expnms, degree=degree,
                       nids=nids, id=id,
                       weights,MCsize=MCsize, hasintercept = hasintercept,
                       ytype,
                       ...){
    if(i==2 & !parallel){
      timeiter = as.numeric(Sys.time() - starttime)
      if((timeiter*B/60)>0.5) message(paste0("Expected time to finish: ", round(B*timeiter/60, 2), " minutes \n"))
    }
    bootids <- data.frame(temp=sort(sample(unique(qdata[,id, drop=TRUE]), nids,
                                           replace = TRUE
    )))
    names(bootids) <- id
    qdata_ <- merge(qdata,bootids, by=id, all.x=FALSE, all.y=TRUE)
    ft = msm_multinomial_fit(f, qdata_, intvals, expnms, main=FALSE, degree, id, weights=weights, bayes, MCsize=MCsize,
                             hasintercept = hasintercept,
                             ...)
    levs = msmfit$fit$lev
    ypredmat=predict(ft$msmfit, type="probs")
    sampmulti = function(len, pp) {
      levs[sample.int(len, prob = pp)[1]]
    }
    #
    ypred = apply(ypredmat, 1, sampmulti, len = length(levs))
    if(ytype=="factor"){
      ypred = factor(ypred, levels(origY))
    }
    if(ytype=="numeric"){
      ypred = as.numeric(ypred)
    }
    
    
    yhatty = data.frame(yhat=ypredmat, psi=ft$msmfit$model[,"psi"])
    
    #
    sumfun <- function(ypredmat, psi, psival){
      cm = colMeans(ypredmat[which(psi == psival),])
      cmt = t(c(cm))
      rownames(cmt) <- psival
      cmt
    }
    # the yhat estimates will be identical across individuals due to this being a marginal model
    margpredlist <- lapply(unique(ft$msmfit$model[,"psi"]), function(x) sumfun(ypredmat, ft$msmfit$model[,"psi"], x))
    
    list(
      cf = coef(ft$msmfit),
      margpreds = do.call(rbind, margpredlist)
    )
  }
  set.seed(seed)
  ytype = "factor"
  if(is.numeric(origY))
    ytype= "numeric"
  if(parallel){
    #Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
    if (parplan) {
      oplan <- future::plan(strategy = future::multisession)
      on.exit(future::plan(oplan), add = TRUE)      
    }
    bootsamps <- future.apply::future_lapply(X=seq_len(B), FUN=psi.only,f=f, qdata=qdata, intvals=intvals,
                                             expnms=expnms, degree=degree, nids=nids, id=id,
                                             weights=qdata$weights,MCsize=MCsize, hasintercept = hasintercept,
                                             future.seed=TRUE, ytype=ytype,
                                             ...)
  }else{
    bootsamps <- lapply(X=seq_len(B), FUN=psi.only,f=f, qdata=qdata, intvals=intvals,
                        expnms=expnms, degree=degree, nids=nids, id=id,
                        weights=weights, MCsize=MCsize, hasintercept = hasintercept,
                        ytype=ytype,
                        ...)
  }
  # bootstrap samples of marginal class probabilities
  hats = do.call("rbind", lapply(bootsamps, function(x) .flatten(x$margpreds, "pred")))
  # bootstrap samples of coefficients
  tcoef = do.call("rbind", lapply(bootsamps, function(x) .flatten(x$cf, "coef")))
  estb = .flatten(estb, "coef")
  psiidx = which(!grepl("[Ii]ntercept", names(estb)))
  cov.yhat = cov(hats)
  seb <- apply(tcoef, 2, sd)
  covmat <- cov(tcoef)
  cnms = c(paste0("psi", seq_len(degree)))
  if(hasintercept)
    cnms = c("(intercept)", cnms)
  colnames(covmat) <- rownames(covmat) <- names(estb)# <- cnms
  tstat <- estb / seb
  df <- nobs - length(attr(terms(f, data = data), "term.labels")) - 1 - degree # df based on obs - gcomp terms - msm terms
  pval <- 2 - 2 * pt(abs(tstat), df = df)
  pvalz <- 2 - 2 * pnorm(abs(tstat))
  ci <- cbind(estb + seb * qnorm(alpha / 2), estb + seb * qnorm(1 - alpha / 2))
  # outcome 'weights' not applicable in this setting, generally (i.e. if using this function for non-linearity,
  #   then weights will vary with level of exposure)
  if (!is.null(oldq)){
    q = oldq
  }
  qx <- qdata[, expnms]
  names(qx) <- paste0(names(qx), "_q")
  res <- .qgcomp_object(
    qx = qx, 
    fit = msmfit$fit, 
    msmfit = msmfit$msmfit,
    labs = msmfit$msmfit$lab,
    nlevels = length(msmfit$msmfit$lab)-1,
    psi = estb[psiidx], 
    var.psi = seb[psiidx] ^ 2, 
    covmat.psi=covmat[psiidx,psiidx, drop=FALSE], 
    ci.psi = ci[psiidx,],
    #
    coef = estb, 
    var.coef = seb ^ 2, 
    covmat.coef=covmat, 
    ci.coef = ci,
    zstat = tstat,
    pval = pvalz,
    #
    expnms=expnms, q=q, breaks=br, degree=degree,
    weights=NULL,
    alpha=alpha,
    call=origcall,
    hasintercept=hasintercept,
    bootstrap=TRUE,
    y.expected=msmfit$Ya, y.expectedmsm=msmfit$Yamsm, index=msmfit$A,
    bootsamps = bootsamps,
    cov.yhat=cov.yhat
  )
  class(res) <- c("qgcompmultfit", class(res))
  res
}

