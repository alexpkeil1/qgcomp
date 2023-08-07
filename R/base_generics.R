# generics from other packages

coef.qgcompfit <- function(object, ...){
  #' @importFrom stats coef
  #' @export
  object$coef
}


df.residual.qgcompfit <- function(object, ...){
  #' @importFrom stats df.residual
  #' @export
  object$fit$df.residual
}


vcov.qgcompfit <- function(object, ...){
  #' @importFrom stats vcov
  #' @export
  object$covmat.coef
}

AIC.qgcompfit <- function(object, ...){
  #' @importFrom stats AIC
  #' @export
  AIC(object$fit)
}

BIC.qgcompfit <- function(object, ...){
  #' @importFrom stats BIC
  #' @export
  BIC(object$fit)
}

logLik.qgcompfit <- function(object, ...){
  #' @importFrom stats logLik
  #' @export
  logLik(object$fit)
}

anova.qgcompfit <- function(object, ...){
  #' @importFrom stats confint
  #' @export
  anova(object$fit)
}

confint.qgcompfit <- function(object, ...){
  #' @importFrom stats anova
  #' @export
  message("not yet implemented")
  anova(object$fit)
}




print.qgcompfit <- function(x, showweights=TRUE, ...){
  #' @title Default printing method for a qgcompfit object
  #' 
  #' @description Gives variable output depending on whether `qgcomp.glm.noboot` or `qgcomp.glm.boot`
  #' is called. For `qgcomp.glm.noboot` will output final estimate of joint exposure
  #' effect (similar to the 'index' effect in weighted quantile sums), as well
  #' as estimates of the 'weights' (standardized coefficients). For `qgcomp.glm.boot`,
  #' the marginal effect is given, but no weights are reported since this approach
  #' generally incorporates non-linear models with interaction terms among exposures,
  #' which preclude weights with any useful interpretation.
  #' 
  #' @param x "qgcompfit" object from `qgcomp`, `qgcomp.glm.noboot` or `qgcomp.glm.boot` 
  #' function
  #' @param showweights logical: should weights be printed, if estimated?
  #' @param ... unused
  #' @seealso \code{\link[qgcomp]{qgcomp.glm.noboot}}, \code{\link[qgcomp]{qgcomp.glm.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @concept variance mixtures
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' obj1 <- qgcomp.glm.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' obj2 <- qgcomp.glm.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=10, seed=125)
  #' # does not need to be explicitly called, but included here for clarity
  #' print(obj1)
  #' print(obj2)
  fam <- x$fit$family$family
  rnm =  c(paste0('psi',1:max(1, length(coef(x)))))
  if(x$hasintercept){
    rnm = c("(Intercept)",rnm[-length(rnm)])
  }
  if(inherits(x, "ziqgcompfit")){
    printZI(x, showweights=showweights, ...)
    return(invisible(x))
  }
  if(!is.null(x$pos.size) & showweights) {
    cat(paste0("Scaled effect size (positive direction, sum of positive coefficients = ", signif(x$pos.size, 3) , ")\n"))
    if (length(x$pos.weights) > 0) {
      print(x$pos.weights, digits = 3)
    } else cat("None\n")
    cat("\n")
  }
  if(!is.null(x$neg.size) & showweights) {
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
  }
  if (fam == "poisson"){
    #message("Poisson family still experimental: use with caution")
    estimand <- 'RR'
    cat(paste0("Mixture log(",estimand,")", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "Z"
  }
  if (fam == "gaussian"){
    cat(paste0("Mixture slope parameters", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "t"
    x$zstat = x$tstat
  }
  if (fam %in% c("cox", "cch")){
    cat(paste0("Mixture log(hazard ratio)", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "Z"
  }
  if (!(fam %in% c("poisson", "binomial", "cox", "cch", "gaussian"))){
    testtype = "Z"
    rnm = c(paste0('psi',1:max(1, length(coef(x)))))
    warning(paste0("The ", fam, " distribution has not been tested with qgcomp! Please use with extreme caution
                   and check results thoroughly with simulated data to ensure it works."))
  }
  plab = ifelse(testtype=="Z", "Pr(>|z|)", "Pr(>|t|)")
  if(is.null(dim(x$ci.coef))){
    pdat <- cbind(Estimate=coef(x), "Std. Error"=sqrt(x$var.coef), "Lower CI"=x$ci.coef[1], "Upper CI"=x$ci.coef[2], "test"=x$zstat, "pval"=x$pval)
  } else{
    pdat <- cbind(Estimate=coef(x), "Std. Error"=sqrt(x$var.coef), "Lower CI"=x$ci.coef[,1], "Upper CI"=x$ci.coef[,2], "test"=x$zstat, "pval"=x$pval)
  }
  colnames(pdat)[which(colnames(pdat)=="test")] = eval(paste(testtype, "value"))
  colnames(pdat)[which(colnames(pdat)=="pval")] = eval(paste(plab))
  rownames(pdat) <- rnm
  printCoefmat(pdat,has.Pvalue=TRUE,tst.ind=5L,signif.stars=FALSE, cs.ind=1L:2)
  invisible(x)
}

summary.qgcompfit <- function(object, ...){
  #' @export
  fam <- object$fit$family$family
  if(inherits(object, "ziqgcompfit")){
    res = summaryZI(object)
    return(res)
  }
  if (fam == "binomial"){
    estimand <- 'OR'
    if(object$bootstrap && object$msmfit$family$link=='log') estimand = 'RR'
    cat(paste0("Mixture log(",estimand,")", ifelse(object$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "Z"
    rnm = c("(Intercept)", c(paste0('psi',1:max(1, length(coef(object))-1))))
  }
  if (fam == "poisson"){
    #message("Poisson family still experimental: use with caution")
    estimand <- 'RR'
    cat(paste0("Mixture log(",estimand,")", ifelse(object$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "Z"
    rnm = c("(Intercept)", c(paste0('psi',1:max(1, length(coef(object))-1))))
  }
  if (fam == "gaussian"){
    cat(paste0("Mixture slope parameters", ifelse(object$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "t"
    rnm = c("(Intercept)", c(paste0('psi',1:max(1, length(coef(object))-1))))
  }
  if (fam %in% c("cox", "cch")){
    cat(paste0("Mixture log(hazard ratio)", ifelse(object$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "Z"
    rnm = c(paste0('psi',1:max(1, length(coef(object)))))
  }
  if (!(fam %in% c("poisson", "binomial", "cox", "cch", "gaussian"))){
    # in development: quasipoisson, Gamma, quasi, quasibinomial, inverse.gaussian
    testtype = "Z"
    rnm = c(paste0('psi',1:max(1, length(coef(object)))))
    warning(paste0("The ", fam, " distribution has not been tested with qgcomp! Please use with extreme caution
                   and check results thoroughly with simulated data to ensure it works."))
  }
  plab = ifelse(testtype=="Z", "Pr(>|z|)", "Pr(>|t|)")
  if(is.null(dim(object$ci.coef))){
    pdat <- cbind(Estimate=coef(object), "Std. Error"=sqrt(object$var.coef), "Lower CI"=object$ci.coef[1], "Upper CI"=object$ci.coef[2], "test"=object$zstat, "pval"=object$pval)
  } else{
    pdat <- cbind(Estimate=coef(object), "Std. Error"=sqrt(object$var.coef), "Lower CI"=object$ci.coef[,1], "Upper CI"=object$ci.coef[,2], "test"=object$zstat, "pval"=object$pval)
  }
  colnames(pdat)[which(colnames(pdat)=="test")] = eval(paste(testtype, "value"))
  colnames(pdat)[which(colnames(pdat)=="pval")] = eval(paste(plab))
  rownames(pdat) <- rnm
  list(coefficients=pdat)
}

family.qgcompfit <- function(object, ...){
  #' @export
  object$fit$family
}



predict.qgcompfit <- function(object, expnms=NULL, newdata=NULL, type="response", ...){
  #' @title Default prediction method for a qgcompfit object (non-survival 
  #' outcomes only)
  #'
  #' @description get predicted values from a qgcompfit object, or make predictions
  #' in a new set of data based on the qgcompfit object. Note that when making predictions
  #' from an object from qgcomp.glm.boot, the predictions are made from the (conditional) g-computation
  #' model rather than the marginal structural model. Predictions from the marginal
  #' structural model can be obtained via \code{\link[qgcomp]{msm.predict}}. Note
  #' that this function accepts non-quantized exposures in "newdata" and automatically
  #' quantizes them according to the quantile cutpoints in the original fit.
  #' 
  #' @param object "qgcompfit" object from `qgcomp.glm.noboot`, `qgcomp.glm.boot`, `qgcomp.zi.noboot`, 
  #' or `qgcomp.zi.boot`functions
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
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' obj1 <- qgcomp.glm.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' obj2 <- qgcomp.glm.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=10, seed=125)
  #' set.seed(52)
  #' dat2 <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' summary(predict(obj1, expnms = c('x1', 'x2'), newdata=dat2))
  #' summary(predict(obj2, expnms = c('x1', 'x2'), newdata=dat2))
  if(is.null(newdata)){
    pred <- predict(object$fit, type=type, ...) 
  }
  if(!is.null(newdata)){
    if(is.null(expnms[1])) expnms = object$expnms # testing
    newqdata <- quantize(newdata, expnms, q=NULL, object$breaks)$data
    pred <- predict(object$fit, newdata=newqdata, type=type, ...) 
  }
  return(pred)
}



