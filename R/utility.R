construction <- function(i="", j=""){
  msg ="This function is not yet fully implemented"
  if(j != "") msg = paste0(msg, ": ", j)
  if(i %in% c("msg", "message", "warning", "wrn", "warn")) warning(msg)
  else stop(msg)
}

# fake cox family function
cox <- function(){
  obj = binomial(link="log")
  obj$family="cox"
  obj
}

# fake zi family function
zi <- function(){
  obj = binomial(link="log")
  obj$family="zi"
  obj
}


#printing of zero-inflated results

printZI <- function(x){
  if(class(x$fit) == "zeroinfl"){
    if(!is.null(x$pos.size$count)) {
      cat("Prob(Y ~ count):\n")
      cat(paste0("Scaled effect size (positive direction, sum of positive coefficients = ", signif(x$pos.size$count, 3) , ")\n"))
      if (length(x$pos.weights$count) > 0) {
        print(x$pos.weights$count, digits = 3)
      } else cat("None\n")
      cat("\n")
    }
    if(!is.null(x$neg.size$count)) {
      cat(paste0("Scaled effect size (negative direction, sum of negative coefficients = ", signif(-x$neg.size$count, 3) , ")\n"))
      if (length(x$neg.weights$count) > 0) {
        print(x$neg.weights$count, digits = 3)
      } else cat("None\n")
      cat("\n")
    }
    if(!is.null(x$pos.size$zero)) {
      cat("Prob(Y ~ zero/count):\n")
      cat(paste0("Scaled effect size (positive direction, sum of positive coefficients = ", signif(x$pos.size$zero, 3) , ")\n"))
      if (length(x$pos.weights$zero) > 0) {
        print(x$pos.weights$zero, digits = 3)
      } else cat("None\n")
      cat("\n")
    }
    if(!is.null(x$neg.size$zero)) {
      cat(paste0("Scaled effect size (negative direction, sum of negative coefficients = ", signif(-x$neg.size$zero, 3) , ")\n"))
      if (length(x$neg.weights$zero) > 0) {
        print(x$neg.weights$zero, digits = 3)
      } else cat("None\n")
      cat("\n")
    }
    
    if(x$fit$dist %in% c("poisson", "negbin")){
      estimand <- 'OR/RR'
      cat(paste0("Mixture log(",estimand,")", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    }
    #if(x$fit$dist=="gaussian"){
    #  estimand <- 'log(OR)/mean diff'
    #  cat(paste0("Mixture ",estimand,"", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    #}
    testtype = "Z"
    
    pdat <- list()
    for(modtype in names(x$psi)){
      pdat[[modtype]] <- cbind(Estimate=coef(x)[[modtype]], "Std. Error"=sqrt(x$var.coef[[modtype]]), 
                               "Lower CI"=x$ci.coef[[modtype]][,1], "Upper CI"=x$ci.coef[[modtype]][,2], 
                               "test"=x$zstat[[modtype]], "Pr(>|z|)"=x$pval[[modtype]])
      colnames(pdat[[modtype]])[5] = eval(paste(testtype, "value"))
      numpsi = length(x$psi[[modtype]])
      if(numpsi>0) rnm = c("(Intercept)", c(paste0('psi',1:max(1, numpsi))))
      if(numpsi==0) rnm = c("(Intercept)")
      rownames(pdat[[modtype]]) <- rnm
      cat(paste0("Prob(Y ~ ", modtype,"):\n"))
      printCoefmat(pdat[[modtype]],has.Pvalue=TRUE,tst.ind=5L,signif.stars=FALSE, cs.ind=1L:2)
    }
  }
}

