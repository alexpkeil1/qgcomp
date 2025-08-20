





.devinstall <- function(...){
  .qgc.require("devtools")
  devtools::install_github("alexpkeil1/qgcomp",...)
}


.fold_list <- function (fold, set1, set2) {
  fold_list <- list(fold = fold, set1 = set1, set2 = set2)
  fold_list
}

.xfitfold_list_from_foldvec <- function (fold, folds, ordermat) {
  nfolds <- length(unique(folds))
  set1 <- which(folds %in% ordermat[1:(nfolds - 1), fold])
  set2 <- which(folds == ordermat[nfolds, fold])
  .fold_list(fold, set1, set2)
}

.make_xfitfolds_iid <- function (n, V = 5) {
  folds <- rep(seq_len(V), length = n)
  folds <- sample(folds)
  combinations <- utils::combn(V, V - 1)
  combinations <- rbind(combinations, apply(combinations, 2, 
                                            function(x) setdiff(1:V, x)))
  if (V > 1) 
    foldobj = lapply(1:V, .xfitfold_list_from_foldvec, 
                     folds = folds, ordermat = combinations)
  if (V == 1) 
    foldobj = list(.fold_list(fold = 1, set1 = 1:n, 
                                           set2 = 1:n))
  foldobj
}


.xfit_grab <- function(foldres, stat,initval=0,whichval=1){
  statmat = matrix(initval, nrow=length(foldres), ncol=2)
  for(i in seq_len(length(foldres))){
    res = foldres[[i]]
    if(!is.null(res$pos.fit)){
      statmat[i,1] = unlist(res$pos.fit[stat])[whichval]
    }
    if(!is.null(res$neg.fit)){
      statmat[i,2] = unlist(res$neg.fit[stat])[whichval]
    }
  }
  statmat
}

.xfit_procfolds <- function(foldres){
  list(
    intercepts = .xfit_grab(foldres, "coef"),
    psis = .xfit_grab(foldres, "coef",whichval=2),
    vars_intercept = .xfit_grab(foldres, "var.coef"),
    vars_psi = .xfit_grab(foldres, "var.coef",whichval=2)
  )
}

.xfit_proclist <- function(proclist,n){
  int_ests = apply(proclist$intercepts, 2, median)
  psi_ests = apply(proclist$psis, 2, median)
  int_resids = t(t(proclist$intercepts)-int_ests)
  psi_resids = t(t(proclist$psis)-psi_ests)
  int_vars = proclist$vars_intercept/n + int_resids^2
  psi_vars = proclist$vars_psi/n + psi_resids^2
  c(proclist, list(int_ests=int_ests,
                   psi_ests=psi_ests,
                   int_resids =int_resids ,
                   psi_resids =psi_resids ,
                   int_vars = apply(int_vars, 2, median),
                   psi_vars = apply(psi_vars, 2, median)
                   )
  )
}
  
.calcxfitci <- function(xdf, alpha=0.05){
  xdf =  cbind(xdf, xdf[,1] + qnorm(alpha/2)* xdf[,2])
  xdf =  cbind(xdf, xdf[,1] + qnorm(1-alpha/2)* xdf[,2])
  xdf =  cbind(xdf, xdf[,1]/xdf[,2])
  xdf =  cbind(xdf, 2 - 2 * pnorm(abs(xdf[,1]/xdf[,2])))
  colnames(xdf) <- c("Estimate", "Std. Err", "Lower CI", "Upper CI", "z value", "Pr(>|t|)")
  rownames(xdf) <- c("(Intercept)", "psi")
  xdf
}

.xfit_coeftab <- function(x){
  x$neg.coeftab = (rbind(
    c(x$int_ests[2],sqrt(x$int_vars[2])),
    c(x$psi_ests[2],sqrt(x$psi_vars[2]))
  ))
  x$pos.coeftab = (rbind(
    c(x$int_ests[1], sqrt(x$int_vars[1])),
    c(x$psi_ests[1], sqrt(x$psi_vars[1]))
  ))
  x$neg.coeftab = .calcxfitci(x$neg.coeftab)
  x$pos.coeftab = .calcxfitci(x$pos.coeftab)
  x
}


.qgcomp_xfitpartials <- function(
    data,
    fun = c("qgcomp.glm.noboot", "qgcomp.cox.noboot", "qgcomp.zi.noboot"),
    V=10,
    expnms=NULL,
    .fixbreaks=TRUE,
    ...
){
  n = nrow(data)
  folds <- .make_xfitfolds_iid(n=n, V=V)
  foldres = list()
  for(f in folds){
    traindata = data[f$set1,,drop=FALSE]
    validdata = data[f$set2,,drop=FALSE]
    foldres[[f$fold]] = qgcomp.partials(
      fun=fun,
      expnms = expnms, 
      traindata = traindata, 
      validdata = validdata,
      .fixbreaks = .fixbreaks,
      ...
    )
  }
  res = .xfit_procfolds(foldres)
  res = .xfit_proclist(res, n)
  res = .xfit_coeftab(res)
  alpha = foldres[[1]]$pos.fit$alpha
  res = c(res, list(foldres = foldres, n=n))
  class(res) <-  "qgcompmultixfit"
  res
}


print.qgcompmultixfit <- function(x,...){
  #' @export
  #cat(paste0("\nVariables with positive effect sizes in training data: ", paste(x$posmix, collapse = ", ")))
  #cat(paste0("\nVariables with negative effect sizes in training data: ", paste(x$negmix, collapse = ", ")))
  cat("\nPartial effect sizes estimated using V-fold cross-fitting\n")
  cat("Positive direction \n")
  if(!is.null(x$pos.coeftab)){
    printCoefmat(x$pos.coeftab, has.Pvalue = TRUE)
  } else{
    cat("\n No positive coefficients in model fit to training data ")
  }
  cat("\nNegative direction \n")
  if(!is.null(x$neg.coeftab)){
    printCoefmat(x$neg.coeftab, has.Pvalue = TRUE)
  } else{
    cat("\n No negative coefficients in model fit to training data ")
  } 
  
}


# in progress plotting function
.plot_boot_multinomial_base <- function(r, p, x, modelband, flexfit, modelfitline, pointwisebars, pointwiseref=1, alpha=0.05){
  #
  x$fit
  x$msmfit
  stop("Not yet implemented")
  p <- p #+ #scale_y_log10()
    p <- p + labs(x = "Joint exposure quantile", y = "P(Y=y)") + lims(x=c(0,1))
    if(modelband) p <- p + .plot.or.mod.bounds(x,alpha)
    if(flexfit)   p <- p + .plot.or.smooth.line(x)
    if(modelfitline) p <- p + .plot.logitlin.line(x)
    if(pointwisebars) p <- p + .plot.or.pw.boot(x,alpha,pointwiseref)
  p
}


#
pointwisebound.coxhr.boot <- function (x, pointwiseref = 1, alpha = 0.05) {
  link = x$msmfit$family$link
  pwr = pointwiseref + 0
  designmat = .makenewdesign(x, seq_len(x$q) - 1)
  refrow <- designmat[pointwiseref, ]
  nrows <- nrow(designmat)
  res = as.data.frame(designmat)
  
  lnhr = designmat %*% coef(x)
  se.lnhr = numeric(nrows)
  se.lnhr[pwr] = 0
  pw.idx = seq_len(nrow(res))[-pwr]
  vc <- vcov(x)
  for (nr in seq_len(nrows)) {
    grad <- as.numeric(designmat[nr, ] - refrow)
    se.lnhr[nr] <- se_comb(names(grad), covmat = vc, grad)
  }
  res = .pointwise.log.boot(x$q, lnhr, se.lnhr, alpha, pointwiseref)
  names(res) <- c(names(res)[c(1,2)], "lnhr", "hr", "se.lnhr", "ll.hr", "ul.hr", "ll.lnhr", "ul.lnhr")
  res
}







