.devinstall <- function(...){
  .qgc.require("devtools")
  devtools::install_github("alexpkeil1/qgcomp",...)
}

# to replace qgcomp partials
.qgcomp.partials <- function(
    fun = c("qgcomp.noboot", "qgcomp.cox.noboot", "qgcomp.zi.noboot"),
    traindata=NULL,
    validdata=NULL,
    expnms=NULL,
    .fixbreaks=FALSE,
    ...
){
  # currently broken
  if(is.null(traindata) | is.null(validdata))
    stop("traindata and validdata must both be specified")
  #
  traincall <- validcall <- match.call(expand.dots = TRUE)
  droppers <- match(c("traindata", "validdata", ".fixbreaks", "fun"), names(traincall), 0L) #index (will need to add names here if more arguments are added)
  traincall[["data"]] <- eval(traincall[["traindata"]], parent.frame())
  validcall[["data"]] <- eval(validcall[["validdata"]], parent.frame())
  traincall <- traincall[-c(droppers)]
  validcall <- validcall[-c(droppers)]
  hasbreaks = ifelse("breaks" %in% names(traincall), TRUE, FALSE)
  if(hasbreaks && .fixbreaks)
    .fixbreaks=FALSE
  #
  if(is.function(fun)){
    traincall[[1L]] <- validcall[[1L]] <- fun
  }else{
    traincall[[1L]] <- validcall[[1L]] <- as.name(fun[1])
  }
  train.fit = eval(traincall, parent.frame())
  #####
  if(.fixbreaks){
    validcall$breaks = train.fit$breaks
    validcall$q = NULL
  }
  ######
  posnms = names(train.fit$pos.weights)
  negnms = names(train.fit$neg.weights)
  if(length(posnms)==1 && all(posnms==c("count", "zero"))){
    posnms = names(train.fit$pos.weights$count)
    negnms = names(train.fit$neg.weights$count)
  }
  res = list(train.fit=train.fit)
  res$negmix <- res$posmix <- "none"
  if(length(posnms)>0){
    res$posmix = posnms
    vc = as.list(validcall)
    vc$expnms = c(posnms)
    res$pos.fit <- eval(as.call(vc), parent.frame())
  }
  if(length(negnms)>0){
    res$negmix = negnms
    vc = as.list(validcall)
    vc$expnms = c(negnms)
    res$neg.fit <- eval(as.call(vc), parent.frame())
    
  }
  class(res) <- "qgcompmultifit"
  res
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
  combinations <- combn(V, V - 1)
  combinations <- rbind(combinations, apply(combinations, 2, 
                                            function(x) setdiff(1:V, x)))
  if (V > 1) 
    foldobj = lapply(1:V, .xfitfold_list_from_foldvec, 
                     folds = folds, ordermat = combinations)
  if (V == 1) 
    foldobj = list(.xfitfold_list_from_foldvec(fold = 1, set1 = 1:n, 
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
  int_resids = t(t(proclist$intercepts)-int_est)
  psi_resids = t(t(proclist$psis)-psi_est)
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
  


.qgcomp_xfitpartials <- function(
    data,
    fun = c("qgcomp.noboot", "qgcomp.cox.noboot", "qgcomp.zi.noboot"),
    V=10,
    expnms=NULL,
    ...
){
  folds <- .make_xfitfolds_iid(n=nrow(data), V=V)
  foldres = list()
  for(f in folds){
    traindata = dat[f$set1,]
    validdata = dat[f$set2,]
    foldres[[f$fold]] = .qgcomp.partials(
      fun=fun,
      expnms = expnms, 
      traindata = traindata, 
      validdata = validdata,
      ...
    )
  }
  res = .xfit_procfolds(foldres)
  res = .xfit_proclist(res, n)
  alpha = foldres[[1]]$pos.fit$alpha
  list(res=res, foldres = foldres, n=nrow(data))
}


