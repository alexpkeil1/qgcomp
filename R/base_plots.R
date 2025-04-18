
######## underlying functions ##########



.plot.md.mod.bounds <- function(x,alpha){
  ymin <- ymax <- v <- w <- y <- NULL
  isboot <- x$bootstrap
  isee <- inherits(x, "eeqgcompfit")
  if(isboot) modbounds = modelbound.boot(x, pwonly = FALSE, alpha = alpha)
  if(isee) modbounds = modelbound.ee(x, pwonly = FALSE, alpha = alpha)

  ret = geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax, 
                        fill="Model confidence band"),
                    data=data.frame(ymin=modbounds$ll.pw, ymax=modbounds$ul.pw, 
                                    x=modbounds$quantile.midpoint))
  list(ribbon=ret)
}


.plot.rr.mod.bounds <- .plot.md.mod.bounds
.plot.or.mod.bounds <- .plot.md.mod.bounds

.plot.linear.smooth.line <- function(x){
  ymin <- ymax <- v <- w <- y <- NULL
  ret = geom_smooth(aes(x=x,y=y, color="Smooth conditional fit"),se = FALSE,
                    method = "gam", formula = y ~ s(x, k=length(table(x$index))-1, bs = "cs"),
                    data=data.frame(y=x$y.expected, x=(x$index+0.5)/max(x$index+1)))
  list(smooth = ret)
}

.plot.rr.smooth.line <- .plot.linear.smooth.line

.plot.or.smooth.line <- function(x){
  ymin <- ymax <- v <- w <- y <- NULL
  ret = geom_smooth(aes(x=x,y=y, color="Smooth conditional fit"),se = FALSE,
                    method = "gam", formula = y ~ s(x, k=length(table(x$index))-1, bs = "cs"),
                    data=data.frame(y=(x$y.expected)/(1-x$y.expected), x=(x$index+0.5)/max(x$index+1)))
  list(smooth=ret)
}



.plot.linear.line <- function(x){
  ymin <- ymax <- v <- w <- y <- NULL
  ret = geom_line(aes(x=x,y=y, color="MSM fit"),
                  data=data.frame(y=x$y.expectedmsm, x=(x$index+0.5)/max(x$index+1)))
  list(line = ret)
}

.plot.loglin.line <- .plot.linear.line

.plot.logitlin.line <- function(x){
  ymin <- ymax <- v <- w <- y <- NULL
  ret = geom_line(aes(x=x,y=y, color="MSM fit"),
                  data=data.frame(y=(x$y.expectedmsm/(1-x$y.expectedmsm)), x=(x$index+0.5)/max(x$index+1)))
  list(line=ret)
}



.plot.md.pw.boot <- function(x, alpha, pointwiseref){
  ymin <- ymax <- v <- w <- y <- NULL
  # plot actual risk or odds bounds
  pwbdat = pointwisebound.boot(x, alpha=alpha, pointwiseref=pointwiseref)
  py = pwbdat$linpred
  ll = pwbdat$ll.linpred
  ul = pwbdat$ul.linpred
  list(
    point = geom_point(aes(x=x,y=y, 
                           color=paste0("Pointwise ",as.character(100*(1-alpha)),"% CI")),
                       data=data.frame(y=py, x=pwbdat$quantile.midpoint)) ,
    errorbar = geom_errorbar(aes(x=x,ymin=ymin,ymax=ymax, 
                                 color=paste0("Pointwise ",as.character(100*(1-alpha)),"% CI")), width = 0.03,
                             data=data.frame(ymin=ll, ymax=ul, x=pwbdat$quantile.midpoint))
  )
}

.plot.rr.pw.boot <- function(x, alpha, pointwiseref){
  ymin <- ymax <- v <- w <- y <- NULL
  # plot actual risk or odds bounds
  pwbdat = pointwisebound.boot(x, alpha=alpha, pointwiseref=pointwiseref)
  py = exp(pwbdat$linpred)
  ll = exp(pwbdat$ll.linpred)
  ul = exp(pwbdat$ul.linpred)
  list(
    point = geom_point(aes(x=x,y=y, 
                           color=paste0("Pointwise ",as.character(100*(1-alpha)),"% CI")),
                       data=data.frame(y=py, x=pwbdat$quantile.midpoint)) ,
    errorbar = geom_errorbar(aes(x=x,ymin=ymin,ymax=ymax, 
                                 color=paste0("Pointwise ",as.character(100*(1-alpha)),"% CI")), width = 0.03,
                             data=data.frame(ymin=ll, ymax=ul, x=pwbdat$quantile.midpoint))
  )
}

.plot.or.pw.boot <- .plot.rr.pw.boot

.plot.zi.pw.boot <- function(x, alpha, pointwiseref){
  ymin <- ymax <- v <- w <- y <- NULL
  # plot actual risk or odds bounds
  pwbdat = pointwisebound.boot(x, alpha=alpha, pointwiseref=pointwiseref)
  py = (pwbdat$ey) # e(Y) 
  ll = py*(pwbdat$ll.rr)
  ul = py*(pwbdat$ul.rr)
  list(
    point = geom_point(aes(x=x,y=y, 
                           color=paste0("Pointwise ",as.character(100*(1-alpha)),"% CI")),
                       data=data.frame(y=py, x=pwbdat$quantile.midpoint)) ,
    errobar = geom_errorbar(aes(x=x,ymin=ymin,ymax=ymax, 
                                color=paste0("Pointwise ",as.character(100*(1-alpha)),"% CI")), width = 0.03,
                            data=data.frame(ymin=ll, ymax=ul, x=pwbdat$quantile.midpoint))
  )
}


.plot.md.pw.noboot <- function(x, alpha, pointwiseref){
  ymin <- ymax <- v <- w <- y <- NULL
  # plot actual risk or odds bounds
  pwbdat = pointwisebound.noboot(x, alpha=alpha, pointwiseref=pointwiseref)
  py = pwbdat$linpred
  ll = pwbdat$ll.linpred
  ul = pwbdat$ul.linpred
  list(
    point = geom_point(aes(x=x,y=y, 
                           color=paste0("Pointwise ",as.character(100*(1-alpha)),"% CI")),
                       data=data.frame(y=py, x=pwbdat$quantile.midpoint)) ,
    errorbar = geom_errorbar(aes(x=x,ymin=ymin,ymax=ymax, 
                                 color=paste0("Pointwise ",as.character(100*(1-alpha)),"% CI")), width = 0.03,
                             data=data.frame(ymin=ll, ymax=ul, x=pwbdat$quantile.midpoint))
  )
}

.plot.rr.pw.noboot <- function(x, alpha, pointwiseref){
  ymin <- ymax <- v <- w <- y <- NULL
  # plot actual risk or odds bounds
  pwbdat = pointwisebound.noboot(x, alpha=alpha, pointwiseref=pointwiseref)
  py = exp(pwbdat$linpred)
  ll = exp(pwbdat$ll.linpred)
  ul = exp(pwbdat$ul.linpred)
  list(
    point = geom_point(aes(x=x,y=y, 
                           color=paste0("Pointwise ",as.character(100*(1-alpha)),"% CI")),
                       data=data.frame(y=py, x=pwbdat$quantile.midpoint)) ,
    errobar = geom_errorbar(aes(x=x,ymin=ymin,ymax=ymax, 
                                color=paste0("Pointwise ",as.character(100*(1-alpha)),"% CI")), width = 0.03,
                            data=data.frame(ymin=ll, ymax=ul, x=pwbdat$quantile.midpoint))
  )
}

.plot.or.pw.noboot <- .plot.rr.pw.noboot




.plfun <- function(plt){ 
  grid::grid.newpage()
  grid::grid.draw(plt) 
}



######## Main plot functions ##########
.plot_noboot_multi_base <- function(r, x, nms, theme_butterfly_r, theme_butterfly_l){
  v <- w <- NULL
  pospsi = x$partial_psi$positive_psi[r]
  negpsi = x$partial_psi$negative_psi[r]
  nms = colnames(x$weights)
  weights = x$weights[r,]
  pos.weights = weights[weights>=0]
  neg.weights = -weights[weights<0]
  varnm = names(pospsi)
  
  poscolwt = 1-pospsi/(pospsi - negpsi)
  #ggplot() +
  pright <-  c(list(
    stat_identity(aes(x=v, y=w), position = "identity", geom="bar", 
                  data=data.frame(w=pos.weights, v=names(pos.weights)),
                  fill=gray(poscolwt)), 
    scale_y_continuous(name="Positive weights", expand=c(0.000,0.000), breaks=c(0.25, 0.5, 0.75)), 
    scale_x_discrete(limits=nms, breaks=nms, labels=nms, drop=FALSE, position="top"),
    geom_hline(aes(yintercept=0)), 
    coord_flip(ylim=c(0,1))),
    theme_butterfly_r
  )
  pleft <-  c(list(
    stat_identity(aes(x=v, y=w), position = "identity", geom="bar", 
                  data=data.frame(w=neg.weights, v=names(neg.weights)),
                  fill=gray(1-poscolwt)),
    scale_y_continuous(trans="reverse", name="Negative weights", expand=c(0.000,0.000), breaks=c(0.25, 0.5, 0.75), limits=c(1,0),labels = waiver()),
    scale_x_discrete(name=varnm, limits=nms, breaks=nms, labels=nms, drop=FALSE),
    geom_hline(aes(yintercept=0)),
    #coord_flip(ylim=c(0,1), expand = FALSE) + 
    coord_flip()), 
    theme_butterfly_l
  )
    
  if((length(neg.weights)>0 || length(pos.weights)>0)){
    #maxstr = max(mapply(nchar, c(names(x$neg.weights), names(x$pos.weights))))
    maxstr = max(nchar(c(names(neg.weights), names(pos.weights))))
    lw = 1+maxstr/20
    pleft = ggplot() + pleft
    pright = ggplot() + pright
    p <- gridExtra::arrangeGrob(grobs=list(pleft, pright), ncol=2, padding=0.0, widths=c(lw,1))
  }
  p
}

.plot_noboot_base <- function(x, nms, theme_butterfly_r, theme_butterfly_l){
  v <- w <- NULL
  # glm
  poscolwt = 1-x$pos.psi/(x$pos.psi - x$neg.psi)
  if(length(x$pos.weights)==0) x$pos.weights = x$neg.weights*0
  if(length(x$neg.weights)==0) x$neg.weights = x$pos.weights*0
  pright <- c(list(
    stat_identity(aes(x=v, y=w), position = "identity", geom="bar", 
                  data=data.frame(w=x$pos.weights, v=names(x$pos.weights)),
                  fill=gray(poscolwt)),
    scale_y_continuous(name="Positive weights", expand=c(0.000,0.000), breaks=c(0.25, 0.5, 0.75)),
    scale_x_discrete(limits=nms, breaks=nms, labels=nms, drop=FALSE, position="top"),
    geom_hline(aes(yintercept=0)),
    coord_flip(ylim=c(0,1))
  ), theme_butterfly_r)
  pleft <-  c(list(
    stat_identity(aes(x=v, y=w), position = "identity", geom="bar", 
                  data=data.frame(w=x$neg.weights, v=names(x$neg.weights)),
                  fill=gray(1-poscolwt)), 
    scale_y_continuous(trans="reverse", name="Negative weights", expand=c(0.000,0.000), breaks=c(0.25, 0.5, 0.75), limits=c(1,0),labels = waiver()),
    scale_x_discrete(name="Variable", limits=nms, breaks=nms, labels=nms, drop=FALSE),
    geom_hline(aes(yintercept=0)),
    #coord_flip(ylim=c(0,1), expand = FALSE),
    coord_flip()
    ), theme_butterfly_l)
  if((length(x$neg.weights)>0 || length(x$pos.weights)>0)){
    #maxstr = max(mapply(nchar, c(names(x$neg.weights), names(x$pos.weights))))
    maxstr = max(nchar(c(names(x$neg.weights), names(x$pos.weights))))
    lw = 1+maxstr/20
    pleft = ggplot() + pleft
    pright = ggplot() + pright
    p <- gridExtra::arrangeGrob(grobs=list(pleft, pright), ncol=2, padding=0.0, widths=c(lw,1))
  }
  p
}



.plot_noboot_zi <- function(x, theme_butterfly_r, theme_butterfly_l){
  v <- w <- NULL
  # zero inflated
  p1 = list()
  maxcidx=1
  for(modtype in names(x$psi)){
    cidx = grep(paste0("^",modtype), names(unlist(x$msmfit$coefficients)))
    maxcidx = max(cidx, maxcidx)
    nms = unique(names(sort(c(x$pos.weights[[modtype]], x$neg.weights[[modtype]]), decreasing = FALSE)))
    poscolwt = 1-x$pos.psi[[modtype]]/(x$pos.psi[[modtype]] - x$neg.psi[[modtype]])
    if(length(x$pos.weights[[modtype]])==0) x$pos.weights[[modtype]] = x$neg.weights[[modtype]]*0
    if(length(x$neg.weights[[modtype]])==0) x$neg.weights[[modtype]] = x$pos.weights[[modtype]]*0
    pright <- list(
      stat_identity(aes(x=v, y=w), position = "identity", geom="bar", 
                    data=data.frame(w=x$pos.weights[[modtype]], v=names(x$pos.weights[[modtype]])),
                    fill=gray(poscolwt)),
      scale_y_continuous(name="Positive weights", expand=c(0.000,0.000), breaks=c(0.25, 0.5, 0.75)),
      scale_x_discrete(limits=nms, breaks=nms, labels=nms, drop=FALSE, position="top"),
      geom_hline(aes(yintercept=0)),
      coord_flip(ylim=c(0,1)),
      theme_butterfly_r
    ) 
    pleft <- list(
      stat_identity(aes(x=v, y=w), position = "identity", geom="bar", 
                    data=data.frame(w=x$neg.weights[[modtype]], v=names(x$neg.weights[[modtype]])),
                    fill=gray(1-poscolwt)), 
      scale_y_reverse(name="Negative weights", expand=c(0.000,0.000), breaks=c(0.25, 0.5, 0.75)),
      scale_x_discrete(name=paste0("Variable (", modtype, " model)"), limits=nms, breaks=nms, labels=nms, drop=FALSE),
      geom_hline(aes(yintercept=0)),
      coord_flip(ylim=c(0,1)),
      theme_butterfly_l
    )  
    if((length(x$neg.weights[[modtype]])>0 & length(x$pos.weights[[modtype]])>0)){
      #maxstr = max(mapply(nchar, c(names(x$neg.weights[[modtype]]), names(x$pos.weights[[modtype]]))))
      maxstr = max(nchar(c(names(x$neg.weights[[modtype]]), names(x$pos.weights[[modtype]]))))
      lw = 1+maxstr/20
      pleft = ggplot() + pleft
      pright = ggplot() + pright
      p1[[modtype]] <- gridExtra::arrangeGrob(grobs=list(pleft, pright), ncol=2, padding=0.0, widths=c(lw,1))
    }
  }
  p1
}


  
.plot_boot_gaussian <- function(p, x, modelband, flexfit, modelfitline, pointwisebars, pointwiseref=1, alpha=0.05){
  if(!(x$msmfit$family$link == "identity")) stop("Plotting not implemented for this link function")
  p <- c(p, list(
    labs(x = "Joint exposure quantile", y = "Y"),
    lims(x=c(0,1))
  ))
  #
  if(modelband)     p <- c(p, .plot.md.mod.bounds(x,alpha=alpha)) # : add alpha to main function
  if(flexfit)       p <- c(p, .plot.linear.smooth.line(x))
  if(modelfitline)  p <- c(p, .plot.linear.line(x))
  if(pointwisebars) p <- c(p, .plot.md.pw.boot(x, alpha, pointwiseref))
  p
}

.plot_ee_gaussian <- .plot_boot_gaussian


.plot_boot_binomial <- function(p, x, modelband, flexfit, modelfitline, pointwisebars, pointwiseref=1, alpha=0.05){
  if(!(x$msmfit$family$link %in% c("log", "logit"))) stop("Plotting not implemented for this link function")
  #
  p <- c(p, list(scale_y_log10()))
  if(x$msmfit$family$link == "logit"){
    p <- c(p, list(labs(x = "Joint exposure quantile", y = "Odds(Y=1)"),  lims(x=c(0,1))))
    if(modelband)     p <- c(p, .plot.or.mod.bounds(x,alpha))
    if(flexfit)       p <- c(p, .plot.or.smooth.line(x))
    if(modelfitline)  p <- c(p, .plot.logitlin.line(x))
    if(pointwisebars) p <- c(p, .plot.or.pw.boot(x,alpha,pointwiseref))
  } else if(x$msmfit$family$link=="log"){
    p <-c(p, list(labs(x = "Joint exposure quantile", y = "Pr(Y=1)"),  lims(x=c(0,1))))
    if(modelband)     p <- c(p, .plot.rr.mod.bounds(x,alpha))
    if(flexfit)       p <- c(p, .plot.rr.smooth.line(x))
    if(modelfitline)  p <- c(p, .plot.loglin.line(x))
    if(pointwisebars) p <- c(p, .plot.rr.pw.boot(x,alpha,pointwiseref))
  }
  p
}

.plot_ee_binomial <- .plot_boot_binomial


.plot_boot_poisson <- function(p, x, modelband, flexfit, modelfitline, pointwisebars, pointwiseref=1, alpha=0.05){
  if(!(x$msmfit$family$link == "log")) stop("Plotting not implemented for this link function")
  p <- c(p, list(scale_y_log10()))
  if(x$msmfit$family$link == "log"){
    p <- c(p, list(labs(x = "Joint exposure quantile", y = "E(Y)"), lims(x=c(0,1))))
    if(modelband)     p <- c(p, .plot.rr.mod.bounds(x,alpha))
    if(flexfit)       p <- c(p, .plot.rr.smooth.line(x))
    if(modelfitline)  p <- c(p, .plot.loglin.line(x))
    if(pointwisebars) p <- c(p, .plot.rr.pw.boot(x,alpha,pointwiseref))
  }
  p
}

.plot_ee_poisson <- .plot_boot_poisson


.plot_boot_cox <- function(p, x, modelband, flexfit, modelfitline, pointwisebars, pointwiseref=1, alpha=0.05){
  # : make the plot more configurable
  surv <- NULL
  scl = qgcomp.survcurve.boot(x)
  cdf0 = scl$cdfq[scl$cdfq$q==1,]
  cdfmax = scl$cdfq[scl$cdfq$q==x$q,]
  mdf0 = scl$mdfq[scl$mdfq$q==1,]
  mdfmax = scl$mdfq[scl$mdfq$q==x$q,]
  p <- c(p, list(
    geom_step(aes(x=time, y=surv, color="MSM", linetype="Average (all quantiles)"), data=scl$mdfpop),
    geom_step(aes(x=time, y=surv, color="Conditional", linetype="Average (all quantiles)"), data=scl$cdfpop),
    geom_step(aes(x=time, y=surv, color="MSM", linetype="Lowest quantile"), data=mdf0),
    geom_step(aes(x=time, y=surv, color="Conditional", linetype="Lowest quantile"), data=cdf0),
    geom_step(aes(x=time, y=surv, color="MSM", linetype="Highest quantile"), data=mdfmax),
    geom_step(aes(x=time, y=surv, color="Conditional", linetype="Highest quantile"), data=cdfmax),
    scale_y_continuous(name="Survival", limits=c(0,1)),
    scale_x_continuous(name="Time"),
    scale_linetype_discrete(name=""),
    theme(legend.position.inside = c(0.01, 0.01), legend.justification = c(0,0)))
  )
  
  return(p)
}


.plot_boot_zi <- function(p, x, modelband, flexfit, modelfitline, pointwisebars, pointwiseref=1, alpha=0.05){
  # zero inflated
  p <- c(p, list(labs(x = "Joint exposure quantile", y = "E(Y)"),  lims(x=c(0,1))))
  if(modelband)     p <- c(p, .plot.rr.mod.bounds(x,alpha=alpha))
  if(flexfit)       p <- c(p, .plot.linear.smooth.line(x))
  if(modelfitline)  p <- c(p, .plot.linear.line(x))
  if(pointwisebars) p <- c(p, .plot.zi.pw.boot(x, alpha=alpha, pointwiseref))
  p
}

.butterfly_themes <- function(){
  theme_butterfly_l <- list(theme(
    legend.position.inside = c(0,0), 
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
  list(theme_butterfly_l, theme_butterfly_r)
}



######## User facing functions ##########
#' @title Default plotting method for a qgcompfit object
#'
#' @description Plot a quantile g-computation object. For qgcomp.glm.noboot, this function will
#' create a butterfly plot of weights. For qgcomp.glm.boot, this function will create
#' a box plot with smoothed line overlaying that represents a non-parametric
#' fit of a model to the expected outcomes in the population at each quantile
#' of the joint exposures (e.g. '1' represents 'at the first quantile for
#' every exposure')
#' 
#' @param x "qgcompfit" object from `qgcomp.glm.noboot`,  `qgcomp.glm.boot`, 
#'   `qgcomp.cox.noboot`,  `qgcomp.cox.boot`, `qgcomp.zi.noboot` or `qgcomp.zi.boot` functions
#' @param suppressprint If TRUE, suppresses the plot, rather than printing it 
#'   by default (it can be saved as a ggplot2 object (or list of ggplot2 objects if x is from a zero-
#'   inflated model) and used programmatically)
#'   (default = FALSE)
#' @param geom_only If TRUE, returns only the geometry (i.e. does not contain the entire plot object). Used for overlays. Only used for `.ee` and `.boot` methods.
#'   (default = FALSE)
#' @param pointwisebars (boot.gcomp only) If TRUE, adds 95%  error bars for pointwise comparisons
#' of E(Y|joint exposure) to the smooth regression line plot
#' @param modelfitline (boot.gcomp only) If TRUE, adds fitted (MSM) regression line
#' of E(Y|joint exposure) to the smooth regression line plot
#' @param modelband If TRUE, adds 95% prediction bands for E(Y|joint exposure) (the MSM fit)
#' @param flexfit (boot.gcomp only) if TRUE, adds flexible interpolation of predictions from 
#' underlying (conditional) model
#' @param pointwiseref (boot.gcomp only) integer: which category of exposure (from 1 to q) 
#' should serve as the referent category for pointwise comparisons? (default=1)
#' @param ... unused
#' @seealso \code{\link[qgcomp]{qgcomp.glm.noboot}}, \code{\link[qgcomp]{qgcomp.glm.boot}}, and \code{\link[qgcomp]{qgcomp}}
#' @import ggplot2 grid gridExtra
#' @importFrom grDevices gray
#' @export
#' @examples
#' set.seed(12)
#' dat <- data.frame(x1=(x1 <- runif(100)), x2=runif(100), x3=runif(100), z=runif(100),
#'                   y=runif(100)+x1+x1^2)
#' ft <- qgcomp.glm.noboot(y ~ z + x1 + x2 + x3, expnms=c('x1','x2','x3'), data=dat, q=4)
#' ft
#' # display weights
#' plot(ft)
#' # examining fit
#' plot(ft$fit, which=1) # residual vs. fitted is not straight line!
#' \dontrun{
#' 
#' # using non-linear outcome model
#' ft2 <- qgcomp.glm.boot(y ~ z + x1 + x2 + x3 + I(x1*x1), expnms=c('x1','x2','x3'), 
#' data=dat, q=4, B=10)
#' ft2
#' plot(ft2$fit, which=1) # much better looking fit diagnostics suggests
#' # it is better to include interaction term for x
#' plot(ft2) # the msm predictions don't match up with a smooth estimate
#' # of the expected outcome, so we should consider a non-linear MSM
#'
#' # using non-linear marginal structural model
#' ft3 <- qgcomp.glm.boot(y ~ z + x1 + x2 + x3 + I(x1*x1), expnms=c('x1','x2','x3'), 
#' data=dat, q=4, B=10, degree=2)
#' # plot(ft3$fit, which=1) - not run - this is identical to ft2 fit
#' plot(ft3) # the MSM estimates look much closer to the smoothed estimates
#' # suggesting the non-linear MSM fits the data better and should be used
#' # for inference about the effect of the exposure
#' 
#' # binary outcomes, logistic model with or without a log-binomial marginal 
#' structural model
#' dat <- data.frame(y=rbinom(100,1,0.5), x1=runif(100), x2=runif(100), z=runif(100))
#' fit1 <- qgcomp.glm.boot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), 
#'          data=dat, q=9, B=100, rr=FALSE)
#' fit2 <- qgcomp.glm.boot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), 
#'          data=dat, q=9, B=100, rr=TRUE)
#' plot(fit1)
#' plot(fit2)
#' # Using survival data ()
#' set.seed(50)
#' N=200
#' dat <- data.frame(time=(tmg <- pmin(.1,rweibull(N, 10, 0.1))), 
#'                 d=1.0*(tmg<0.1), x1=runif(N), x2=runif(N), z=runif(N))
#' expnms=paste0("x", 1:2)
#' f = survival::Surv(time, d)~x1 + x2
#' (fit1 <- survival::coxph(f, data = dat))
#' # non-bootstrap method to get a plot of weights
#' (obj <- qgcomp.cox.noboot(f, expnms = expnms, data = dat))
#' plot(obj)
#' 
#' # bootstrap method to get a survival curve
#' # this plots the expected survival curve for the underlying (conditional) model
#' # as well as the expected survival curve for the MSM under the following scenarios:
#' #  1) highest joint exposure category
#' #  2) lowest joint exposure category
#' #  3) average across all exposure categories 
#' # differences between the MSM and conditional fit suggest that the MSM is not flexible
#' # enough to accomodate non-linearities in the underlying fit (or they may simply signal that
#' # MCSize should be higher). Note that if linearity
#' # is assumed in the conditional model, the MSM will typically also appear linear and
#' # will certainly appear linear if no non-exposure covariates are included in the model
#' # not run (slow when using boot version to proper precision)
#' (obj2 <- qgcomp.cox.boot(f, expnms = expnms, data = dat, B=10, MCsize=2000))
#' plot(obj2)
#' }
plot.qgcompfit <- function(x, 
                           suppressprint=FALSE, 
                           geom_only=FALSE,
                           pointwisebars=TRUE, 
                           modelfitline=TRUE, 
                           modelband=TRUE, 
                           flexfit=TRUE, 
                           pointwiseref = ceiling(x$q/2),
                           ...){
  
  requireNamespace("ggplot2")
  requireNamespace("grid")
  requireNamespace("gridExtra")
  ymin <- ymax <- w <- v <- NULL # appease R CMD check
  isesteq = inherits(x, "eeqgcompfit")
  iszi = is.null(x$fit$family)
  isboot = x$bootstrap
  if(geom_only & !suppressprint){
    message("geomonly=TRUE, suppressing output by setting suppressprint=TRUE")
    suppressprint=TRUE
  }  #vpl <- grid::viewport(width=0.525, height=1, x=0, y=0, just=c("left", "bottom"))
  #vpr <- grid::viewport(width=0.475, height=1, x=0.525, y=0, just=c("left", "bottom"))
  if(!isboot & !isesteq){
    themes = .butterfly_themes()
    theme_butterfly_l = themes[[1]]
    theme_butterfly_r = themes[[2]]
    # zero inflated
    if(iszi){
      p <- .plot_noboot_zi(x, theme_butterfly_r, theme_butterfly_l) # arrangeGrob object
      if(!suppressprint) {
        lapply(p, .plfun)
      }
    } else {
      nms = unique(names(sort(c(x$pos.weights, x$neg.weights), decreasing = FALSE)))
      p <- .plot_noboot_base(x, nms, theme_butterfly_r, theme_butterfly_l) # arrangeGrob object
      if(!suppressprint) {
        .plfun(p)
      }
    }
  }
  if(isesteq){
    #p <- ggplot() 
    p <- list()
    if(x$msmfit$family$family=='gaussian') temppfun <- .plot_ee_gaussian
    if(x$msmfit$family$family=='binomial') temppfun <- .plot_ee_binomial
    if(x$msmfit$family$family=='poisson')  temppfun <- .plot_ee_poisson
    p <- temppfun(p, x, modelband, flexfit, modelfitline, pointwisebars, pointwiseref, alpha=0.05)
    p <- c(p, list(
      scale_fill_grey(name="", start=.9), 
      scale_colour_grey(name="", start=0.0, end=0.6),
      theme_classic()
    ))
    if (!suppressprint){
      print(ggplot2::ggplot() + p)
    }
  }
  
  if(isboot & !isesteq){
    # variance based on delta method and knowledge that non-linear
    #functions will always be polynomials in qgcomp
    # default plot for bootstrap results (no weights obtained)
    p = list()
    if(iszi) {
      temppfun <- .plot_boot_zi
    } else{
      # Cox model or standard GLM
      if(x$msmfit$family$family=='cox')      temppfun <- .plot_boot_cox
      if(x$msmfit$family$family=='gaussian') temppfun <- .plot_boot_gaussian
      if(x$msmfit$family$family=='binomial') temppfun <- .plot_boot_binomial
      if(x$msmfit$family$family=='poisson')  temppfun <- .plot_boot_poisson
    }
    p <- temppfun(p, x, modelband, flexfit, modelfitline, pointwisebars, pointwiseref, alpha=0.05)
    p <- c(p, list(
      scale_fill_grey(name="", start=.9),
      scale_colour_grey(name="", start=0.0, end=0.6),
      theme_classic()
    ))
    if (!suppressprint){
      print(ggplot2::ggplot() + p)
    }
  }
  if (suppressprint & geom_only ) return(p)
  if (suppressprint & !geom_only ) return(ggplot2::ggplot() + p)
}


#' @export
#' @describeIn plot.qgcompfit Plot method for qgcomp multinomial fits
plot.qgcompmultfit <- function(
    x,                           
    suppressprint=FALSE, 
    pointwisebars=TRUE, 
    modelfitline=TRUE, 
    modelband=TRUE, 
    flexfit=TRUE, 
    pointwiseref = ceiling(x$q/2),
    ...
){
  if(x$bootstrap){
    warning("The default plot function is only functional for qgcomp.multinomial.noboot currently")
  } else{
    nms = list()
    for (r in seq_len(nrow(x$weights))){
      nms[[r]] = names(sort(-abs(x$weights[r,])))
    }
    themes = .butterfly_themes()
    theme_butterfly_l = themes[[1]]
    theme_butterfly_r = themes[[2]]
    plist <- list()
    for(r in seq_len(nrow(x$weights))){
      plist[[r]] <- .plot_noboot_multi_base(r, x, nms, theme_butterfly_r, theme_butterfly_l)
    }
    
    p <- gridExtra::arrangeGrob(grobs=plist)
    
    if(!suppressprint) {
      .plfun(p) # arrangeGrob object
    }
    
  }
}
