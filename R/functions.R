
se_comb <- function(weightvec, covmat){
  #' qgcomp::se_comb
  #'
  #' calculate standard error of weighted linear combination of random variables
  #'  given a vector of weights and a covariance matrix
  #' @param weightvec a logical index of which columns are used to create
  #' composite effect estimate
  #' @param covmat covariance matrix from a glm fit
  #' @keywords variance, mixtures
  #' @examples
  #' dat = data.frame(y=runif(10), x=runif(10))
  #' lmfit = glm(y ~ x, data=dat, family='gaussian')
  #' qgcomp:::se_comb(weightvec=c(0,1), covmat=summary(lmfit)$cov.scaled)

  #calculate standard error of weighted linear combination of random variables
  wcovmat <- weightvec %*% t(weightvec) * covmat
  var <- sum(wcovmat)
  sqrt(var)
}

quantize <- function (data, expcoefs, expnms, q) {
  #' qgcomp::quantize
  #'
  #' create variables representing indicator functions with cutpoints defined
  #' by quantiles
  #' @param data a data frame
  #' @param expcoefs a logical index of which columns contain the columns to be
  #' quantized
  #' @param expnms a character vector with the names of  the columns to be
  #' quantized
  #' @param q integer, number of quantiles used in creating quantized variables
  #' @keywords variance, mixtures
  #' @import stats
  #' @export
  #' @examples
  #' dat = data.frame(y=runif(10), x1=runif(10), x2=runif(10), z=runif(10))
  #' qdata = quantize(data=dat, expcoefs=c(0,1,1,0), expnms=c("x1", "x2"), q=2)
    Xnms <- names(data[, expnms[which(as.logical(expcoefs[-1]))]])
    for (i in 1:length(Xnms)) {
        dat_num <- as.numeric(unlist(data[, Xnms[i]]))
        data[[Xnms[i]]] <- cut(dat_num, breaks = unique(quantile(dat_num,
             probs = seq(0, 1, by = 1 / q), na.rm = TRUE)), labels = FALSE,
             include.lowest = TRUE) - 1
    }
    return(data)
}


qgcomp.noboot <- function(f, data, expcoefs=NULL, q=4, alpha=0.05, ...){
  #' qgcomp::qgcomp.noboot
  #'
  #' create variables representing indicator functions with cutpoints defined
  #' by quantiles
  #' @param f R style formula
  #' @param data data frame
  #' @param expcoefs logical/numeric vector that points to variables in the
  #' formula that should be considered part of the 'mixture' (not counting intercept!)
  #' @param q number of quantiles used to create quantile indicator variables
  #' representing the exposure variables
  #' @param alpha alpha level for confidence limit calculation
  #' @param ... arguments to glm (e.g. family)
  #' @keywords variance, mixtures
  #' @import stats
  #' @export
  #' @examples
  #' dat = data.frame(y=runif(10), x1=runif(10), x2=runif(10), z=runif(10))
  #' qgcomp.noboot(y ~ z + x1 + x2, expcoefs=c(0,1,1), data=dat, q=2)
   expnms <- attr(terms(f, data = data), "term.labels")
    if (is.null(expcoefs)) {
      expcoefs <- c(0, rep(1, length(expnms)))
    } else expcoefs <- c(0, expcoefs) # append term for intercept
    if (!is.null(q)){
      qdata <- quantize(data, expcoefs, expnms, q)
    } else qdata <- data
    fit <- glm(f, data = qdata, ...)
    mod <- summary(fit)
    estb <- sum(mod$coefficients[which(as.logical(expcoefs)), 1])
    seb <- se_comb(expcoefs, covmat = mod$cov.scaled)
    tstat <- estb / seb
    df <- mod$df.null - sum(expcoefs)
    pval <- 2 - 2 * pt(abs(tstat), df = df)
    pvalz <- 2 - 2 * pnorm(abs(tstat))
    ci <- c(estb + seb * qnorm(alpha / 2), estb + seb * qnorm(1 - alpha / 2))
    # 'weights'
    wcoef <- fit$coefficients[which(as.logical(expcoefs))]
    names(wcoef) <- gsub("_q", "", names(wcoef))
    poscoef <- which(wcoef > 0)
    pweights <- abs(wcoef[poscoef]) / sum(abs(wcoef[poscoef]))
    nweights <- abs(wcoef[-poscoef]) / sum(abs(wcoef[-poscoef]))
    # 'post-hoc' positive and negative estimators 
    # similar to constrained gWQS
    pos.gamma <- sum(wcoef[poscoef])
    neg.gamma <- sum(wcoef[-poscoef])
    se.pos.gamma <- se_comb(expcoefs*(fit$coefficients>0), covmat = mod$cov.scaled)
    se.neg.gamma <- se_comb(expcoefs*(fit$coefficients<0), covmat = mod$cov.scaled)
    qx <- qdata[, expnms]
    names(qx) <- paste0(names(qx), "_q")
    res <- list(
      qx = qx, fit = fit, gamma = estb, var.gamma = seb ^ 2, ci = ci,
      pos.gamma = pos.gamma, var.pos.gamma = se.pos.gamma^2,
      neg.gamma = neg.gamma, var.neg.gamma = se.neg.gamma^2,
                pweights = sort(pweights, decreasing = TRUE),
                nweights = sort(nweights, decreasing = TRUE), 
                psize = sum(abs(wcoef[poscoef])),
                nsize = sum(abs(wcoef[-poscoef]))
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

print.qgcompfit <- function(x, ...){
  #' qgcomp::print.qgcompfit
  #'
  #' carry out quantile g-computation
  #' by quantiles
  #' @param x "qgcompfit" object from `qgcomp.noboot` function
  #' @param ... unused
  #' @keywords variance, mixtures
  #' @export
  #' @examples
  #' dat = data.frame(y=runif(10), x1=runif(10), x2=runif(10), z=runif(10))
  #' qgcomp.noboot(y ~ z + x1 + x2, expcoefs=c(0,1,1), data=dat, q=2)
  fam <- x$fit$family$family
  cat(paste0("Scaled effect size (positive direction, sum of positive coefficients = ", signif(x$psize, 3) , ")\n"))
  if (length(x$pweights) > 0) {
    print(x$pweights, digits = 3)
  } else cat("None\n")
  cat("\n")
  cat(paste0("Scaled effect size (negative direction, sum of negative coefficients = ", signif(-x$nsize, 3) , ")\n"))
  if (length(x$nweights) > 0) {
    print(x$nweights, digits = 3)
  } else cat("None\n")
  cat("\n")
  if (fam == "binomial"){
    cat("Mixture log(OR):\n")
    cat(paste0("gamma (CI): ", signif(x$gamma, 3), " (",
             signif(x$ci[1], 3), ",", signif(x$ci[2], 3), "), z=",
             signif(x$zstat, 3), ", p=",
             signif(x$pval, 3), "\n"))
  }
  if (fam == "gaussian"){
    cat("Mixture slope:\n")
    cat(paste0("gamma (CI): ", signif(x$gamma, 3), " (",
             signif(x$ci[1], 3), ",", signif(x$ci[2], 3), "), t=",
             signif(x$tstat, 3), ", df=", x$df, ", p=",
             signif(x$pval, 3), "\n"))
  }
}

plot.qgcompfit <- function(x, ...){
  #' qgcomp::plot.qgcompfit
  #'
  #' plot quantile g-computation object
  #' @param x "qgcompfit" object from `qgcomp.noboot` function
  #' @param ... unused
  #' @keywords variance, mixtures
  #' @import ggplot2 grid gridExtra
  #' @export
  #' @examples
  #' dat = data.frame(y=runif(10), x1=runif(10), x2=runif(10), z=runif(10))
  #' qgcomp.noboot(y ~ z + x1 + x2, expcoefs=c(0,1,1), data=dat, q=2)

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
  
    maxstr = max(mapply(nchar, c(names(x$nweights), names(x$pweights))))
    lw = 1+maxstr/20
    p1 <- gridExtra::arrangeGrob(grobs=list(pleft, pright), ncol=2, padding=0.0, widths=c(lw,1))
    grid::grid.newpage()
    grid::grid.draw(p1)
  #grid.text("Density", x=0.55, y=0.1, gp=gpar(fontsize=14, fontface="bold", fontfamily="Helvetica"))
}
