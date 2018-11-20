## ---- echo=FALSE, results='markup', message=FALSE------------------------
library(gWQS)
library(qgcomp)
library(knitr)
head(wqs_data[, c(37, 36, 35, 1:34)], 10)

## ---- results='markup', fig.show='hold', fig.height=5, fig.width=5, cache=TRUE----

# we save the names of the mixture variables in the variable "Xnm"
data("wqs_data", package="gWQS")

Xnm <- c("log_LBX074LA", "log_LBX099LA", "log_LBX105LA", "log_LBX118LA",
"log_LBX138LA", "log_LBX153LA", "log_LBX156LA", "log_LBX157LA", "log_LBX167LA",
"log_LBX170LA", "log_LBX180LA", "log_LBX187LA", "log_LBX189LA", "log_LBX194LA",
"log_LBX196LA", "log_LBX199LA", "log_LBXD01LA", "log_LBXD02LA", "log_LBXD03LA",
"log_LBXD04LA", "log_LBXD05LA", "log_LBXD07LA", "log_LBXF01LA", "log_LBXF02LA",
"log_LBXF03LA", "log_LBXF04LA", "log_LBXF05LA", "log_LBXF06LA", "log_LBXF07LA",
"log_LBXF08LA", "log_LBXF09LA", "log_LBXPCBLA", "log_LBXTCDLA", "log_LBXHXCLA")


# we run the model and save the results in the variable "results"
results <- qgcomp.noboot(y~.,dat=wqs_data[,c(Xnm, 'y')], family=gaussian())


# we compare a qgcomp.noboot fit:
gcompmod <- qgcomp.noboot(disease_state~., expnms=Xnm, data = wqs_data[,c(Xnm, 'disease_state')], family=binomial(), q=4)

# and a qgcomp.boot fit:
gcompmod2 <- qgcomp.boot(disease_state~., expnms=Xnm, data = wqs_data[,c(Xnm, 'disease_state')], family=binomial(), q=4, B=200)


# with a gwqs fit:
suppressWarnings(wqsmod <- gwqs(disease_state ~ 1, mix_name = Xnm, data = wqs_data, q = 4, validation = 0.6, b = 3, b1_pos = F, b1_constr = F, family='binomial', seed=125))
     

# WQS fit (with reduced number of bootstraps to save time)
wqsmod$final_weights
wqsmod$fit
# conditional OR fit, similar to WQS
gcompmod
# population average RR with bootstrap confidence intervals
gcompmod2



## ---- results='markup', fig.show='hold', fig.height=5, fig.width=5, cache=TRUE----

gcompmod <- qgcomp.noboot(y~sex+log_LBX074LA+log_LBX099LA+log_LBX105LA+log_LBX118LA+log_LBX138LA+log_LBX153LA+
                           log_LBX156LA+log_LBX157LA+log_LBX167LA+log_LBX170LA+log_LBX180LA+log_LBX187LA+log_LBX189LA+
                           log_LBX194LA+log_LBX196LA+log_LBX199LA+log_LBXD01LA+log_LBXD02LA+log_LBXD03LA+log_LBXD04LA+
                           log_LBXD05LA+log_LBXD07LA+log_LBXF01LA+log_LBXF02LA+log_LBXF03LA+log_LBXF04LA+log_LBXF05LA+
                           log_LBXF06LA+log_LBXF07LA+log_LBXF08LA+log_LBXF09LA+log_LBXPCBLA+log_LBXTCDLA+log_LBXHXCLA,
                         expnms=Xnm,
                         wqs_data, family=gaussian(), q=4)
gcompmod2 <- qgcomp.boot(y~sex+log_LBX074LA+log_LBX099LA+log_LBX105LA+log_LBX118LA+log_LBX138LA+log_LBX153LA+
                           log_LBX156LA+log_LBX157LA+log_LBX167LA+log_LBX170LA+log_LBX180LA+log_LBX187LA+log_LBX189LA+
                           log_LBX194LA+log_LBX196LA+log_LBX199LA+log_LBXD01LA+log_LBXD02LA+log_LBXD03LA+log_LBXD04LA+
                           log_LBXD05LA+log_LBXD07LA+log_LBXF01LA+log_LBXF02LA+log_LBXF03LA+log_LBXF04LA+log_LBXF05LA+
                           log_LBXF06LA+log_LBXF07LA+log_LBXF08LA+log_LBXF09LA+log_LBXPCBLA+log_LBXTCDLA+log_LBXHXCLA,
                         expnms=Xnm,
                         wqs_data, family=gaussian(), q=4, B=200)
wqsmod <- gwqs(y ~ sex, mix_name = Xnm, data = wqs_data, q = 4, 
     validation = 0.6, b = 3, b1_pos = TRUE, b1_constr = F, family='gaussian', seed=125, plots=FALSE)
     


## ---- results='markup', message=FALSE, fig.height=7, fig.width=6---------
plot(gcompmod)
plot(gcompmod2)

## ---- echo=F, results='markup', message=FALSE, fig.height=7, fig.width=6----
wqsmod <- gwqs(y ~ sex, mix_name = Xnm, data = wqs_data, q = 4, 
     validation = 0.6, b = 3, b1_pos = TRUE, b1_constr = F, family='gaussian', seed=125, plots=TRUE)


