## ---- echo=FALSE, results='asis', message=FALSE--------------------------
library(gWQS)
library(qgcomp)
library(knitr)
knitr::kable(head(wqs_data[, c(37, 36, 35, 1:34)], 10))

## ---- results='asis', fig.show='hold', fig.height=5, fig.width=5, cache=TRUE----

# we save the names of the mixture variables in the variable "Xnm"
data("wqs_data", package="gWQS")

Xnm = c("log_LBX074LA", "log_LBX099LA", "log_LBX105LA", "log_LBX118LA",
"log_LBX138LA", "log_LBX153LA", "log_LBX156LA", "log_LBX157LA", "log_LBX167LA",
"log_LBX170LA", "log_LBX180LA", "log_LBX187LA", "log_LBX189LA", "log_LBX194LA",
"log_LBX196LA", "log_LBX199LA", "log_LBXD01LA", "log_LBXD02LA", "log_LBXD03LA",
"log_LBXD04LA", "log_LBXD05LA", "log_LBXD07LA", "log_LBXF01LA", "log_LBXF02LA",
"log_LBXF03LA", "log_LBXF04LA", "log_LBXF05LA", "log_LBXF06LA", "log_LBXF07LA",
"log_LBXF08LA", "log_LBXF09LA", "log_LBXPCBLA", "log_LBXTCDLA", "log_LBXHXCLA")


# we run the model and save the results in the variable "results"
results = qgcomp.noboot(y~.,dat=wqs_data[,c(Xnm, 'y')], family=gaussian())


# we compare...

gcompmod = qgcomp.noboot(disease_state~.,wqs_data[,c(Xnm, 'disease_state')], family=binomial(), q=4)

wqsmod = gwqs(disease_state ~ 1, mix_name = Xnm, data = wqs_data, q = 4, 
     validation = 0.6, b = 3, b1_pos = F, b1_constr = F, family='binomial', seed=125)
     

wqsmod$final_weights
summary(wqsmod$fit)
gcompmod



## ---- echo=FALSE, results='asis', message=FALSE--------------------------
knitr::kable(results$pweights, digits = 3, row.names = F)
knitr::kable(results$nweights, digits = 3, row.names = F)

## ---- results='asis', message=FALSE, eval=F------------------------------
#  print(results)

## ---- results='asis', fig.show='hold', fig.height=5, fig.width=5, cache=TRUE----

gcompmod = qgcomp.noboot(disease_state~sex+log_LBX074LA+log_LBX099LA+log_LBX105LA+log_LBX118LA+log_LBX138LA+log_LBX153LA+
                           log_LBX156LA+log_LBX157LA+log_LBX167LA+log_LBX170LA+log_LBX180LA+log_LBX187LA+log_LBX189LA+
                           log_LBX194LA+log_LBX196LA+log_LBX199LA+log_LBXD01LA+log_LBXD02LA+log_LBXD03LA+log_LBXD04LA+
                           log_LBXD05LA+log_LBXD07LA+log_LBXF01LA+log_LBXF02LA+log_LBXF03LA+log_LBXF04LA+log_LBXF05LA+
                           log_LBXF06LA+log_LBXF07LA+log_LBXF08LA+log_LBXF09LA+log_LBXPCBLA+log_LBXTCDLA+log_LBXHXCLA,
                         expnms=Xnm,
                         wqs_data[,c(Xnm, "sex", 'disease_state')], family=binomial(), q=4)
wqsmod = gwqs(disease_state ~ sex, mix_name = Xnm, data = wqs_data, q = 4, 
     validation = 0.6, b = 3, b1_pos = F, b1_constr = F, family='binomial', seed=125)
     


## ---- results='asis', message=FALSE, eval=F------------------------------
#  print(gcompmod)
#  summary(gcompmod$fit)

## ---- echo=F, results='asis', message=FALSE------------------------------
wqsmod$final_weights
summary(wqsmod$fit)


