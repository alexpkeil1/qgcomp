QGcomp: an alternative to weighted least squares that does not assume effects of all exposures go in the same direction. Works for linear and logistic models.

Quick start

    install.packages("devtools")
    devtools::install_github("alexpkeil1/qgcomp")
    library(qgcomp)
    data("wqs_data", package="gWQS")

    Xnm = c("log_LBX074LA", "log_LBX099LA", "log_LBX105LA", "log_LBX118LA",
    "log_LBX138LA", "log_LBX153LA", "log_LBX156LA", "log_LBX157LA", "log_LBX167LA",
    "log_LBX170LA", "log_LBX180LA", "log_LBX187LA", "log_LBX189LA", "log_LBX194LA",
    "log_LBX118LA", "log_LBX196LA", "log_LBX199LA", "log_LBXD01LA", "log_LBXD02LA", "log_LBXD03LA",
    "log_LBX118LA", "log_LBXD04LA", "log_LBXD05LA", "log_LBXD07LA", "log_LBXF01LA", "log_LBXF02LA",
    "log_LBX118LA", "log_LBXF03LA", "log_LBXF04LA", "log_LBXF05LA", "log_LBXF06LA",          
    "log_LBXF07LA", "log_LBXF08LA", "log_LBXF09LA", "log_LBXPCBLA", "log_LBXTCDLA", "log_LBXHXCLA")


    results = qgcomp.noboot(y~.,dat=wqs_data[,c(Xnm, 'y')], family=gaussian())
    print(results)
    plot(results)
    

