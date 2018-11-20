---
title: "How to use qgcomp package"
author: "Alexander Keil"
date: "2018-11-19"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{How to use qgcomp package}
  %\VignetteEngine{knitr::knitr}
  \usepackage[utf8]{inputenc}
---

## Introduction

## How to use the `qgcomp` package

### Example 1



|          y| disease_state| sex| log_LBX074LA| log_LBX099LA| log_LBX105LA| log_LBX118LA| log_LBX138LA| log_LBX153LA| log_LBX156LA| log_LBX157LA| log_LBX167LA| log_LBX170LA| log_LBX180LA| log_LBX187LA| log_LBX189LA| log_LBX194LA| log_LBX196LA| log_LBX199LA| log_LBXD01LA| log_LBXD02LA| log_LBXD03LA| log_LBXD04LA| log_LBXD05LA| log_LBXD07LA| log_LBXF01LA| log_LBXF02LA| log_LBXF03LA| log_LBXF04LA| log_LBXF05LA| log_LBXF06LA| log_LBXF07LA| log_LBXF08LA| log_LBXF09LA| log_LBXPCBLA| log_LBXTCDLA| log_LBXHXCLA|
|----------:|-------------:|---:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|------------:|
| -0.4224332|             1|   0|    0.0680716|   -0.1608614|    0.9648013|   -0.1815684|    0.7506797|    0.9363278|   -0.0364088|    0.1971499|    1.1192873|    0.2972489|    0.1920384|   -0.5257044|   -0.1659253|    0.4017490|   -0.0143133|    0.8692975|    0.0846309|    0.0012511|   -0.2847784|    0.6327454|    0.4078422|    0.1661874|   -1.2905580|   -1.0359049|    1.6893966|   -0.1942655|    0.4322813|   -1.8049708|   -0.7816081|    0.2459645|   -1.0917383|   -0.3476487|   -0.0358878|    1.2711118|
|  0.7613784|             0|   0|    0.5276676|    1.4403102|    0.7753093|    0.1055764|   -0.9015101|    0.1677693|    0.7643587|   -0.6001477|    0.6175614|   -0.3577213|    0.2879546|   -0.1557095|   -0.6755189|    0.0519114|    0.3351334|   -1.1385368|    0.5852100|   -1.4539998|    1.1116691|    0.0767604|   -0.0295292|   -1.2578515|   -0.0246641|    0.5600564|    0.1143720|    1.1365974|    0.8777775|   -0.1214907|    0.7902525|    0.2255480|   -1.0510728|    0.5891812|   -0.6304355|    0.7265865|
| -2.1113324|             1|   1|   -1.8650139|    1.4885628|   -1.5420050|   -0.5412270|   -1.3624492|   -0.4177779|    0.1414001|    0.1333590|   -0.9901708|    0.1024100|   -0.0622932|    1.0824292|   -0.6920183|   -1.7204760|   -1.2079225|    0.3543356|   -0.8116185|   -0.8790250|    0.4616429|   -2.0101411|   -1.2859585|    0.8281122|    0.0921212|   -0.6980090|    0.2987547|   -0.0450820|    0.4917081|    0.3754392|   -1.6193861|    0.0812163|    0.2236327|   -0.5211774|   -1.0174183|   -0.4092153|
|  2.2222394|             0|   0|    1.2585209|    1.6110700|    2.2069536|    1.3286183|    1.8871201|    1.8393058|    2.2751535|    1.1307767|    1.2204913|    0.1958512|    0.9461497|   -0.0325942|   -0.3940764|    0.2984753|    0.7837346|    1.9352485|    0.1973576|    0.8376586|   -1.2715081|   -1.2577145|   -0.5786619|    0.6304414|    1.3148279|    0.2485757|   -0.4723713|    0.7557946|    0.3307601|   -0.0342781|   -0.2875164|    0.4176680|    0.0977815|   -1.1273511|    1.0095760|    0.6639172|
| -0.2190805|             1|   0|   -1.0138713|   -1.9903436|   -0.6009229|   -0.8907467|   -0.4789443|   -2.8190387|   -0.0095402|   -1.5477137|    0.5180853|   -2.1491031|   -2.7951863|   -1.5755779|   -0.5270997|   -0.2450904|   -1.1242993|   -0.6902691|   -0.5791061|   -1.9495504|   -1.3426006|    0.6178038|   -1.2587666|   -0.3310535|   -1.5882632|   -1.3967154|   -1.2128921|   -0.8581706|   -0.6420270|   -1.6946821|   -0.2143440|   -1.4720725|   -0.5997930|   -1.0896972|   -0.3474791|   -0.7175427|
| -1.7282770|             0|   0|   -0.9752573|   -0.6864727|    0.3389021|    0.5984246|    0.0270300|   -0.4798327|   -0.3867123|   -1.0377681|   -0.7329465|   -0.7148301|    0.4060471|   -0.5076881|    0.2215813|    0.3212146|   -0.1809994|    1.4622912|    0.0165252|   -0.6026058|   -0.0865081|   -0.7574191|   -0.0918236|   -0.3804055|   -0.0150627|    0.5944020|    0.4030363|   -0.1233755|   -0.7816779|   -0.3795485|    0.4316160|    1.7915524|   -0.9154308|    0.2882417|    0.6230951|    1.4571064|
|  0.4219639|             0|   1|   -0.6225983|    1.0113598|   -0.2162429|    0.7575008|   -0.2571729|    0.9804511|    0.4458043|    1.3962681|   -0.1979610|    0.1995122|    1.4669089|    1.0634414|    0.8339600|    1.5867141|    0.7249490|    0.7584681|    0.2701363|    0.9317193|    0.0963273|    0.3944993|    0.6827900|   -0.5437556|    0.9504820|    0.5142011|    0.5474847|    0.4452620|   -0.4822166|    0.7318831|    0.0323481|   -0.9661976|    0.7297474|   -0.1620114|    0.5754746|    0.9274046|
| -0.4762419|             0|   1|   -0.1801051|   -0.5805092|   -0.8038376|   -0.7786034|   -1.5457539|   -1.3117169|   -0.9734300|   -0.0261900|    0.5308421|   -0.2485865|    0.2103118|   -0.0912344|   -0.1011264|   -0.4217051|    0.4682131|   -0.9008740|   -0.4390236|   -1.0384994|   -1.0510064|   -0.1617107|   -1.5014766|   -1.1384429|    1.7696655|    1.9844850|   -0.4584083|   -1.8996102|   -1.0736679|    1.2719959|    1.2918048|   -0.0684242|    0.4754351|   -1.0453685|    0.0445513|   -1.8944198|
| -0.5126990|             1|   0|    0.4502187|    0.7892304|    0.7107290|    0.2754549|    0.0613806|    0.3399312|    1.0176518|    1.1807174|    1.2440091|    1.5589123|   -0.8405278|   -0.2002233|    0.2201809|    0.2284441|   -0.1742271|   -0.1887278|    0.3027553|   -0.3797491|   -0.4089264|   -0.5428011|    1.2964970|   -0.2374201|    0.1206623|    0.4179228|    1.3594435|    1.6920543|   -0.3431321|   -0.2684203|    0.9664860|   -0.0108648|   -0.9556556|    2.0269212|    0.0600171|    0.4541168|
| -0.8479847|             1|   0|   -0.4461221|   -1.0556249|    1.0676569|    0.0681564|   -0.9582167|   -0.3967903|    0.0196100|    0.0467119|    0.5172333|    0.2781621|    0.1674821|    1.1879319|   -1.0221745|   -1.3879976|   -0.7295127|   -1.7526492|   -0.7668207|    0.5923798|    1.0463819|    0.9856947|   -1.1673793|   -0.7855519|   -0.8638427|   -0.7379376|   -0.0984674|   -0.7446359|   -0.2253229|   -0.9085362|   -0.8710364|    1.6918778|   -0.8597897|   -0.4343029|    0.4885121|   -1.5592026|

qgcomp with a continuous outcome:
This script calls a wqs model for a continuous outcome using the function `gwqs`. 


```r
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
```

Including all model terms as exposures of interest

```r
# we compare...

gcompmod = qgcomp.noboot(disease_state~.,wqs_data[,c(Xnm, 'disease_state')], family=binomial(), q=4)
```

Including all model terms as exposures of interest

```r
wqsmod = gwqs(disease_state ~ 1, mix_name = Xnm, data = wqs_data, q = 4, 
     validation = 0.6, b = 3, b1_pos = F, b1_constr = F, family='binomial', seed=125)
```

```
## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function

## Warning in .safefunx(tmpv, .solnp_fun, .env, ...): 
## solnp-->warning: NaN detected in function call...check your function
```

[1] "The optimization function did not converge 0 times"

```r
wqsmod$final_weights
```

                 mix_name  mean_weight
log_LBXF01LA log_LBXF01LA 1.594500e-01
log_LBX167LA log_LBX167LA 1.512775e-01
log_LBX180LA log_LBX180LA 1.389033e-01
log_LBX170LA log_LBX170LA 1.293262e-01
log_LBXF04LA log_LBXF04LA 1.112770e-01
log_LBX196LA log_LBX196LA 5.989367e-02
log_LBXD02LA log_LBXD02LA 4.308660e-02
log_LBXPCBLA log_LBXPCBLA 4.238708e-02
log_LBX118LA log_LBX118LA 2.748439e-02
log_LBX189LA log_LBX189LA 2.641851e-02
log_LBXHXCLA log_LBXHXCLA 2.616894e-02
log_LBX153LA log_LBX153LA 1.611002e-02
log_LBXF02LA log_LBXF02LA 1.332884e-02
log_LBXF07LA log_LBXF07LA 1.323860e-02
log_LBXF06LA log_LBXF06LA 1.167540e-02
log_LBXD01LA log_LBXD01LA 1.165293e-02
log_LBX199LA log_LBX199LA 9.822605e-03
log_LBXD07LA log_LBXD07LA 6.926303e-03
log_LBX194LA log_LBX194LA 1.571371e-03
log_LBX138LA log_LBX138LA 2.458828e-07
log_LBXD04LA log_LBXD04LA 1.714727e-07
log_LBX099LA log_LBX099LA 6.504292e-08
log_LBX187LA log_LBX187LA 5.880036e-08
log_LBXTCDLA log_LBXTCDLA 5.419637e-08
log_LBX105LA log_LBX105LA 3.803308e-08
log_LBX156LA log_LBX156LA 3.404966e-08
log_LBX157LA log_LBX157LA 2.928430e-08
log_LBXF08LA log_LBXF08LA 2.903891e-08
log_LBXF09LA log_LBXF09LA 1.954859e-08
log_LBX074LA log_LBX074LA 1.656317e-08
log_LBXD03LA log_LBXD03LA 1.393948e-08
log_LBXF05LA log_LBXF05LA 1.106955e-08
log_LBXD05LA log_LBXD05LA 1.053564e-08
log_LBXF03LA log_LBXF03LA 8.367659e-09

```r
summary(wqsmod$fit)
```


Call:
glm(formula = y ~ ., family = binomial(link = "logit"), data = new_data)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-0.8938  -0.8801  -0.8741   1.5017   1.5292  

Coefficients:
            Estimate Std. Error z value Pr(>|z|)  
(Intercept) -0.80027    0.35111  -2.279   0.0227 *
wqs          0.03074    0.21696   0.142   0.8873  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 376.12  on 299  degrees of freedom
Residual deviance: 376.10  on 298  degrees of freedom
AIC: 380.1

Number of Fisher Scoring iterations: 4

```r
gcompmod
```

Scaled effect size (positive direction, sum of positive coefficients = 1.69)
log_LBX156LA log_LBX194LA log_LBX105LA log_LBX138LA log_LBXF09LA 
    0.098960     0.094117     0.090794     0.073728     0.071542 
log_LBXHXCLA log_LBXD02LA log_LBXF08LA log_LBXD04LA log_LBXD01LA 
    0.071025     0.067497     0.066298     0.055614     0.052980 
log_LBXF05LA log_LBXPCBLA log_LBX118LA log_LBX153LA log_LBX189LA 
    0.051612     0.039568     0.037921     0.028775     0.026268 
log_LBXF01LA log_LBX157LA log_LBX180LA log_LBXF03LA log_LBXD05LA 
    0.024350     0.024323     0.013647     0.008858     0.001437 
log_LBX187LA 
    0.000685 

Scaled effect size (negative direction, sum of negative coefficients = -1.74)
log_LBX170LA log_LBXD07LA log_LBX199LA log_LBXF07LA log_LBX167LA 
      0.1666       0.1455       0.1442       0.1280       0.1041 
log_LBXF04LA log_LBXF06LA log_LBX099LA log_LBX196LA log_LBXF02LA 
      0.1021       0.0426       0.0406       0.0302       0.0264 
log_LBXD03LA log_LBX074LA log_LBXTCDLA 
      0.0247       0.0241       0.0210 

Mixture log(OR) (Delta method CI):
gamma (CI): -0.0581 (-0.472,0.356), z=-0.275, p=0.784

This model ...



|     x|
|-----:|
| 0.160|
| 0.147|
| 0.124|
| 0.113|
| 0.093|
| 0.071|
| 0.057|
| 0.045|
| 0.038|
| 0.036|
| 0.035|
| 0.019|
| 0.017|
| 0.012|
| 0.010|
| 0.007|
| 0.007|
| 0.006|



|     x|
|-----:|
| 0.215|
| 0.141|
| 0.097|
| 0.088|
| 0.087|
| 0.068|
| 0.063|
| 0.054|
| 0.043|
| 0.040|
| 0.038|
| 0.034|
| 0.019|
| 0.005|
| 0.005|
| 0.002|

In the second plot...

To see the estimates of the mixture effect, we run the following


```r
print(results)
```

This last table tells us...  




The `qgcomp.noboot` function gives back other outputs...

### Example 2

In the following code we run a logistic regression (`family = "binomial"`)...


```r
gcompmod = qgcomp.noboot(disease_state~sex+log_LBX074LA+log_LBX099LA+log_LBX105LA+log_LBX118LA+log_LBX138LA+log_LBX153LA+
                           log_LBX156LA+log_LBX157LA+log_LBX167LA+log_LBX170LA+log_LBX180LA+log_LBX187LA+log_LBX189LA+
                           log_LBX194LA+log_LBX196LA+log_LBX199LA+log_LBXD01LA+log_LBXD02LA+log_LBXD03LA+log_LBXD04LA+
                           log_LBXD05LA+log_LBXD07LA+log_LBXF01LA+log_LBXF02LA+log_LBXF03LA+log_LBXF04LA+log_LBXF05LA+
                           log_LBXF06LA+log_LBXF07LA+log_LBXF08LA+log_LBXF09LA+log_LBXPCBLA+log_LBXTCDLA+log_LBXHXCLA,
                         expnms=Xnm,
                         wqs_data[,c(Xnm, "sex", 'disease_state')], family=binomial(), q=4)
wqsmod = gwqs(disease_state ~ sex, mix_name = Xnm, data = wqs_data, q = 4, 
     validation = 0.6, b = 3, b1_pos = F, b1_constr = F, family='binomial', seed=125)
```

[1] "The optimization function did not converge 0 times"

From the first plot we see...


```r
print(gcompmod)
summary(gcompmod$fit)
```

                 mix_name  mean_weight
log_LBX167LA log_LBX167LA 1.773694e-01
log_LBXF01LA log_LBXF01LA 1.579708e-01
log_LBX180LA log_LBX180LA 1.408461e-01
log_LBX170LA log_LBX170LA 1.360615e-01
log_LBXF04LA log_LBXF04LA 1.153746e-01
log_LBXPCBLA log_LBXPCBLA 4.562147e-02
log_LBX196LA log_LBX196LA 4.351125e-02
log_LBXD02LA log_LBXD02LA 3.211508e-02
log_LBX118LA log_LBX118LA 3.147274e-02
log_LBXHXCLA log_LBXHXCLA 3.028686e-02
log_LBX189LA log_LBX189LA 2.805234e-02
log_LBX153LA log_LBX153LA 2.624729e-02
log_LBXF06LA log_LBXF06LA 1.656195e-02
log_LBXD01LA log_LBXD01LA 1.007162e-02
log_LBX199LA log_LBX199LA 3.430755e-03
log_LBXF02LA log_LBXF02LA 3.173771e-03
log_LBXD07LA log_LBXD07LA 1.831359e-03
log_LBX194LA log_LBX194LA 3.583818e-07
log_LBX138LA log_LBX138LA 1.417820e-07
log_LBX099LA log_LBX099LA 9.857947e-08
log_LBXF07LA log_LBXF07LA 9.418940e-08
log_LBXD04LA log_LBXD04LA 7.637178e-08
log_LBX187LA log_LBX187LA 6.021228e-08
log_LBXTCDLA log_LBXTCDLA 4.845691e-08
log_LBX156LA log_LBX156LA 4.835792e-08
log_LBXF08LA log_LBXF08LA 3.349906e-08
log_LBX105LA log_LBX105LA 3.113432e-08
log_LBX157LA log_LBX157LA 3.062864e-08
log_LBXF09LA log_LBXF09LA 2.353027e-08
log_LBX074LA log_LBX074LA 1.447614e-08
log_LBXD03LA log_LBXD03LA 1.282351e-08
log_LBXF05LA log_LBXF05LA 1.230270e-08
log_LBXF03LA log_LBXF03LA 1.102888e-08
log_LBXD05LA log_LBXD05LA 1.037097e-08

Call:
glm(formula = y ~ ., family = binomial(link = "logit"), data = new_data)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-0.9532  -0.9377  -0.8148   1.4305   1.6030  

Coefficients:
            Estimate Std. Error z value Pr(>|z|)   
(Intercept) -0.96314    0.36516  -2.638  0.00835 **
wqs          0.02236    0.21507   0.104  0.91721   
sex          0.34491    0.24915   1.384  0.16626   
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 376.12  on 299  degrees of freedom
Residual deviance: 374.16  on 297  degrees of freedom
AIC: 380.16

Number of Fisher Scoring iterations: 4


## References


## Acknowledgements

The development of this package was supported by NIH Grant RO1ES02953101
