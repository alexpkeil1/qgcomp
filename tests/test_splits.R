cat("sample splits testing")
library(qgcomp)

data(metals)
set.seed(1231124)
spl = split_data(metals)
Xnm <- c(
  'arsenic','barium','cadmium','calcium','chromium','copper',
 'iron','lead','magnesium','manganese','mercury','selenium','silver',
 'sodium','zinc'
)
dim(spl$traindata) # 181 observations = 40% of total
dim(spl$validdata) # 271 observations = 60% of total
splitres <- qgcomp.partials(fun="qgcomp.noboot", f=y~., q=4, 
  traindata=spl$traindata,validdata=spl$validdata, expnms=Xnm)
splitres

spl = split_data(metals)
