cat("sample splits testing")
library(qgcomp)

data(metals)
metals$clust = sample(seq_len(floor(nrow(metals)/8)), nrow(metals), replace=TRUE)

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

# are clusters allocated equally across training/testing?
margdist = as.numeric(prop.table(table(metals$clust))) # 70/30 split
# distance between marginal distribution of cluster and split specific clustering
#sqrt(sum((as.numeric(prop.table(table(spl$traindata$clust))) - margdist)^2)) #  invalid because this doesnt contain all clusters
#sqrt(sum((as.numeric(prop.table(table(spl$validdata$clust))) - margdist)^2)) # 0.04049317
# do all clusters show up in both datasets?
length(table(spl$traindata$clust)) - length(margdist)
length(table(spl$validdata$clust)) - length(margdist)


spl2 = split_data(metals, cluster="clust")
dim(spl2$traindata) # 181 observations = 40% of total
dim(spl2$validdata) # 271 observations = 60% of total
# distance between marginal distribution of cluster and split specific clustering

sqrt(sum((as.numeric(prop.table(table(spl2$traindata$clust))) - margdist)^2)) #  0.0116399
sqrt(sum((as.numeric(prop.table(table(spl2$validdata$clust))) - margdist)^2)) # 0.007774251
# do all clusters show up in both datasets?
stopifnot(length(table(spl2$traindata$clust)) == length(margdist))
stopifnot(length(table(spl2$validdata$clust)) == length(margdist))


