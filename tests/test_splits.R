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
  traindata=spl$traindata,validdata=spl$validdata, expnms=Xnm, .fixbreaks = FALSE, .globalbreaks = TRUE)
splitres

# check for break preservation
posbr = splitres$pos.fit$breaks[[1]]
posnm = splitres$pos.fit$expnms[[1]]
negbr = splitres$neg.fit$breaks[[1]]
negnm = splitres$neg.fit$expnms[[1]]
posidx = which(splitres$train.fit$expnms == posnm)
negidx = which(splitres$train.fit$expnms == negnm)
stopifnot(all.equal(splitres$train.fit$breaks[[posidx]], posbr))
stopifnot(all.equal(splitres$train.fit$breaks[[negidx]], negbr))

splitres2 <- qgcomp.partials(fun="qgcomp.noboot", f=y~., q=4, 
                            traindata=spl$traindata,validdata=spl$validdata, expnms=Xnm, .fixbreaks = TRUE, .globalbreaks = FALSE)
splitres2

# check for break preservation
posbr2 = splitres2$pos.fit$breaks[[1]]
posnm2 = splitres2$pos.fit$expnms[[1]]
negbr2 = splitres2$neg.fit$breaks[[1]]
negnm2 = splitres2$neg.fit$expnms[[1]]
posidx2 = which(splitres2$train.fit$expnms == posnm)
negidx2 = which(splitres2$train.fit$expnms == negnm)
stopifnot(all.equal(splitres2$train.fit$breaks[[posidx2]], posbr2))
stopifnot(all.equal(splitres2$train.fit$breaks[[negidx2]], negbr2))


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
length(table(spl2$traindata$clust)) - length(margdist)
length(table(spl2$validdata$clust)) - length(margdist)


