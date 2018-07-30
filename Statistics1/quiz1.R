# Q8
library(DESeq2)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata = pData(bm)
edata = exprs(bm)

par(mfrow = c(1,2))
# log2 transform
mm = log2(edata[,1]+1) - log2(edata[,2]+1)
aa = log2(edata[,1]+1) + log2(edata[,2]+1)
plot(aa,mm,col=2)

# rlog transform
mm = rlog(edata)[,1] - rlog(edata)[,2]
aa = rlog(edata)[,1] + rlog(edata)[,2]
plot(aa,mm,col=2)


# Q9
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

library(rafalib)

dist1 = dist(t(edata))
dist2 = dist(t(dplyr::filter(edata,rowMeans(edata)>100)))
dist3 = dist(t(log2(edata+1)))
hclust1 = hclust(dist1)
hclust2 = hclust(dist2)
hclust3 = hclust(dist3)

par(mfrow=c(1,3))

myplclust(hclust1,labels=pdata$sample.id,lab.col = as.numeric(pdata$study),hang=0.1,xlab="",sub="")
myplclust(hclust2,labels=pdata$sample.id,lab.col = as.numeric(pdata$study),hang=0.1,xlab="",sub="")
myplclust(hclust3,labels=pdata$sample.id,lab.col = as.numeric(pdata$study),hang=0.1,xlab="",sub="")

# Q10
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

library(rafalib)
set.seed(1235)

# k-means clustering
k_mean = kmeans(t(log2(edata+1)),centers = 2)
#matplot(k_mean$centers)
# Hierarchical clustering
dist1 = dist(t(log2(edata+1)))
hclust1 = hclust(dist1)
plot(hclust1,hang = -1)

# Compare the results
par(mfrow=c(1,3))
myplclust(hclust1,labels=pdata$sample.id,lab.col = as.numeric(pdata$study),hang=0.1,xlab="",sub="")
hclust2 = cutree(hclust1,2)
myplclust(hclust1,labels=pdata$sample.id,lab.col = as.numeric(k_mean$cluster),hang=0.1,xlab="",sub="")
myplclust(hclust1,labels=pdata$sample.id,lab.col = as.numeric(hclust2),hang=0.1,xlab="",sub="")

# dend = as.dendrogram(hclust1)
# dend = color_labels(hclust1,2,1:2)