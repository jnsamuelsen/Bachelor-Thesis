library(affy)

# Loading 10 sample subset gliablastoma microarray data matrix
load(file = "TCGA10SampleSubset.rda")

# Quality control on probe level data

MAplot(TCGA10SampleSubset,pairs=TRUE,plot.method="smoothScatter")

Index <- c(1,2,3,100,1000,2000)
pm(TCGA10SampleSubset)[Index,]

# Quality control on expression level data

affy::hist(TCGA10SampleSubset[,1:2])

par(mfrow=c(2,2))
image(TCGA10SampleSubset)

par(mfrow=c(1,1))
boxplot(TCGA10SampleSubset, col=c(2,3,4))

# Analyzing RNA degradation plots

deg <- AffyRNAdeg(TCGA10SampleSubset)
