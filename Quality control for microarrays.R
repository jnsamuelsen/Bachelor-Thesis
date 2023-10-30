# Loading R-packages
library(Matrix)
library(lattice)
library(fdrtool)
library(rpart)
library(ggplot2)

# Installing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

BiocManager::install(c("oligo", "limma", "Biobase", "Biostrings", "genefilter", "ath1121501cdf"))

# Loading Bioconductor packages
library(oligo)
library(limma)
library(Biobase)
library(Biostrings)
library(genefilter)
library(ath1121501cdf)

#--------------------------------------------------------------------------------------------------------------------------

# Look at the data

# Retrieve intensities using oligo

# Retrieve sample annotation of the data using oligo

# Retrieve probe annotation of the data using oligo

# Retrieve experiment annotation of the data using oligo

# Retrieve other information about the data?


# Quality control of microarray data

# Add meaningful annotation to the phenoData


# Plots

# Microarray pictures and Chip pseudo-images - can show large inconsistencies on individual arrays

# Histograms - plot the distribution of log base 2 intensities (log2(PMij) for array i and probe j) of perfect match probes for comparison of probe intensity behavior between different arrays. If you see differences in shape or center of the distributions, it means that normalization is required

# Box plots - boxplots and histograms show the same differences in probe intensity behavior between arrays. In order to perform meaningful statistical analysis and inferences from the data, you need to ensure that all the samples are comparable. To examine and compare the overall distribution of log transformed PM intensities between the samples you can use a histogram but you will get a clearer view with a box plot

# MA plots - people started using them to compare each Affymetrix array to a pseudo-array. The pseudo array consists of the median intensity of each probe over all arrays
# The MA plot shows to what extent the variability in expression depends on the expression level (more variation on high expression values?)


# Calculate quality measures to assess the quality of the arrays

# Normalization - systematic differences between the samples that are due to noise rather than true biological variability should be removed in order to make biologically meaningfull conclusions about the data
# Normalization using RMA
# Checking the effect of the background correction by plotting raw versus background corrected data, can also be compared using boxplots or MA plots

# PCA plot - To check whether the overall variability of the samples reflects their grouping, we can perform a Principal Component Analysis. In other words, we can check if replicates are homogenous and distinguishable from samples of (the) other group(s)

# Other sanity checks?!


#--------------------------------------------------------------------------------------------------------------------------

# A workled example of quality control of gene ST affymetrix arrays (using the oligo package)

# Loading the affyGeneFS objects through the oligoData package

# Installing and loading the oligoData package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("oligoData")

library(oligoData)

# First some preprocessing work #

# Loading affyGeneFS feature set
data(affyGeneFS)
affyGeneFS

# RMA at probeset level
genePS <- rma(affyGeneFS, target = "probeset")

# RMA at transcript (expression) level
geneCore <- rma(affyGeneFS, target = "core")

# Retrieving intensity data
expr <- exprs(affyGeneFS)
expr[1:5,]

# Looking at pm probes for row 1 to 5 (reordered)
pm(affyGeneFS)[1:5,]

### Plotting intensities of all the probes in a probe set

# Retrieving sample annotation of the data
ph <- affyGeneFS@phenoData
ph

# Looking at the sample names
ph@data$index

pData(affyGeneFS)

# Retrieving probe annotation of the data
feat <- affyGeneFS@featureData
feat
feat@data

# Retrieving experiment annotation
exp <- affyGeneFS@experimentData
exp

# Retrieving IDs of probe sets that are represented on the arrays
featureNames(affyGeneFS)

# Retrieving number of probe sets represented on the arrays
length(featureNames(affyGeneFS))

# Retrieving number of probes represented on the arrays
length(probeNames(affyGeneFS))

### Now to the quality control #################################################

## Adding shorter descriptive anntotation to the phenoData ##

ph@data[ ,1] <- c("Brain_1_WT", "Brain_2_WT", "Brain_3_WT", "Breast_1_WT", "Breast_2_WT", "Breast_3_WT", "Heart_1_WT", "Heart_2_WT", "Heart_3_WT", "Kidney_1_WT", "Kidney_2_WT", "Kidney_3_WT", "Liver_1_WT", "Liver_2_WT", "Liver_3_WT", "Pancreas_1_WT", "Pancreas_2_WT", "Pancreas_3_WT", "Prostate_1_WT", "Prostate_2_WT", "Prostate_3_WT", "SkMus_1_WT", "SkMus_2_WT", "SkMus_3_WT", "Spleen_1_WT", "Spleen_2_WT", "Spleen_3_WT", "Testis_1_WT", "Testis_2_WT", "Testis_3_WT", "Thyroid_1_WT", "Thyroid_2_WT", "Thyroid_3_WT")
ph@data$index

## Creating microarray pictures ##
for (i in 1:33)
{
  name = paste("image",i,".jpg",sep="")
  jpeg(name)
  image(affyGeneFS[,i],main=ph@data$index[i])
  dev.off()
}

## Creating chip pseudo-images ##

# First, fitting a probe-level model to the data that assumes that all probes of a probe set behave the same in the different sample
# It fits robust probe level linear models to all the probe sets in an FeatureSet
Pset = fitProbeLevelModel(affyGeneFS)

# Creating pseudo-images based on weights (can also be done based on residuals) for the first 5 samples
for (i in 1:5)
{
  name = paste("pseudoimage",i,".jpg",sep="")
  jpeg(name)
  image(Pset,which=i,main=ph@data$index[i])
  dev.off()
}

## Now we're gong to create a histogram for the samples using ggplot ##

# First, the data should get in the correct format for ggplot: a data frame with log intensities in one column and sample names in the second column

# Using only the PM intensities
pmexp <- pm(affyGeneFS)

# Creating the data frame
sampleNames = vector()
logs = vector()
for (i in 1:33)
{
  sampleNames = c(sampleNames,rep(ph@data[i,1],dim(pmexp)[1]))
  logs = c(logs,log2(pmexp[,i]))
}

logData <- data.frame(logInt=logs,sampleName=sampleNames)

# Plotting histogram with ggplot
dataHist2 <- ggplot(logData, aes(logInt, colour = sampleName)) 
dataHist2 <- dataHist2 + geom_density()


# Plotting histogram with oligo on normalized data
color=c('green','green','green','red','red','red')
oligo::hist(geneCore[,1:33],lwd=2,col=color,ylab='Density',xlab='Log2 intensities',main='Histogram of raw data')

# Plotting boxplot with oligo on normalized data
name = "boxplot.jpg"
jpeg(name)
boxplot(geneCore,which='pm',col='red',names=ph@data$index) 
dev.off()

# Creating MA plots for the first 5 arrays in the data set

for (i in 1:5)
{
  name = paste("MAplot",i,".jpg",sep="")
  jpeg(name)
  MAplot(geneCore,which=i)
  dev.off()
}


# Creating PCA plot of data
color=c('green','green','green','red','red','red','blue','blue','blue', 'purple', 'purple', 'purple', 'gray', 'gray', 'gray','black', 'black', 'black', 'brown', 'brown', 'brown', 'yellow', 'yellow', 'yellow','magenta', 'magenta', 'magenta', 'orange', 'orange', 'orange','violetred', 'violetred', 'violetred')
data.PC = prcomp(t(expr),scale.=TRUE)
plot(data.PC$x,col=color)



