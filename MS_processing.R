library(oligo)
library(reshape2)
library(ggplot2)
library(arrayQualityMetrics)
library(limma)
library(pd.hugene.2.0.st)
library(miscTools)

## Read phenotype
pheno <- read.table("/Users/lro/Dropbox/Research/Projects/Ongoing/MS_microarray_RH/pheno.txt", fill = TRUE, sep="\t", header=TRUE)

setwd("~/Dropbox/Research/Projects/Ongoing/MS_microarray_RH/raw/")

## Read and normalize data
celFiles <- list.celfiles(full.names=TRUE)
matrix <- read.celfiles(celFiles)
matrix_rma <- rma(matrix)
matrix_exprs <- exprs(matrix_rma)

## QC
arrayQualityMetrics(matrix_rma, outdir = "../QC/")

## Check distribution
normdist <- melt(as.matrix(matrix_exprs))
ggplot(data=normdist, aes(x=value, color=Var2)) + 
  geom_density() +
  theme(legend.position = "none")
# remove probes with mean intensity below the peak of the distributions
matrix_exprs <- matrix_exprs[rowMedians(matrix_exprs) > 2.5,]

## Read probe to symbol translation table
probe2symbol <- read.csv("../HuGene-2_0-st-v1.na36.hg19.transcript.csv")
probe2symbol <- probe2symbol[,c(1,8)]
temp <- strsplit(as.character(probe2symbol[,2]), "//")
temp2 <- lapply(temp, function(x) x[2])
temp3 <- unlist(temp2)
temp3 <- gsub(" ", "", temp3)
probe2symbol[,2] <- temp3

## Remove probes without assignment
matrix_exprs <- matrix_exprs[rownames(matrix_exprs) %in% probe2symbol[complete.cases(probe2symbol),1],]
probe2symbol <- probe2symbol[complete.cases(probe2symbol),]

## Translate to gene symbols
rownames(matrix_exprs) <- probe2symbol[match(rownames(matrix_exprs), probe2symbol$transcript_cluster_id),2]

## Split matrix_exprs into 1:1 matches and 1:several
matrix_exprs_1 <- matrix_exprs[rownames(matrix_exprs) %in% names(which(table(probe2symbol$gene_assignment)==1)), ]
matrix_exprs_2 <- matrix_exprs[!rownames(matrix_exprs) %in% rownames(matrix_exprs_1),]

## Find largest median expression probes
matrix_exprs_2 <- matrix_exprs_2[order(rowMedians(matrix_exprs_2), decreasing = TRUE),]
matrix_exprs_2 <- matrix_exprs_2[!duplicated(rownames(matrix_exprs_2)),]

## Combine
matrix_exprs <- rbind(matrix_exprs_1, matrix_exprs_2)

## Save
save(matrix_exprs, pheno)