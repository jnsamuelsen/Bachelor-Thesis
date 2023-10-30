### Pre-processing data

library(oligo)
library(reshape2)
library(ggplot2)
library(arrayQualityMetrics)
library(limma)
library(pd.hugene.2.0.st)
library(miscTools)
library(sva)
library(gridExtra)
library(xCell)

## Read phenotype
pheno <- read.table("/Users/lro/Dropbox/Research/Projects/Ongoing/MS_microarray_RH/pheno.txt", fill = TRUE, sep="\t", header=TRUE)

setwd("~/Dropbox/Research/Projects/Ongoing/MS_microarray_RH/raw/")

## Read and normalize data
celFiles <- list.celfiles(full.names=TRUE)
matrix <- read.celfiles(celFiles)
matrix_rma <- rma(matrix)
matrix_exprs <- exprs(matrix_rma)
matrix_exprs <- matrix_exprs[,pheno$Filename]

## QC
arrayQualityMetrics(matrix_rma, outdir = "../QC/")

## Check for batch effects and other confounders
pca <- prcomp(t(matrix_exprs), scale=TRUE)
df <- data.frame(pc1=pca$x[,1], pc2=pca$x[,2],batch=as.factor(pheno$Affy.Run_nr), age=pheno$Age, treatment=pheno$R.RRMS..S.HC..G.Gilenya.behandlet., previous_treatment=pheno$Tidligere.Beh)

ggplot(data=df, aes(x=pc1, y=pc2, color=batch)) +
  geom_point()

ggplot(data=df, aes(x=pc1, y=pc2, color=age)) +
  geom_point()

ggplot(data=df, aes(x=pc1, y=pc2, color=treatment)) +
  geom_point()

ggplot(data=df, aes(x=pc1, y=pc2, color=previous_treatment)) +
  geom_point()
# Batch correction not possible

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

## Adjust for sex
matrix_exprs_sex_adj <- NULL
for(j in 1:nrow(matrix_exprs)) {
  print(j)
  fit <- lm(matrix_exprs[j,] ~ pheno$sex.male1)
  matrix_exprs_sex_adj <- rbind(matrix_exprs_sex_adj, fit$residuals + fit$coefficients["(Intercept)"])
}
rownames(matrix_exprs_sex_adj) <- rownames(matrix_exprs); colnames(matrix_exprs_sex_adj) <- colnames(matrix_exprs)

## Adjust for sex and age
matrix_exprs_sex_age_adj <- NULL
for(j in 1:nrow(matrix_exprs)) {
  print(j)
  fit <- lm(matrix_exprs[j,] ~ pheno$sex.male1)
  adj <- fit$residuals + fit$coefficients["(Intercept)"]
  fit <- lm(adj ~ pheno$Age)
  matrix_exprs_sex_age_adj <- rbind(matrix_exprs_sex_age_adj, fit$residuals + fit$coefficients["(Intercept)"])
}
rownames(matrix_exprs_sex_age_adj) <- rownames(matrix_exprs); colnames(matrix_exprs_sex_age_adj) <- colnames(matrix_exprs)


## Find most severe example of sex adjustment
r2s <- c()
cor <- c()
slopes <- c()
for(j in 1:nrow(matrix_exprs)) {
  fit <- lm(matrix_exprs_sex_adj[j,] ~ matrix_exprs[j,])
  r2s <- append(r2s, summary(fit)$adj.r.squared)
  cor <- append(cor, cor(matrix_exprs_sex_adj[j,], matrix_exprs[j,]))
  slopes <- append(slopes, fit$coef[[2]])
}
rownames(matrix_exprs)[which.min(r2s)]
rownames(matrix_exprs)[which.min(cor)]
rownames(matrix_exprs)[which.max(1-abs(slopes))]
gene <- which.min(r2s)


## Plot most severe case
fit <- lm(matrix_exprs[gene,] ~ pheno$sex.male1)
data <- data.frame(sex=as.factor(pheno$sex.male1), expr=matrix_exprs[gene,])
g1 <- ggplot(data, aes(x=sex, y=expr, group=sex)) +
  geom_boxplot(aes(color=sex)) +
  annotate(geom="text", x=2, y=max(data$expr), label = paste0("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), "\nSlope = ",signif(fit$coef[[2]], 5)), hjust=0) +
  ggtitle(paste0("sex vs expr of ", rownames(matrix_exprs)[gene]))

fit <- lm(matrix_exprs_sex_adj[gene,] ~ matrix_exprs[gene,])
data <- data.frame(expr_adj=matrix_exprs_sex_age_adj[gene,], expr=matrix_exprs[gene,], sex=as.factor(pheno$sex.male1))
g2 <- ggplot(data, aes(x=expr_adj, y=expr)) +
  geom_point(aes(color=sex)) +
  geom_smooth(method=lm) +
  annotate(geom="text", x=max(data$expr_adj), y=max(data$expr), label = paste0("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), "\nSlope = ",signif(fit$coef[[2]], 5)), hjust=1) +
  ggtitle(paste0("expr vs sex adjusted expr of ", rownames(matrix_exprs_sex_age_adj)[gene]))

fit <- lm(matrix_exprs_sex_adj[gene,] ~ pheno$sex.male1)
data <- data.frame(sex=as.factor(pheno$sex.male1), expr=matrix_exprs_sex_adj[gene,])
g3 <- ggplot(data, aes(x=sex, y=expr, group=sex)) +
  geom_boxplot(aes(color=sex)) +
  # geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  # geom_smooth(method=lm) +
  annotate(geom="text", x=2, y=max(data$expr), label = paste0("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), "\nSlope = ",signif(fit$coef[[2]], 5)), hjust=0) +
  ggtitle(paste0("sex vs sex adjusted expr of ", rownames(matrix_exprs)[gene]))

grid.arrange(g1,g2,g3, ncol=3, nrow=1)


## Find most severe example of age adjustment
r2s <- c()
cor <- c()
slopes <- c()
for(j in 1:nrow(matrix_exprs)) {
  fit <- lm(matrix_exprs_sex_age_adj[j,] ~ matrix_exprs_sex_adj[j,])
  r2s <- append(r2s, summary(fit)$adj.r.squared)
  cor <- append(cor, cor(matrix_exprs_sex_age_adj[j,], matrix_exprs_sex_adj[j,]))
  slopes <- append(slopes, fit$coef[[2]])
}
rownames(matrix_exprs)[which.min(r2s)]
rownames(matrix_exprs)[which.min(cor)]
rownames(matrix_exprs)[which.max(1-abs(slopes))]
gene <- which.min(r2s)


## Plot most severe case of age adjustment
fit <- lm(matrix_exprs_sex_adj[gene,] ~ pheno$Age)
data <- data.frame(age=pheno$Age, expr=matrix_exprs_sex_adj[gene,])
g1 <- ggplot(data, aes(x=age, y=expr)) +
  geom_point() +
  geom_smooth(method=lm) +
  annotate(geom="text", x=max(data$age), y=max(data$expr), label = paste0("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), "\nSlope = ",signif(fit$coef[[2]], 5)), hjust=1) +
  ggtitle(paste0("age vs expr of ", rownames(matrix_exprs)[gene]))

fit <- lm(matrix_exprs_sex_age_adj[gene,] ~ matrix_exprs_sex_adj[gene,])
data <- data.frame(expr_adj=matrix_exprs_sex_age_adj[gene,], expr=matrix_exprs_sex_adj[gene,])
g2 <- ggplot(data, aes(x=expr_adj, y=expr)) +
  geom_point() +
  geom_smooth(method=lm) +
  annotate(geom="text", x=max(data$expr_adj), y=max(data$expr), label = paste0("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), "\nSlope = ",signif(fit$coef[[2]], 5)), hjust=1) +
  ggtitle(paste0("expr vs age adjusted expr of ", rownames(matrix_exprs_sex_age_adj)[gene]))

fit <- lm(matrix_exprs_sex_age_adj[gene,] ~ pheno$Age)
data <- data.frame(age=pheno$Age, expr=matrix_exprs_sex_age_adj[gene,])
g3 <- ggplot(data, aes(x=age, y=expr)) +
  geom_point() +
  geom_smooth(method=lm) +
  annotate(geom="text", x=max(data$age), y=max(data$expr), label = paste0("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), "\nSlope = ",signif(fit$coef[[2]], 5)), hjust=1) +
  ggtitle(paste0("age vs age adjusted expr of ", rownames(matrix_exprs)[gene]))

grid.arrange(g1,g2,g3, ncol=3, nrow=1)


### Transcriptomic immunophenotyping analysis

## xCell on unadjusted
xCell_matrix <- xCellAnalysis(matrix_exprs)

## xCell on age adjusted genes
xCell_matrix_adj <- xCellAnalysis(matrix_exprs_sex_age_adj)

## Save data
save(matrix_exprs, matrix_exprs_sex_adj, matrix_exprs_sex_age_adj, xCell_matrix, xCell_matrix_adj, pheno, file="/Users/lro/Dropbox/Research/Projects/Ongoing/MS_microarray_RH/data.Rdata")


### Statistical tests for difference in phenotypes
## RRMS vs HC controlling for age+sex


## Active vs inactive MS/fingolimod controlling for age+sex


## RRMS male vs female (not controling for age)



