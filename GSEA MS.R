# Clearing global environment
rm(list = ls())

# -------------------------------------------------------------------------------------------------------------------
# Loading libraries
library(GSEABase)
library(singscore)

# Loading data
# Loading data into R
load("C:/Users/jeppe/Desktop/Bachelorprojekt/Data/MS data/data.RData")

# Loading gene sets
gene_sets <- getGmt("C:/Users/jeppe/Desktop/Bachelorprojekt/c5.bp.v7.1.symbols.gmt")

scoredf <- multiScore(rankData = rankGenes(ExpressionSet(matrix_exprs)), upSetColc = gene_sets, knownDirection = TRUE)

ExpressionSet(matrix_exprs), phenoData = AnnotatedDataFrame(pheno))
AnnotatedDataFrame(pheno)