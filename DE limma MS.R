# Clearing global environment
rm(list = ls())

# -------------------------------------------------------------------------------------------------------------------
# Loading libraries

library(tidyverse)
library(ggplot2)
library(dplyr)
library(matrixTests)
library(broom)
library(purrr)
library(modelr)
library(ggpubr)
library(limma)
library(mtconnectR)

# Loading data into R
load("C:/Users/jeppe/Desktop/Bachelorprojekt/Data/MS data/data.RData")

# Renaming pheno labels
pheno <- pheno %>% 
  rename(status = R.RRMS..S.HC..G.Gilenya.behandlet.,
         gender = sex.male1,
         past_treatment = Tidligere.Beh,
         stability = X0Stabil_1ustabil)

# Sorting phenodata by Filename column
pheno <- pheno[order(pheno$Filename),]

### Do DE anaylysis for RRMS vs HC ###

# Choosing only the R and S columns
subset_columns <- grep("-R-|-S-", names(as.data.frame(matrix_exprs)), value = TRUE)

matrix_exprs_R_S <- matrix_exprs[, subset_columns]

# Subsetting the pheno data for only RRMS patients and healthy controls
pheno_R_S <- pheno %>%
  filter(status == "R" | status == "S")

# Creating factor for the searched levels
status_factor <- factor(pheno_R_S$status, levels = c("S", "R"))
gender_factor <- factor(pheno_R_S$gender, levels = c(1,0))

design <- model.matrix(~gender_factor+status_factor)

fit <- lmFit(matrix_exprs_R_S, design)

fit <- eBayes(fit, trend=TRUE, robust = TRUE)

results <- decideTests(fit)

summary(results)

# Differentially expressed genes for RRMS
toptable_R_S <- topTable(fit, coef = "status_factorR", n = 20)

# We can plot the fold-changes for the significant genes
plotMD(fit)

# Heatmap of the top 20 genes
colnames(matrix_exprs_R_S) <- pheno_R_S$status
matrix_exprs_R_S <- matrix_exprs_R_S[, order(colnames(matrix_exprs_R_S))]
heatmap(matrix_exprs_R_S[rownames(toptable_R_S), ], scale = "none")


### Do DE anaylysis for active vs inactive on fingolimod ###

# Choosing only the G columns
subset_columns_1 <- grep("-G-", names(as.data.frame(matrix_exprs)), value = TRUE)

matrix_exprs_G <- matrix_exprs[, subset_columns_1]

# Subsetting the pheno data for only MS patients on fingolimod treatment
pheno_G <- pheno %>%
  filter(status == "G")

# Creating factor for the searched levels - 0 = inactive MS, 1 = active mS
status_factor_1 <- factor(pheno_G$stability, levels = c(0, 1))
gender_factor_1 <- factor(pheno_G$gender, levels = c(1,0))

# Creating design matrix
design_1 <- model.matrix(~gender_factor_1+status_factor_1)

# Fit models
fit_1 <- lmFit(matrix_exprs_G, design_1)

fit_1 <- eBayes(fit_1, trend=TRUE, robust = TRUE)

results_1 <- decideTests(fit_1)

summary(results_1)

# Differentially expressed genes for active MS
toptable_G <- topTable(fit_1, coef = "status_factor_11", n = 20)

# We can plot the fold-changes for the significant genes
plotMD(fit_1)

# Heatmap of the top 20 genes
colnames(matrix_exprs_G) <- pheno_G$stability
matrix_exprs_G <- matrix_exprs_G[, order(colnames(matrix_exprs_G))]
heatmap(matrix_exprs_G[rownames(toptable_G), ], scale = "none")

