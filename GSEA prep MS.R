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

# Preparing data for ssgsea

# RRMS vs HC
# Creating matrix with only RRMS patients 
# Choosing only the R samples
subset_R_columns <- grep("-R-", names(as.data.frame(matrix_exprs)), value = TRUE)

matrix_exprs_R <- matrix_exprs[, subset_R_columns]

# Choosing only the S samples
subset_S_columns <- grep("-S-", names(as.data.frame(matrix_exprs)), value = TRUE)

matrix_exprs_S <- matrix_exprs[, subset_S_columns]

# Saving tab-separated files
write.table(matrix_exprs_R, file = "matrix_exprs_RRMS.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

write.table(matrix_exprs_S, file = "matrix_exprs_HC.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
