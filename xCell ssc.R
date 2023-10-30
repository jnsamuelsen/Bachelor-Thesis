# Clearing global environment
rm(list = ls())

# -------------------------------------------------------------------------------------------------------------------
# Loading libraries

library(tidyverse)
library(ggplot2)
library(dplyr)
library(xCell)

# Loading data into R
# data <- read.table("C:/Users/jeppe/Desktop/Bachelorprojekt/CIBERSORT/Mixture files/whole_blood_ssc_mixture_file.txt", sep = "\t", header = TRUE)

exprMatrix = read.table("C:/Users/jeppe/Desktop/Bachelorprojekt/CIBERSORT/Mixture files/whole_blood_ssc_mixture_file.txt",header=TRUE,row.names=1, as.is=TRUE)

xCellAnalysis(exprMatrix)
