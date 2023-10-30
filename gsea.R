# Clearing global environment
rm(list = ls())

# Loading R and Bioconductor packages
library(tidyverse)
library(oligo)
library(arrayQualityMetrics)
library(fgsea)

# -------------------------------------------------------------------------------------------------------
# Loading and preprocessing of the data
# -------------------------------------------------------------------------------------------------------

# Data is taken from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145120

# Path to CEL-files of whole blood ssc samples
celpath <- "C:/Users/jeppe/Desktop/Bachelorprojekt/Data/tmp data/Whole-blood SSC subsets"

# import CEL files containing raw probe-level data into an R AffyBatch object
list = list.files(celpath,full.names=TRUE)
data = read.celfiles(list)

# RMA on the raw data level
data_rma <- rma(data)

# Retrieving intensity data
data_exprs <- exprs(data_rma)

# Retrieving sample annotation of the data
ph <- data@phenoData
ph

# Looking at the sample names
ph@data

# Retrieving probe annotation of the data
feat <- data@featureData
feat
feat@data

# Retrieving experiment annotation
exp <- data@experimentData
exp

# Retrieving IDs of probe sets that are represented on the arrays
featureNames(data)

# Retrieving number of probe sets represented on the arrays
length(featureNames(data))

# Retrieving number of probes represented on the arrays
length(probeNames(data))

# Creating and adding shorter descriptive anntotation to the phenoData ("s" for sick individuals, "h" for healthy individuals)
new_pheno <- c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11", "s12", "s13", "s14", "s15", "s16", "s17", "s18", "s19", "s20", "s21", "s22", "s23", "s24", "s25", "s26", "s27", "s28", "s29", "s30", "s31", "s32", "s33", "s34", "s35", "h1", "h2", "h3", "h4", "h5", "h6", "h7", "h8", "h9", "h10", "h11", "h12", "s36", "s37", "s38", "s39", "s40", "s41", "s42", "s43", "s44", "s45", "s46", "s47", "s48", "s49", "s50", "s51", "s52", "s53", "s54", "s55", "s56", "s57", "s58", "s59", "s60", "s61", "s62", "s63", "s64", "s65", "s66", "s67", "s68", "s69", "s70", "s71", "h13", "s72", "s73", "s74", "s75", "s76", "s77", "s78", "s79", "s80", "s81", "s82", "s83", "s84", "s85", "s86", "s87", "s88", "s89", "s90", "s91", "s92", "s93", "s94", "s95", "s96", "s97", "s98", "s99", "s100", "s101", "s102", "s103", "s104", "s105", "h14", "s106", "s107", "s108", "s109", "s110", "h15", "s111", "s112", "s113", "s114", "s115", "s116", "s117", "s118", "s119", "s120", "s121", "s122", "s123", "h16", "h17", "h18", "h19", "h20", "h21", "s124", "h22", "h23", "h24", "h25", "h26", "h27", "h28", "h29", "s125", "s126", "s127", "s128", "s129", "s130", "s131", "s132", "s133", "s134", "s135", "s136", "s137", "s138", "s139", "s140", "s141", "h30", "h31", "s142", "h32", "h33", "h34", "h35", "h36", "h37", "h38", "s143", "h39", "h40", "s144", "s145", "s146", "s147", "s148", "s149", "s150")

ph@data[ ,1] <- new_pheno
ph@data

# -------------------------------------------------------------------------------------------------------
# Gene set enrichment analysis
# -------------------------------------------------------------------------------------------------------