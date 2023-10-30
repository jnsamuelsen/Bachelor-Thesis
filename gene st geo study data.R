### gene st geo study
celpath <- "C:/Users/jeppe/Desktop/Bachelorprojekt/Data/CEL files/CELfiles1"

# import CEL files containing raw probe-level data into an R AffyBatch object
list = list.files(celpath,full.names=TRUE)
data = read.celfiles(list)

# Retrieving intensity data
expr <- exprs(data)
expr[1:5,]

# Looking at pm probes for row 1 to 5 (reordered)
pm(data)[1:5,]

# Looking at one specific probe set ID:
pm(data, "1056")

### Plotting intensities of all the probes in a probe set

# Retrieving sample annotation of the data
ph <- data@phenoData
ph

# Looking at the sample names
ph@data

pData(affyGeneFS)

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

# Now to the quality control #

# Adding shorter descriptive anntotation to the phenoData

ph@data[ ,1] <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25")
ph
data$index

# Creating microarray pictures
for (i in 1:26)
{
  name = paste("image",i,".jpg",sep="")
  jpeg(name)
  image(data[,i],main=ph@data$index[i])
  dev.off()
}


MAplot(data, pairs=TRUE,plot.method = "smoothScatter")

hist(data[,1:2])
