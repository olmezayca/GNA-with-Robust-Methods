
install.packages("pcaPP")
library(pcaPP)
install.packages("chemometrics")
library(chemometrics)
install.packages("BiocManager")
BiocManager::install("WGCNA")
library(WGCNA)
install.packages("flashClust")
library(flashClust)


#Datasets:

getwd()
workingDir = "" # Enter your working directory
setwd(workingDir)
options(stringsAsFactors = FALSE)

femData = read.csv("LiverFemale3600.csv")

datExpr0 = as.data.frame(t(femData[, -c(1:8)]))
names(datExpr0) = femData$substanceBXH
rownames(datExpr0) = names(femData)[-c(1:8)]

# If the statement returns TRUE, there were no missing values in samples.
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK


# Outlier detection
sampleTree = hclust(dist(datExpr0), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 15, col = "red");

row.names(datExpr0)  #F2_221 is the 81.sample

#Final dataset
set.seed(5)
sample<-sample(1:3600,200)


datExpr = datExpr0[80:95,sample]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#########
traitData = read.csv("ClinicalTraits.csv");
allTraits = traitData[, -c(31, 16)]
allTraits = allTraits[, c(2, 11:36) ]

femaleSamples = rownames(datExpr)
traitRows = match(femaleSamples, allTraits$Mice)
datTraits = allTraits[traitRows, -1]       #clinical traits
rownames(datTraits) = allTraits[traitRows, 1]
collectGarbage()

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

save(datExpr, datTraits, file = "datExpr_toydata.RData")


