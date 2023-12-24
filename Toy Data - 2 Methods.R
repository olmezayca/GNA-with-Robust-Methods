# Choose a set of soft-thresholding powers

powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

### NETWORK CONSTRUCTIONS

softPower = 8;

#################### WGCNA

adjacency = adjacency(datExpr, power = softPower)
TOM = TOMsimilarity(adjacency)
wgcna = 1-TOM


#################### BICOR

adjacency = abs(bicor(datExpr, use="pairwise.complete.obs"))^softPower
TOM = TOMsimilarity(adjacency)
bicor = 1-TOM


#################### PLSR

adjacency = abs(plsr_s_ik)^softPower
TOM = TOMsimilarity(adjacency)
plsr = 1-TOM


#################### PRM

adjacency = abs(prm_s_ik)^softPower
TOM = TOMsimilarity(adjacency)
prm = 1-TOM


