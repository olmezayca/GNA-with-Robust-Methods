
#####NA 

is.na(datExpr)

sum(is.na(datExpr))

for (j in 1:ncol(datExpr)) {
  for (i in 1:nrow(datExpr)) {
    if (is.na(datExpr[i,j])!= F) {
      datExpr[i,j] <- median(datExpr[-i,j], na.rm = T)
    }
  }
}


#### PCA
length(datExpr)

library("factoextra")
library("FactoMineR")

pca <- prcomp(t(datExpr))
eigs <- pca$sdev^2
eigs[1] / sum(eigs)
summary(pca)

rbind(
  SD = sqrt(eigs),
  Proportion = eigs/sum(eigs),
  Cumulative = cumsum(eigs)/sum(eigs))


res.pca <- PCA(t(datExpr), graph = FALSE)
eig.val <- get_eigenvalue(pca)
eig.val
dev.new()
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
dev.off()
collectGarbage()

Cumulative = cumsum(eigs)/sum(eigs)
plot(1:30,Cumulative[1:30])



pca_result <- prcomp(t(datExpr), scale. = TRUE)

# Calculate cumulative explained variance
cumulative_variance <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

# Plot cumulative explained variance
plot(cumulative_variance, type = "b", xlab = "Number of Components", 
     ylab = "Cumulative Explained Variance", 
     main = "Cumulative Explained Variance by PCA Components",
     ylim = c(0, 100))

# Add a horizontal line at 95% for reference
abline(h = 95, col = "red", lty = 2)

index_95_perce <- which(cumulative_variance >= 95)[1]
index_95_perce

