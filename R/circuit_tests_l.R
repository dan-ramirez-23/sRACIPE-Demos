rm(list=ls())
library(sRACIPE)
library(ggplot2)


topo <- data.frame(Source=c("So", "St", "Lo", "Lo", "De", "So", "De"),
                   Target=c("St", "Lo", "St", "De", "Lo", "De", "So"),
                   Type=c(2, 1, 1, 2, 2, 1, 1))


racipe <- sracipeSimulate(circuit = topo, numModels = 100, nIC = 100)


racipeNorm <- sracipeNormalize(racipe)

exprMat <- as.data.frame(t(assay(racipeNorm)))

pca <- prcomp(exprMat)


pca_df <- pca$x


ggplot(pca_df, aes(x=PC1, y=PC2)) +
  geom_point(size=3, aes(color=exprMat$Lo))
ggplot(pca_df, aes(x=PC1, y=PC2)) +
  geom_point(size=3, aes(color=exprMat$De))
ggplot(pca_df, aes(x=PC1, y=PC2)) +
  geom_point(size=3, aes(color=exprMat$So))
ggplot(pca_df, aes(x=PC1, y=PC2)) +
  geom_point(size=3, aes(color=exprMat$St))
