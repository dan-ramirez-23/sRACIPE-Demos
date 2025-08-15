library(sRACIPE)
library(ggplot2)

## Here we use some synthetic circuits to compare how distributions of features
## changes when depending on whether Transcription Factor (TF) or Protein 
## Degradation (PD) regulation occurs.
## First we compare the robustness of bistability in the toggle switch when the
## inhibition is either via TF or PD regulation. To begin, we pull both forms of
## the circuit from the topos file in this directory.

toggleSwitch_tf <- read.table(file.path("topos",paste0("toggleSwitch",".tpo")), header=T)
toggleSwitch_pd <- read.table(file.path("topos",paste0("toggleSwitch_pd",".tpo")), header=T)


## Next we run simulations of both to get the steady states.

ts_tfSet <- sracipeSimulate(circuit = toggleSwitch_tf, numModels = 1000, nIC = 5,
                            uniqueDigits = 3)
ts_pdSet <- sracipeSimulate(circuit = toggleSwitch_pd, numModels = 1000, nIC = 5,
                            uniqueDigits = 3)

## The number of unique steady states per model in the ensemble is automatically
## stored in the output metadata when nIC > 1

ts_tfStateCounts <- ts_tfSet@metadata[["uniqueStateCounts"]][["UniqueStableStateNo"]]
ts_pdStateCounts <- ts_pdSet@metadata[["uniqueStateCounts"]][["UniqueStableStateNo"]]

## Finally, we can make a bar plot to compare the proportion of bistability in 
## each model

ts_tfFraction <- sum(ts_tfStateCounts == 2) / length(ts_tfStateCounts)
ts_pdFraction <- sum(ts_pdStateCounts == 2) / length(ts_pdStateCounts)

tsFractions <- c(TF = ts_tfFraction, PD = ts_pdFraction)

barplot(tsFractions,
        ylim = c(0, 0.2),
        ylab = "Bistable Fraction",
        main = "Proportion of Bistability in Different TS Models",
        col = "skyblue")


## We can do a further comparison by comparing clustering data for each model. 
## We can use the sracipeUniqueStates() function to only examine unique steady 
## states in our ensembles by filtering out repeats

ts_tf_SS <- sracipeUniqueStates(ts_tfSet)
ts_pd_SS <- sracipeUniqueStates(ts_pdSet)

## sracipeUniqueStates() returns the steady states in a list so that we can still
## connect them to their original model parameter set, but to use the states for
## clustering, we need to cbind the list into a single dataframe, while filtering
## out the models with no steady states

ts_tf_SS <- Filter(function(x) nrow(x) > 0 && ncol(x) > 0, ts_tf_SS)
ts_tf_SS <- t(do.call(cbind, ts_tf_SS))

ts_pd_SS <- Filter(function(x) nrow(x) > 0 && ncol(x) > 0, ts_pd_SS)
ts_pd_SS <- t(do.call(cbind, ts_pd_SS))


## Next we perform normalization on the data by log-transforming the states and 
## normalizing by mean and standard deviation by gene. Because we want to compare the 
## resulting clusters, it is important to normalize by the same mean and standard
## deviation for both data sets. Here, we have chosen to normalize the PD data
## using the TF means and standard deviations, which "projects" the PD data onto
## the TF space.

tmpMeans <- rowMeans(t(log2(1+ts_tf_SS)))
tmpSds <- apply(t(log2(1+ts_tf_SS)),1,sd)
tf_SS_norm <- log2(1+ts_tf_SS)
tf_SS_norm <- sweep(tf_SS_norm, 2, tmpMeans, FUN = "-") # scale
tf_SS_norm <- sweep(tf_SS_norm, 2, tmpSds, FUN = "/") # scale

pd_SS_norm <- log2(1+ts_pd_SS) # Log transform
pd_SS_norm <- sweep(ts_pd_norm, 2, tmpMeans, FUN = "-") # scale
pd_SS_norm <- sweep(ts_pd_norm, 2, tmpSds, FUN = "/") # scale

## Next, PCA is performed on the TF data, and we continue to project the PD data
## into the PCA space for the TF data
  
pca_tf <- prcomp(tf_SS_norm)
pca_pd <- scale(ts_pd_norm, pca_tf$center, pca_tf$scale) %*% pca_tf$rotation


## Now that we have finished preparing our data for comparison, we can perform
## hierarchical clustering on each data set. Ideally, one should use a statistical
## method to determine the ideal number of clusters, but in this simple case, 
## we know that there should be roughly two clusters.

k_use <- 2
dist_mat_tf <- dist(pca_tf$x)
hc_tf <- hclust(dist_mat_tf, method = "ward.D2")
cluster_labels_tf <- cutree(hc_tf, k = k_use)

dist_mat_pd <- dist(pca_pd)
hc_pd <- hclust(dist_mat_pd, method = "ward.D2")
cluster_labels_pd <- cutree(hc_pd, k = k_use)



## The last step is plot the results of the clustering, first as projected onto
## the PCA space, and then as heat maps
pc1_weight <- round(100*summary(pca_tf)$importance[2,1],2)
pc2_weight <- round(100*summary(pca_tf)$importance[2,2],2)
plot_xlab <- paste("PC1 (",pc1_weight,"%)",sep="")
plot_ylab <- paste("PC2 (",pc2_weight,"%)",sep="")

image_tf <- ggplot(pca_tf$x, aes(x=PC1, y=PC2, color=as.factor(cluster_labels_tf))) +
  geom_point(size=2) +
  scale_color_manual(values=c("blue", "red")) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  labs(color="Cluster") + 
  ggtitle("Toggle Switch TF") +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        panel.background = element_rect("white"),
        plot.background = element_rect("white"),
        axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        plot.title = element_text(size = 22, hjust = 0.5))
image_tf


image_pd <- ggplot(pca_pd, aes(x=PC1, y=PC2, color=as.factor(cluster_labels_pd))) +
  geom_point(size=2) +
  scale_color_manual(values=c("blue", "red")) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  labs(color="Cluster") + 
  ggtitle("Toggle Switch PD") +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        panel.background = element_rect("white"),
        plot.background = element_rect("white"),
        axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        plot.title = element_text(size = 22, hjust = 0.5))
image_pd



