##### SETUP & GLOBALS #####
rm(list=ls())
#devtools::install_github("lusystemsbio/sRACIPE")
library(sRACIPE)
library(ggplot2)
library(microbenchmark)

# seed for reproducibility & color palette for plots
set.seed(1234)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(3,2,4:8)]


topo_list <- c("toggleSwitch", "CTS", "TS_rep_2", "cellcycle")
k_list <- c(2, 4, 3, 6) # expected number of clusters for each topo
names(k_list) <- topo_list

##### SIMULATION & EDA ON UNPERTURBED NETWORK #####
for(topo_name in topo_list) {
  ## SETUP
  topoDir <- file.path(getwd(),topo_name)
  plotDir <- file.path(topoDir,"plots")
  dataDir <- file.path(topoDir,"data")

  if(!dir.exists(topoDir)) {
    dir.create(topoDir)
  }
  if(!dir.exists(dataDir)) {
    dir.create(dataDir)
  }
  if(!dir.exists(plotDir)) {
    dir.create(plotDir)
  }

  topo <- read.table(file.path("topos",paste0(topo_name,".tpo")), header=T)


  ## SIMULATION (w/ time cost benchmarking)
  mb_rs_fname <- file.path(dataDir,"microbenchmark_simulation_results.Rds")
  if(!file.exists(mb_rs_fname)) {
    mb_rs <- microbenchmark({
      racipe <- sracipeSimulate(topo, numModels = 1000, nIC = 100)
    }, times = 10)
    saveRDS(mb_rs, mb_rs_fname)
    saveRDS(racipe, racipe_fname)
  } else {
    mb_rs <- readRDS(mb_rs_fname)
    racipe <- readRDS(racipe_fname)
  }


  ## GET UNIQUE STEADY STATES

  ## NORMALIZATION
  genes <- rownames(racipe)
  unnormData <- t(assay(racipe))
  simExp <- assay(racipe, 1)
  simExp <- log2(1+simExp)
  tmpMeans <- rowMeans(simExp)
  tmpSds <- apply(simExp,1,sd)

  racipeNorm <- sracipeNormalize(racipe)
  racipeData <- as.data.frame(t(assay(racipeNorm)))
  exprMat_norm <- racipeData


  ## PCA
  pca_fname <- file.path(dataDir,"pca.Rds")
  if(!file.exists(pca_fname)) {
    saveRDS(pca, pca_fname)
  } else {
    pca <- readRDS(pca_fname)
  }

  ## CLUSTERING
  clust_fname <- file.path(dataDir,"cluster_labels.Rds")
  if(!file.exists(clust_fname)) {
    k_use <- k_list[topo_name]
    pca_data <- pca$x
    dist_mat <- dist(pca_data)
    hc <- hclust(dist_mat, method = "ward.D2")
    cluster_labels <- cutree(hc, k = k)
    saveRDS(cluster_labels, clust_fname)
  } else {
    cluster_labels <- readRDS(clust_fname)
  }

  ## PLOTS

  # PCA

  # Heatmap
}


##### IDENTIFY BISTABLE MODELS #####
# def: bistable models are those with multiple steady states across different clusters


##### OU VS WHITE NOISE - ENSEMBLE #####
## Simulate ensembles with different noise type & level, compare summary statistics


##### OU VS WHITE NOISE - SINGLE MODEL #####
## Select a bistable model & simulate long stochastic trajectories, compare transition rate vs noise type & level


##### CLAMPING-DRIVEN TRANSITIONS #####
## Subset bistable models & simulate transitions between clusters


##### EXTRINSIC SIGNAL-DRIVEN TRANSITIONS #####
## Subset bistable models & simulate transitions between clusters












