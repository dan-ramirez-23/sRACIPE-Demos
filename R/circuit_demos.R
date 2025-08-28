##### SETUP & GLOBALS #####
rm(list=ls())
#remotes::install_github("lusystemsbio/sRACIPE")
#remotes::install_github("dan-ramirez-23/sRACIPE_clamps@patch-1")
library(sRACIPE)                 # GRN simulation
library(ggplot2)                 # plotting
library(microbenchmark)          # time cost benchmarking
library(dplyr)                   # data management
library(tidyr)                   # data management
library(ComplexHeatmap)          # plotting


# seed for reproducibility & color palette for plots
set.seed(1234)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(3,2,4:8)]


topo_list <- c("toggleSwitch", "CTS", "TS_rep_2") #, "cellcycle")
k_list <- c(2, 4, 3)#, 6) # expected number of clusters for each topo
names(k_list) <- topo_list

nModels <- 1000
nICs <- 100
n_places_rounded <- 4


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
  racipe_fname <- file.path(dataDir,"racipe.Rds")
  if(!file.exists(mb_rs_fname)) {
    mb_rs <- microbenchmark({
      racipe <- sracipeSimulate(topo, numModels = nModels, nIC = nICs)
    }, times = 10)
    saveRDS(mb_rs, mb_rs_fname)
    saveRDS(racipe, racipe_fname)
  } else {
    mb_rs <- readRDS(mb_rs_fname)
    racipe <- readRDS(racipe_fname)
  }


  ## GET UNIQUE STEADY STATES
  ss_unique_fname <- file.path(dataDir,"ss_unique.Rds")
  racipe_summary_fname <- file.path(dataDir,"racipe_summary.Rds")
  if(!file.exists(ss_unique_fname)) {
    states <- round(as.data.frame(t(assay(racipe))), n_places_rounded)
    states$Model <- rep(1:nModels, each=nICs)
    ss_unique <- states %>%
      group_by(Model) %>%
      distinct(across(-Model), .keep_all = TRUE) %>%
      mutate(State = row_number()) %>%
      ungroup()

    racipe_summary_df <- ss_unique %>%
      group_by(Model) %>%
      summarise(NumStates = n(), .groups = "drop")

    saveRDS(ss_unique, ss_unique_fname)
    saveRDS(racipe_summary_df, racipe_summary_fname)
  } else {
    ss_unique <- readRDS(ss_unique_fname)
    racipe_summary_df <- readRDS(racipe_summary_fname)
  }

  ## NORMALIZATION
  ## TODO: save a RACIPE object with unique steady states as ICs for later sims?
  genes <- rownames(racipe)
  expr_mat <- ss_unique[,genes]
  tmpMeans <- rowMeans(t(log2(1+expr_mat)))
  tmpSds <- apply(t(log2(1+expr_mat)),1,sd)
  expr_mat_norm <- log2(1+expr_mat)
  expr_mat_norm[,genes] <- sweep(expr_mat_norm[,genes], 2, tmpMeans, FUN = "-") # scale
  expr_mat_norm[,genes] <- sweep(expr_mat_norm[,genes], 2, tmpSds, FUN = "/") # scale

  ## PCA
  pca_fname <- file.path(dataDir,"pca.Rds")
  if(!file.exists(pca_fname)) {

    pca <- prcomp(expr_mat_norm[,genes])

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
    cluster_labels <- cutree(hc, k = k_use)
    saveRDS(cluster_labels, clust_fname)
  } else {
    cluster_labels <- readRDS(clust_fname)
  }

  ## PLOTS
  # PCA
  pc1_weight <- round(100*summary(pca)$importance[2,1],2)
  pc2_weight <- round(100*summary(pca)$importance[2,2],2)
  plot_xlab <- paste("PC1 (",pc1_weight,"%)",sep="")
  plot_ylab <- paste("PC2 (",pc2_weight,"%)",sep="")
  
  image <- ggplot(pca$x, aes(x=PC1, y=PC2, color=as.factor(cluster_labels))) +
    geom_point(size=2) +
    scale_color_manual(values=cbPalette) +
    xlab(plot_xlab) +
    ylab(plot_ylab) +
    labs(color="Cluster") +
    theme(axis.line = element_line(linewidth = 1, color = "black"), 
          axis.ticks = element_line(linewidth = 1, color="black"),
          panel.background = element_rect("white"),
          plot.background = element_rect("white"),
          axis.title = element_text(size=28),
          axis.text = element_text(size=24),
          legend.title = element_text(size=18),
          legend.text = element_text(size=14))
  image

  # Heatmap
  ha_df <- data.frame(Cluster=cluster_labels)
  
  # Create an annotation object for the columns
  column_annotation <- HeatmapAnnotation(df = ha_df, 
                                         col=list(Cluster=c("1"=unname(cbPalette[1]),"2"=unname(cbPalette[2]),
                                                            "3"=unname(cbPalette[3]),"4"=unname(cbPalette[4]))))
  # Create the heatmap with annotation
  image <- Heatmap(as.matrix(t(expr_mat_norm)), 
                   name = "Expression", 
                   top_annotation = column_annotation,
                   row_names_gp=gpar(fontsize=12),
                   clustering_method_columns = "ward.D2",
                   show_column_names = F)
  image
  

}


##### DEMO NEW FEATURES ON toggleSwitch CIRCUIT #####
topo_name <- "toggleSwitch"

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


mb_rs_fname <- file.path(dataDir,"microbenchmark_simulation_results.Rds")
racipe_fname <- file.path(dataDir,"racipe.Rds")
ss_unique_fname <- file.path(dataDir,"ss_unique.Rds")
mb_rs <- readRDS(mb_rs_fname)
racipe <- readRDS(racipe_fname)
ss_unique <- readRDS(ss_unique_fname)

genes <- rownames(racipe)
expr_mat <- ss_unique[,genes]
tmpMeans <- rowMeans(t(log2(1+expr_mat)))
tmpSds <- apply(t(log2(1+expr_mat)),1,sd)
expr_mat_norm <- log2(1+expr_mat)
expr_mat_norm[,genes] <- sweep(expr_mat_norm[,genes], 2, tmpMeans, FUN = "-") # scale
expr_mat_norm[,genes] <- sweep(expr_mat_norm[,genes], 2, tmpSds, FUN = "/") # scale

ss_unique_fname <- file.path(dataDir,"ss_unique.Rds")
racipe_summary_fname <- file.path(dataDir,"racipe_summary.Rds")
ss_unique <- readRDS(ss_unique_fname)
racipe_summary_df <- readRDS(racipe_summary_fname)

pca_fname <- file.path(dataDir,"pca.Rds")
pca <- readRDS(pca_fname)

clust_fname <- file.path(dataDir,"cluster_labels.Rds")
cluster_labels <- readRDS(clust_fname)


##### IDENTIFY BISTABLE MODELS #####
# def: bistable models are those with multiple steady states across different clusters
ss_unique$Cluster <- cluster_labels
ss_unique <- as.data.frame(ss_unique)
racipe_summary_df$Stability <- NA
for(model in racipe_summary_df$Model) {
  numModelClusters <- length(unique(ss_unique[which(ss_unique$Model == model),"Cluster"]))
  if(numModelClusters == 2) {
    racipe_summary_df[which(racipe_summary_df$Model == model), "Stability"] <- "Bistable"
  }
}
bistable_models <- which(racipe_summary_df$Stability == "Bistable" & racipe_summary_df$NumStates == 2)




##### CLAMPING-DRIVEN TRANSITIONS #####
## Subset bistable models & simulate transitions between clusters by clamping gene values
keepIdx <- c()
for(model in bistable_models) {
  # add steady states for cluster 1 and 2
  addIdx <- which(ss_unique$Model == model)
  keepIdx <- c(keepIdx, addIdx)
  
}
clamp_df <- pivot_longer(ss_unique[keepIdx,], cols = all_of(genes),
                         names_to = "Gene", values_to = "Expression")
clamp_df$ModelIndex <- as.numeric(factor(clamp_df$Model))

# set up racipe object with params & ICs from bistable models (paired initialization)
racipe_bistable <- sracipeSimulate(topo, numModels = (2*length(bistable_models)), nIC = 1,
                                   genIC = T, genParams = T, integrate = F)
sracipeIC(racipe_bistable) <- t(ss_unique[keepIdx,genes])
sracipeParams(racipe_bistable) <- sracipeParams(racipe)[rep(bistable_models, each=2),]

# clamp B to drive toward cluster 2 values
clamp_df_filt_A <- clamp_df[which(clamp_df$Cluster == 2 & clamp_df$Gene == "A"),]
clamp_df_filt_B <- clamp_df[which(clamp_df$Cluster == 2 & clamp_df$Gene == "B"),]
clamp_df_use <- data.frame(#A_NULL=rep(clamp_df_filt_A$Expression, each=2), 
                           B=rep(clamp_df_filt_B$Expression, each=2))

# simulate with clamping
## TODO: continuing bug with gene clamping, needs a PR
undebug(sracipeSimulate)
racipe_clamped <- sracipeSimulate(racipe_bistable, integrate=T, genIC = F, genParams = F,
                                  geneClamping = clamp_df_use, numModels = 326, nIC = 1, stepper = "EM",
                                  initialNoise = 0.2, scaledNoise = T, nNoise = 1, simDet = F)
sracipeConvergeDist(racipe_clamped)

assay(racipe_clamped)

racipe_clamp_relax <- racipe_clamped
sracipeIC(racipe_clamp_relax) <- assay(racipe_clamped)
racipe_clamp_relax <- sracipeSimulate(racipe_clamp_relax, integrate=T, genIC=F, genParams = F,
                                      numModels = 326, nIC=1, stepper = "EM", initialNoise = 0, nNoise = 0, simDet = T)


# plot initial & final states from clamping
ics_raw <- as.data.frame(t(sracipeIC(racipe_bistable)))
ics_norm <- log2(1+ics_raw)
ics_norm[,genes] <- sweep(ics_norm[,genes], 2, tmpMeans, FUN = "-") # scale
ics_norm[,genes] <- sweep(ics_norm[,genes], 2, tmpSds, FUN = "/") # scale
ics_pca <- as.data.frame(as.matrix(ics_norm[, c(1, 2)]) %*% pca$rotation)
ics_norm$Time <- "Initial"

final_states_raw <- as.data.frame(t(sracipeIC(racipe_clamp_relax)))
final_states_norm <- log2(1+final_states_raw)
final_states_norm[,genes] <- sweep(final_states_norm[,genes], 2, tmpMeans, FUN = "-") # scale
final_states_norm[,genes] <- sweep(final_states_norm[,genes], 2, tmpSds, FUN = "/") # scale
final_states_pca <- as.data.frame(as.matrix(final_states_norm[, c(1, 2)]) %*% pca$rotation)
final_states_norm$Time <- "Final"


states_comb <- rbind(ics_norm, final_states_norm)


image <- ggplot(states_comb[,],aes(x=A, y=B, color=Time)) +
  geom_point(alpha=0.6, size=3) +
  labs(x = "A", y = "B", color="Time") +
  theme_minimal() +
  theme(axis.text = element_text(size=24),
        axis.title = element_text(size=28),)

pdf(file = file.path(plotDir,"clamping_before_vs_after.pdf"), width = 10, height = 10)
print(image)
dev.off()





##### EXTRINSIC SIGNAL-DRIVEN TRANSITIONS #####
## Subset bistable models & simulate transitions between clusters using parameter modifications
# As opposed to clamping B, here we will drive toward cluster 2 by increasing G_B 20x
param_signaling_df <- data.frame(Time=c(0, 25, 50, 200), 
                                 G_B=c(1,200,1,1))
#debug(sracipeSimulate)
racipe_signaling <- sracipeSimulate(racipe_bistable, integrate=T, genIC = F, genParams = F,
                                    paramSignalVals = param_signaling_df, simulationTime = 200,
                                    initialNoise = 0.2, nNoise = 1, simDet = T, stepper = "EM", scaledNoise = T)


racipe_signaling_relax <- racipe_signaling
sracipeIC(racipe_signaling_relax) <- assay(racipe_signaling,"deterministic")
racipe_signaling_relax <- sracipeSimulate(racipe_signaling_relax, integrate=T, 
                                          genIC = F, genParams = F,
                                          simulationTime = 100,
                                          initialNoise = 0, nNoise = 0, simDet = T, stepper="EM")


# sanity checks
sracipeIC(racipe_bistable)[1:2,1:12] # ICs
assay(racipe_signaling, "deterministic")[1:2,1:12] # det signaling
assay(racipe_signaling, "0.2")[1:2,1:12] # stoch signaling
assay(racipe_signaling_relax)[1:2,1:12] # after relaxation

## tally how many models transitioned between clusters
initial_states <- t(sracipeIC(racipe_signaling))
final_states <- t(assay(racipe_signaling, "0.2"))
diffs <- final_states - t(sracipeIC(racipe_signaling))
colSums(diffs)

relaxed_states <- t(assay(racipe_signaling_relax))
diffs <- relaxed_states - t(sracipeIC(racipe_signaling))
colSums(diffs)

param_signaling_transition_summary <- data.frame(Model=rep(bistable_models, each=2),
                                                 Init_Clust=cluster_labels[keepIdx],
                                                 Final_Clust=NA
                                                 )
for(i in 1:nrow(relaxed_states)) {
  model <- bistable_models[round((i+0.1) / 2)]
  model_states <- ss_unique[which(ss_unique$Model == model),]
  
  final_state <- relaxed_states[c(i, i),]
  diffs = final_state - model_states[,genes]
  diffs$Cluster <- model_states$Cluster
  diff_sums <- rowSums(diffs[,genes])
  names(diff_sums) <- model_states$Cluster
  
  if(abs(diff_sums[which(names(diff_sums) == 1)]) < 1e-4) {
    param_signaling_transition_summary[i,"Final_Clust"] <- 1
  } else if(diff_sums[which(names(diff_sums) == 2)] < 1e-4) {
    param_signaling_transition_summary[i,"Final_Clust"] <- 2
  } else {
    param_signaling_transition_summary[i,"Final_Clust"] <- "New"
  }
    
}
param_signaling_transition_summary$Transition <- "None"
param_signaling_transition_summary[which(param_signaling_transition_summary$Init_Clust == 1 &
                                        param_signaling_transition_summary$Final_Clust == 2   ),"Transition"] <- "1-->2"
param_signaling_transition_summary[which(param_signaling_transition_summary$Init_Clust == 2 &
                                           param_signaling_transition_summary$Final_Clust == 1   ),"Transition"] <- "2-->1"
table(param_signaling_transition_summary$Transition)





# plot initial & final states from clamping
ics_raw <- as.data.frame(t(sracipeIC(racipe_bistable)))
ics_norm <- log2(1+ics_raw)
ics_norm[,genes] <- sweep(ics_norm[,genes], 2, tmpMeans, FUN = "-") # scale
ics_norm[,genes] <- sweep(ics_norm[,genes], 2, tmpSds, FUN = "/") # scale
ics_pca <- as.data.frame(as.matrix(ics_norm[, c(1, 2)]) %*% pca$rotation)
ics_norm$Time <- "Initial"

final_states_raw <- as.data.frame(t(assay(racipe_signaling)))
final_states_norm <- log2(1+final_states_raw)
final_states_norm[,genes] <- sweep(final_states_norm[,genes], 2, tmpMeans, FUN = "-") # scale
final_states_norm[,genes] <- sweep(final_states_norm[,genes], 2, tmpSds, FUN = "/") # scale
final_states_pca <- as.data.frame(as.matrix(final_states_norm[, c(1, 2)]) %*% pca$rotation)
final_states_norm$Time <- "Final"


states_comb <- rbind(ics_norm, final_states_norm)


image <- ggplot(states_comb,aes(x=A, y=B, color=Time, shape=Time)) +
  geom_point(alpha=0.8, size=3) +
  labs(x = "A", y = "B", color="Time") +
  theme_minimal() +
  theme(axis.text = element_text(size=24),
        axis.title = element_text(size=28),)
image

pdf(file = file.path(plotDir,"signaling_before_vs_after.pdf"), width = 10, height = 10)
print(image)
dev.off()






##### OU VS WHITE NOISE - ENSEMBLE #####
## Simulate ensembles with different noise type & level, compare summary statistics

# long-term variance of WN: d^2 * h  
  # 0.2^2 * 0.2 = 0.04 * 0.2 = 0.008
# of OU: d^2 / (2*a), where a is 1/tau
  # 0.04^2 / (2*0.1) = 0.0016 / 0.2 = 0.008


0.8^2 * 0.2
0.16^2 / (2*0.1)

## Loop over param sets and plot side-by-side ensembles for diff parameter sets

param_sets <- data.frame(WN_level=c(0.05, 0.1, 0.15, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.08, 0.08, 0.08),
                         OU_level=c(0.05, 0.1, 0.15, 0.4, 0.4, 0.4, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08),
                         OU_tau=c(10, 10, 10, 1, 10, 100, 1, 10, 100, 1, 10, 100))

for(pset in 1:nrow(param_sets)) {
  wn_noise <- param_sets[pset,"WN_level"]
  ou_noise <- param_sets[pset,"OU_level"]
  ou_tau   <- param_sets[pset,"OU_tau"]
  
  
  racipe_wn_ensemble <- sracipeSimulate(racipe_bistable, integrate=T, genIC = F, genParams = F,
                                        simulationTime = 100, initialNoise = 0.2, nNoise = 1,
                                        simDet = F, integrateStepSize = 0.2, stepper = "EM")
  # racipe_wn_ensemble_relax <- racipe_wn_ensemble
  # sracipeIC(racipe_wn_ensemble_relax) <- assay(racipe_wn_ensemble)
  # racipe_wn_ensemble_relax <- sracipeSimulate(racipe_wn_ensemble_relax, integrate=T, genIC = F, genParams = F,
  #                                             initialNoise = 0, nNoise = 0, simDet = T, stepper = "EM")
  
  
  racipe_ou_ensemble <- sracipeSimulate(racipe_bistable, integrate=T, genIC = F, genParams = F,
                                        simulationTime = 100, initialNoise = 0.04, nNoise = 1,
                                        simDet = F, ouNoise_t = 10, integrateStepSize = 0.2, stepper = "EM")
  # racipe_ou_ensemble_relax <- racipe_ou_ensemble
  # sracipeIC(racipe_ou_ensemble_relax) <- assay(racipe_ou_ensemble)
  # racipe_ou_ensemble_relax <- sracipeSimulate(racipe_ou_ensemble_relax, integrate=T, genIC = F, genParams = F,
  #                                             initialNoise = 0, nNoise = 0, simDet = T, stepper = "EM")
  
  
  
  ou_ensemble_raw <- as.data.frame(t(assay(racipe_ou_ensemble)))
  ou_ensemble_norm <- log2(1+ou_ensemble_raw)
  ou_ensemble_norm[,genes] <- sweep(ou_ensemble_norm[,genes], 2, tmpMeans, FUN = "-") # scale
  ou_ensemble_norm[,genes] <- sweep(ou_ensemble_norm[,genes], 2, tmpSds, FUN = "/") # scale
  ou_ensemble_pca <- as.data.frame(as.matrix(ou_ensemble_norm[, c(1, 2)]) %*% pca$rotation)
  ou_ensemble_pca$Type <- paste0("OU (tau=",ou_tau,")")
  
  wn_ensemble_raw <- as.data.frame(t(assay(racipe_wn_ensemble)))
  wn_ensemble_norm <- log2(1+wn_ensemble_raw)
  wn_ensemble_norm[,genes] <- sweep(wn_ensemble_norm[,genes], 2, tmpMeans, FUN = "-") # scale
  wn_ensemble_norm[,genes] <- sweep(wn_ensemble_norm[,genes], 2, tmpSds, FUN = "/") # scale
  wn_ensemble_pca <- as.data.frame(as.matrix(wn_ensemble_norm[, c(1, 2)]) %*% pca$rotation)
  wn_ensemble_pca$Type <- "White"
  
  traj_ensemble_combined_norm <- rbind(ou_ensemble_norm, wn_ensemble_norm)
  traj_ensemble_combined_norm$Type <- c(rep(paste0("OU (tau=",ou_tau,")"), nrow(ou_ensemble_raw)), rep("White", nrow(wn_ensemble_raw)))
  traj_ensemble_combined <- rbind(ou_ensemble_pca, wn_ensemble_pca)
  
  
  
  image <- ggplot(traj_ensemble_combined_norm, aes(x=A, y=B, color = Type)) +
    geom_point(alpha = 0.6, size=3) +
    labs(x = "A", y = "B", color="Type") +
    theme_minimal() +
    theme(axis.text = element_text(size=24),
          axis.title = element_text(size=28),)
  
  pdf(file = file.path(plotDir,paste0("WN=",wn_noise,"_vs_OU=",ou_noise,"_tau=",ou_tau,"_ensemble_AB.pdf")), width = 10, height = 10)
  print(image)
  dev.off()
  
  # plot on PCA
  image <- ggplot(traj_ensemble_combined,aes(x=PC1, y=PC2, color=Type)) +
    geom_point(alpha=0.6, size=3) +
    labs(x = "PC1", y = "PC2", color="Type") +
    theme_minimal() +
    theme(axis.text = element_text(size=24),
          axis.title = element_text(size=28),)
  
  pdf(file = file.path(plotDir,paste0("WN=",wn_noise,"_vs_OU=",ou_noise,"_tau=",ou_tau,"_ensemble_PCA.pdf")), width = 10, height = 10)
  print(image)
  dev.off()
  
  
}





##### OU VS WHITE NOISE - SINGLE MODEL #####
## Select a bistable model & simulate long stochastic trajectories, compare transition rate vs noise type & level
long_sim_one_model_time <- 500
model_use <- bistable_models[1]
wn_noise <- 0.2
ou_noise <- 0.04
ou_tau <- 10#1e-5


ic_use <- t(ss_unique[which(ss_unique$Model == model_use & 
                              ss_unique$Cluster == 1),genes])
params_use_placeholder <- sracipeParams(racipe)[model_use,]
params_use <- params_use_placeholder
params_use <- c(70, 70, 0.2, 0.2, 5, 5, 3, 3, 80, 80)



racipe_wn_test <- sracipeSimulate(topo, integrate = F, numModels = 1, nIC = 100,
                                  genIC = T, genParams = T, integrateStepSize = 0.2,
                                  nNoise = 0, initialNoise = 0)
sracipeParams(racipe_wn_test) <- t(params_use)
racipe_wn_test <- sracipeSimulate(racipe_wn_test, integrate=T, numModels = 100, nIC = 1, genParams = F)
model_states <- distinct(as.data.frame(t(assay(racipe_wn_test))))
# steady states at (211, 4) and (4, 211)

racipe_wn_one_model <- sracipeSimulate(topo, integrate = F, numModels = 1, nIC = 1,
                                       timeSeries = T, simulationTime = long_sim_one_model_time, 
                                       genIC = T, genParams = T, 
                                       printStart = 1, printInterval = 1, integrateStepSize = 0.2,
                                       nNoise = 1, initialNoise = wn_noise,
                                       simDet = F)
sracipeIC(racipe_wn_one_model) <- ic_use
sracipeParams(racipe_wn_one_model) <- params_use
racipe_wn_one_model <- sracipeSimulate(racipe_wn_one_model, integrate = T, 
                                       genIC = F, genParams = F, simulationTime = long_sim_one_model_time,
                                       timeSeries = T, nNoise = 1, initialNoise = wn_noise,
                                       simDet = F, printStart = 1, printInterval = 1, 
                                       integrateStepSize = 0.2)
trajectory_wn <- as.data.frame(t(racipe_wn_one_model@metadata$timeSeries))



racipe_ou_one_model <- sracipeSimulate(topo, integrate = F, numModels = 1, nIC = 1,
                                       timeSeries = T, simulationTime = long_sim_one_model_time, 
                                       genIC = T, genParams = T, 
                                       printStart = 1, printInterval = 1, integrateStepSize = 0.2,
                                       nNoise = 1, initialNoise = ou_noise, ouNoise_t = ou_tau,
                                       simDet = F, stepper = "EM_OU")
sracipeIC(racipe_ou_one_model) <- ic_use
sracipeParams(racipe_ou_one_model) <- params_use
racipe_ou_one_model <- sracipeSimulate(racipe_ou_one_model, integrate = T, 
                                       genIC = F, genParams = F, simulationTime = long_sim_one_model_time,
                                       timeSeries = T, nNoise = 1, initialNoise = ou_noise,
                                       simDet = F, printStart = 1, printInterval = 1, 
                                       integrateStepSize = 0.2, ouNoise_t = ou_tau, stepper = "EM_OU")
trajectory_ou <- as.data.frame(t(racipe_ou_one_model@metadata$timeSeries))
for(gene in genes) {
  trajectory_wn[,gene] <- as.numeric(trajectory_wn[,gene])
  trajectory_ou[,gene] <- as.numeric(trajectory_ou[,gene])
}



# ggplot() +
#   geom_path(data=trajectory_wn, aes(x=A, y=B), color="red") +
#   geom_path(data=trajectory_ou, aes(x=A, y=B), color="blue")

ggplot() +
  geom_line(aes(x=1:nrow(trajectory_ou), y=trajectory_ou$A)) 



##### OU VS WHITE NOISE - SINGLE MODEL DETAILED #####

# ---- simulation hyperparameters -----------------------------------------------------
genes          <- c("A","B")
sim_time_long  <- 5000
h              <- 0.2               # integrateStepSize you already use
print_start    <- 1
print_interval <- 1

wn_noise <- 0.40                    # pick by eye so you get occasional switches
ou_noise <- 0.08                    # pick by eye per tau so in-well wiggle is similar
ou_tau   <- 20                      # try a few: c(0.5, 2, 10, 50)



model_states$Cluster <- NA
model_states[which(model_states$A > model_states$B),"Cluster"] <- 1
model_states[which(model_states$B > model_states$A),"Cluster"] <- 2
ic_use <- model_states[model_states$Cluster == 1, genes]
colnames(ic_use) <- genes

params_use <- c(70, 70, 0.2, 0.2, 5, 5, 3, 3, 80, 80)

model_states
params_use
ic_use



model_use <- bistable_models[1]
P1 <- as.vector(unlist(model_states[2,c("A","B")]))#as.numeric(ss_unique[ss_unique$Model == model_use & ss_unique$Cluster == 1, genes][1, ])
P2 <- as.vector(unlist(model_states[1,c("A","B")]))#as.numeric(ss_unique[ss_unique$Model == model_use & ss_unique$Cluster == 2, genes][1, ])
names(P1) <- genes
names(P2) <- genes

ic_use     <- t(P1)                 # start in P1 basin; you can also start at P2
#params_use <- sracipeParams(racipe)[model_use, ]

# ---- helpers ----------------------------------------------------------------
q_of <- function(df) df$A - df$B

sum_to_q <- function(sum_x, q) {
  x1 <- 0.5 * (sum_x + q); x2 <- 0.5 * (sum_x - q)
  c(A = max(x1, 0), B = max(x2, 0))
}

# naive dividing value + hysteresis (no Jacobians)
make_thresholds <- function(P1, P2, frac = 0.10, max_frac = 0.30) {
  q1 <- P1["A"] - P1["B"]; q2 <- P2["A"] - P2["B"]
  q_mid <- 0.5 * (q1 + q2)
  dQ <- abs(q2 - q1)
  delta <- min(frac * dQ, max_frac * dQ) # 10% of separation, capped
  list(q_lo = q_mid - delta, q_hi = q_mid + delta,
       q_mid = q_mid, dQ = dQ, q1 = q1, q2 = q2)
}

# count switches with tiny debounce & hysteresis
count_switches <- function(q, q_lo, q_hi, min_stay_steps = 3) {
  # states: -1 (P1 core), +1 (P2 core), 0 (middle)
  state <- ifelse(q <= q_lo, -1L, ifelse(q >= q_hi, +1L, 0L))
  n <- length(state)
  
  # compress consecutive runs
  runs <- rle(state)
  # keep only core runs with enough stay
  core_idx <- which(runs$values != 0 & runs$lengths >= min_stay_steps)
  
  if (length(core_idx) <= 1) {
    return(list(n_switch = 0L, switch_indices = integer(0)))
  }
  
  core_states <- runs$values[core_idx]
  # switches are sign changes between consecutive qualifying core runs
  n_switch <- sum(core_states[-1] != core_states[-length(core_states)])
  # approximate indices where each qualifying run starts
  run_starts <- cumsum(c(1L, runs$lengths[-length(runs$lengths)]))
  switch_pos <- run_starts[core_idx[-1]]
  
  list(n_switch = as.integer(n_switch), switch_indices = switch_pos)
}

# build correct time vector from printStart/printInterval
make_time <- function(n_rows, print_start, print_interval) {
  seq(from = print_start, by = print_interval, length.out = n_rows)
}

# simple “runner” wrappers that return A,B dataframe for a given sim_time
make_white_runner <- function(initialNoise, topo, params_use, ic_use, h,
                              print_start, print_interval) {
  function(sim_time = 200) {
    r <- sracipeSimulate(topo, integrate = FALSE, numModels = 1, nIC = 1,
                         timeSeries = TRUE, simulationTime = sim_time,
                         genIC = TRUE, genParams = TRUE,
                         printStart = print_start, printInterval = print_interval,
                         integrateStepSize = h, nNoise = 1,
                         initialNoise = initialNoise, simDet = FALSE, scaledNoise = T)
    sracipeIC(r) <- ic_use
    sracipeParams(r) <- t(params_use)
    r <- sracipeSimulate(r, integrate = TRUE, genIC = FALSE, genParams = FALSE,
                         simulationTime = sim_time, timeSeries = TRUE,
                         nNoise = 1, initialNoise = initialNoise, simDet = FALSE,
                         printStart = print_start, printInterval = print_interval,
                         integrateStepSize = h, scaledNoise = T)
    df <- as.data.frame(t(r@metadata$timeSeries))[, c("A","B")]
    for (j in seq_len(ncol(df))) df[, j] <- as.numeric(df[, j])
    df
  }
}

make_ou_runner <- function(initialNoise, tau, topo, params_use, ic_use, h,
                           print_start, print_interval) {
  function(sim_time = 200) {
    r <- sracipeSimulate(topo, integrate = FALSE, numModels = 1, nIC = 1,
                         timeSeries = TRUE, simulationTime = sim_time,
                         genIC = TRUE, genParams = TRUE,
                         printStart = print_start, printInterval = print_interval,
                         integrateStepSize = h, nNoise = 1,
                         initialNoise = initialNoise, ouNoise_t = tau,
                         simDet = FALSE, stepper = "EM_OU", scaledNoise = T)
    sracipeIC(r) <- ic_use
    sracipeParams(r) <- t(params_use)
    r <- sracipeSimulate(r, integrate = TRUE, genIC = FALSE, genParams = FALSE,
                         simulationTime = sim_time, timeSeries = TRUE,
                         nNoise = 1, initialNoise = initialNoise, simDet = FALSE,
                         printStart = print_start, printInterval = print_interval,
                         integrateStepSize = h, ouNoise_t = tau, stepper = "EM_OU", scaledNoise = T)
    df <- as.data.frame(t(r@metadata$timeSeries))[, c("A","B")]
    for (j in seq_len(ncol(df))) df[, j] <- as.numeric(df[, j])
    df
  }
}


# ---- run sims ---------------------------------------------------------------
thr <- make_thresholds(P1, P2, frac = 0.10)

#undebug(make_white_runner)
runner_white <- make_white_runner(wn_noise, topo, params_use, t(ic_use), h,
                                  print_start, print_interval)
runner_ou    <- make_ou_runner(ou_noise, ou_tau, topo, params_use, t(ic_use), h,
                               print_start, print_interval)

trajW <- runner_white(sim_time_long)
trajO <- runner_ou(sim_time_long)

qW <- q_of(trajW); qO <- q_of(trajO)
tW <- make_time(nrow(trajW), print_start, print_interval)
tO <- make_time(nrow(trajO), print_start, print_interval)

# ---- count crude switches ---------------------------------------------------
swW <- count_switches(qW, thr$q_lo, thr$q_hi, min_stay_steps = 3)
swO <- count_switches(qO, thr$q_lo, thr$q_hi, min_stay_steps = 3)

rateW <- swW$n_switch / (tail(tW, 1) - head(tW, 1))
rateO <- swO$n_switch / (tail(tO, 1) - head(tO, 1))

cat(sprintf("White: %d switches (%.3f per time unit)\n", swW$n_switch, rateW))
cat(sprintf("OU(τ=%.2f): %d switches (%.3f per time unit)\n", ou_tau, swO$n_switch, rateO))

# ---- quick plots ------------------------------------------------------------
dfW <- data.frame(t = tW, q = qW, type = "White")
dfO <- data.frame(t = tO, q = qO, type = paste0("OU (tau=", ou_tau, ")"))
both <- rbind(dfW, dfO)

ggplot(both, aes(t, q, color = type)) +
  geom_line(linewidth = 0.4, alpha = 0.9) +
  geom_hline(yintercept = c(thr$q_lo, thr$q_hi), linetype = "dashed") +
  labs(x = "Time", y = "q = A - B", title = "White vs OU noise: q(t)") +
  theme_minimal()


##### OU VS WHITE NOISE - SINGLE MODEL MULTIPLE SIMS #####

genes          <- c("A","B")
sim_time_long  <- 5000
h              <- 0.2             
print_start    <- 1
print_interval <- 1

model_states <- distinct(as.data.frame(t(assay(racipe_wn_test))))
model_states$Cluster <- NA
model_states[which(model_states$A > model_states$B),"Cluster"] <- 1
model_states[which(model_states$B > model_states$A),"Cluster"] <- 2

params_use <- c(70, 70, 0.2, 0.2, 5, 5, 3, 3, 80, 80)

P1 <- as.vector(unlist(model_states[2,c("A","B")]))
P2 <- as.vector(unlist(model_states[1,c("A","B")]))
names(P1) <- genes
names(P2) <- genes

ic_use     <- t(P1)               
colnames(ic_use) <- genes

param_sets <- data.frame(WN_level=c(0.05, 0.1, 0.15),#, 0.08, 0.08, 0.08, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
                         OU_level=c(0.05, 0.1, 0.15),#, 0.08, 0.08, 0.08, 0.4, 0.4, 0.4, 0.08, 0.08, 0.08),
                         OU_tau=c(10, 10, 10 ))#, 1, 10, 100, 1, 10, 100, 1, 10, 100))

for(pset in 1:nrow(param_sets)) {
  wn_noise <- param_sets[pset,"WN_level"]
  ou_noise <- param_sets[pset,"OU_level"]
  ou_tau   <- param_sets[pset,"OU_tau"]
  
  runner_white <- make_white_runner(wn_noise, topo, params_use, t(ic_use), h,
                                    print_start, print_interval)
  runner_ou    <- make_ou_runner(ou_noise, ou_tau, topo, params_use, t(ic_use), h,
                                 print_start, print_interval)
  
  trajW <- runner_white(sim_time_long)
  trajO <- runner_ou(sim_time_long)
  
  qW <- q_of(trajW); qO <- q_of(trajO)
  tW <- make_time(nrow(trajW), print_start, print_interval)
  tO <- make_time(nrow(trajO), print_start, print_interval)
  
  swW <- count_switches(qW, thr$q_lo, thr$q_hi, min_stay_steps = 3)
  swO <- count_switches(qO, thr$q_lo, thr$q_hi, min_stay_steps = 3)
  
  rateW <- swW$n_switch / (tail(tW, 1) - head(tW, 1))
  rateO <- swO$n_switch / (tail(tO, 1) - head(tO, 1))
  
  #cat(sprintf("White: %d switches (%.3f per time unit)\n", swW$n_switch, rateW))
  #cat(sprintf("OU(τ=%.2f): %d switches (%.3f per time unit)\n", ou_tau, swO$n_switch, rateO))
  
  dfW <- data.frame(t = tW, q = qW, type = "White")
  dfO <- data.frame(t = tO, q = qO, type = paste0("OU (tau=", ou_tau, ")"))
  both <- rbind(dfW, dfO)
  
  image <- ggplot(both, aes(t, q, color = type)) +
    geom_line(linewidth = 0.4, alpha = 0.9) +
    geom_hline(yintercept = c(thr$q_lo, thr$q_hi), linetype = "dashed") +
    labs(x = "Time", y = "q = A - B", color="Type") +
    theme_minimal() +
    theme(axis.text = element_text(size=24),
          axis.title = element_text(size=28),)
  
  pdf(file = file.path(plotDir,paste0("WN=",wn_noise,"_vs_OU=",ou_noise,"_tau=",ou_tau,"_qplot.pdf")), width = 10, height = 10)
  print(image)
  dev.off()
  
  # plot on PCA
  wn_traj_norm <- log2(1+trajW[which(trajW$A != 0 & trajW$B != 0),])
  wn_traj_norm[,genes] <- sweep(wn_traj_norm[,genes], 2, tmpMeans, FUN = "-") # scale
  wn_traj_norm[,genes] <- sweep(wn_traj_norm[,genes], 2, tmpSds, FUN = "/") # scale
  wn_traj_pca <- as.data.frame(as.matrix(wn_traj_norm[, c(1, 2)]) %*% pca$rotation)
  wn_traj_pca$Type <- "White"
  
  ou_traj_norm <- log2(1+trajO[which(trajO$A != 0 & trajO$B != 0),])
  ou_traj_norm[,genes] <- sweep(ou_traj_norm[,genes], 2, tmpMeans, FUN = "-") # scale
  ou_traj_norm[,genes] <- sweep(ou_traj_norm[,genes], 2, tmpSds, FUN = "/") # scale
  ou_traj_pca <- as.data.frame(as.matrix(ou_traj_norm[, c(1, 2)]) %*% pca$rotation)
  ou_traj_pca$Type <- paste0("OU (tau=",ou_tau,")")
  
  traj_pca_combined <- rbind(wn_traj_pca, ou_traj_pca)
  
  
  model_states_norm <- log2(1+model_states[,genes])
  model_states_norm[,genes] <- sweep(model_states_norm[,genes], 2, tmpMeans, FUN = "-") # scale
  model_states_norm[,genes] <- sweep(model_states_norm[,genes], 2, tmpSds, FUN = "/") # scale
  model_states_pca <- as.matrix(model_states_norm[, c(1, 2)]) %*% pca$rotation
  
  
  image <- ggplot(traj_pca_combined,aes(x=PC1, y=PC2, color=Type)) +
    geom_point(alpha=0.5) +
    geom_point(data=model_states_pca, aes(x=PC1, y=PC2), color="red", size=5) +
    labs(x = "PC1", y = "PC2", color="Type") +
    theme_minimal() +
    theme(axis.text = element_text(size=24),
          axis.title = element_text(size=28),)
  
  pdf(file = file.path(plotDir,paste0("WN=",wn_noise,"_vs_OU=",ou_noise,"_tau=",ou_tau,"_PCA.pdf")), width = 10, height = 10)
  print(image)
  dev.off()
  
  
  
}




