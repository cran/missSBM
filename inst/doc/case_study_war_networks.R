## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE, 
  warning = FALSE,
  fig.width = 5, fig.height = 4)
set.seed(777) # set seed for reproducibility

## ----package requirement, message=FALSE---------------------------------------
library(igraph)
library(ggplot2)
library(corrplot)
library(magrittr)
library(missSBM)
library(future)

## ----future-plan--------------------------------------------------------------
future::plan("multisession", workers = 2)

## ----load data set------------------------------------------------------------
data("war")

## ----war network plot, fig.width=7,  fig.height=7-----------------------------
par(mar = c(0,0,0,0))
plot(war$belligerent,
     vertex.shape="none", vertex.label=V(war$belligerent)$name,
     vertex.label.color = "steel blue", vertex.label.font=1.5,
     vertex.label.cex=.6, edge.color="gray70", edge.width = 2)

## ----belligerent network------------------------------------------------------
belligerent_adjacency <- as_adj(war$belligerent, sparse = FALSE)
belligerent_power     <- war$belligerent$power
belligerent_trade     <- war$belligerent$trade

## ----sampling node------------------------------------------------------------
partlyObservedNet_war <- missSBM::observeNetwork(belligerent_adjacency, sampling = "node", parameters = .8)
corrplot(partlyObservedNet_war, 
  is.corr      = FALSE,
  tl.pos       = "n",
  method       = "color", 
  cl.pos       = "n",
  na.label.col = "grey",
  mar          = c(0,0,1,0)
  )

## ----inference node, results='hide'-------------------------------------------
vBlocks <- 1:5
collection_sbm <- estimateMissSBM(partlyObservedNet_war, vBlocks, sampling = "node")

## ----inference full, results='hide'-------------------------------------------
collection_sbm_full <- 
  estimateMissSBM(belligerent_adjacency, vBlocks, sampling = "node", control = list(iterates = 2))

## ----plot comparison full-----------------------------------------------------
rbind(
  data.frame(ICL = collection_sbm_full$ICL, nbBlocks = vBlocks, type = "full"),
  data.frame(ICL = collection_sbm$ICL, nbBlocks = vBlocks, type = "missing")
) %>% 
  ggplot(aes(x = nbBlocks, y = ICL, group = type, color = type)) + 
  labs(title = "Model selection", x = "#blocks", y = "Integrated Classification Likelihood") +
  geom_line() + theme_bw()

## ----clustering comparison----------------------------------------------------
table(
  collection_sbm$bestModel$fittedSBM$memberships,
  collection_sbm_full$bestModel$fittedSBM$memberships
  )

## ----plot, fig.width=7, fig.height=7------------------------------------------
par(mfrow = c(1,2))
plot(collection_sbm$bestModel, type = "expected")
plot(collection_sbm_full$bestModel, type = "expected")

## ----war network with covariates full, results = 'hide'-----------------------
vBlocks <- 1:3
collection_sbm_power_full <- estimateMissSBM(belligerent_adjacency, vBlocks = vBlocks, sampling = "node", covariates = list(belligerent_power)) 

## ----power_effect-------------------------------------------------------------
collection_sbm_power_full$bestModel$fittedSBM$covarParam

## ----power_sampling-----------------------------------------------------------
nWar <- nrow(belligerent_adjacency)
parameters_sample <- 600
sampleNet_power_miss <- missSBM::observeNetwork(
   belligerent_adjacency,
   sampling = "covar-node",
   parameters = parameters_sample, covariates = list(belligerent_power), intercept = -2
  )
observedNodes <- !is.na(rowSums(sampleNet_power_miss))
boxplot(1/(1 + exp(-cbind(1,belligerent_power) %*% c(-2, parameters_sample))) ~ observedNodes, ylab="mil power",xlab = "observed node")
corrplot(sampleNet_power_miss, 
  is.corr      = FALSE,
  tl.pos       = "n",
  method       = "color", 
  cl.pos       = "n",
  na.label.col = "grey",
  mar          = c(0,0,1,0)
  )

## ----fit power missing,results = 'hide'---------------------------------------
collection_sbm_power_miss <- estimateMissSBM(sampleNet_power_miss, vBlocks = vBlocks, sampling = "covar-node", covariates = list(belligerent_power))

## ----estimated parameters sample----------------------------------------------
collection_sbm_power_miss$bestModel$fittedSampling$parameters

## ----estimated parameters SBM-------------------------------------------------
collection_sbm_power_miss$bestModel$fittedSBM$covarParam

## ----trade--------------------------------------------------------------------
trade <- belligerent_trade 
trade[is.na(trade)] <- 0
trade <- trade + t(trade)
trade <- log(trade + 1)
diag(trade) = 0

## ----samptrade----------------------------------------------------------------
parameters_sample <- 1
sampleNet_trade_miss <- missSBM::observeNetwork(belligerent_adjacency, sampling = "covar-dyad", parameters = parameters_sample, covariates = list(trade), intercept = -2)
corrplot(sampleNet_trade_miss, 
  is.corr      = FALSE,
  tl.pos       = "n",
  method       = "color", 
  cl.pos       = "n",
  na.label.col = "grey",
  mar          = c(0,0,1,0)
  )

## ----estimate trade-----------------------------------------------------------
collection_sbm_trade_miss <- estimateMissSBM(sampleNet_trade_miss, vBlocks = vBlocks, sampling = "covar-dyad", covariates = list(trade))
collection_sbm_trade_miss$bestModel$fittedSampling$parameters
collection_sbm_trade_miss$bestModel$fittedSBM$covarParam

## ----future-plan-unset--------------------------------------------------------
future::plan("sequential")

