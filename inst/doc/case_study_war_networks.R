## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE, 
  warning = FALSE,
  fig.width = 5, fig.height = 4)
set.seed(9877) # set seed for reproducibility

## ----package requirement, message=FALSE----------------------------------
library(igraph)
library(ggplot2)
library(magrittr)
library(missSBM)

## ----load data set-------------------------------------------------------
data("war")

## ----war network plot, fig.width=7,  fig.height=7------------------------
par(mar = c(0,0,0,0))
plot(war$beligerent, 
     vertex.shape="none", vertex.label=V(war$beligerent)$name,
     vertex.label.color = "steel blue", vertex.label.font=1.5,
     vertex.label.cex=.6, edge.color="gray70", edge.width = 2)

## ----beligenrent network-------------------------------------------------
beligerent_adjacency <- as_adj(war$beligerent, sparse = FALSE)
beligerent_power     <- war$beligerent$power
beligerent_trade     <- war$beligerent$trade

## ----sampling node-------------------------------------------------------
sampledNet_war <- missSBM::sample(beligerent_adjacency, sampling = "node", parameters = .8)
plot(sampledNet_war)

## ----inference node, results='hide'--------------------------------------
vBlocks <- 1:5
collection_sbm <- missSBM::estimate(sampledNet_war, vBlocks = vBlocks, sampling = "node")
res_unsmoothed <- data.frame(
  ICL     = collection_sbm$ICL,
  nBlocks = vBlocks, 
  type    = "raw"
)

## ----smoothed node, results='hide'---------------------------------------
smooth(collection_sbm, "both")
res_smoothed <- data.frame(
  ICL     = collection_sbm$ICL,
  nBlocks = vBlocks, 
  type    = "smoothed"
)

## ----smoothing effect plot, fig.width = 7, fig.height = 5----------------
rbind(res_unsmoothed, res_smoothed) %>% 
  ggplot(aes(x = nBlocks, y = ICL, group = type, color = type)) + 
    geom_line() + theme_bw()

## ----inference full, results='hide'--------------------------------------
collection_sbm_full <- 
  missSBM::estimate(
    sampledNet  = prepare_data(beligerent_adjacency), 
    vBlocks     = vBlocks, 
    sampling    = "node"
  )
smooth(collection_sbm_full, "forward", control = list(iterates = 3))

## ----plot comparison full------------------------------------------------
res_missing <- res_smoothed
res_missing$type <- "missing"
res_full <- data.frame(
  ICL     = collection_sbm_full$ICL,
  nBlocks = vBlocks, 
  type    = "full"
)
rbind(res_missing, res_full) %>% 
  ggplot(aes(x = nBlocks, y = ICL, group = type, color = type)) + 
    geom_line() + theme_bw()


## ----clustering comparison-----------------------------------------------
table(
  collection_sbm$bestModel$fittedSBM$memberships,
  collection_sbm_full$bestModel$fittedSBM$memberships
  )

## ----plot, fig.width=7, fig.height=7-------------------------------------
par(mfrow = c(2,2))
plot(collection_sbm$bestModel$fittedSBM, type = "network")
plot(collection_sbm$bestModel$fittedSBM, type = "connectivity")
plot(collection_sbm_full$bestModel$fittedSBM, type = "network")
plot(collection_sbm_full$bestModel$fittedSBM, type = "connectivity")

## ----war network with power----------------------------------------------
sampleNet_cov <- prepare_data(beligerent_adjacency, list(beligerent_power)) 

## ----war network with covariates full, results = 'hide'------------------
vBlocks <- 1:4
collection_sbm_power_full <- estimate(sampleNet_cov, vBlocks = vBlocks, sampling = "node", useCovariates = TRUE)

## ----power_effect--------------------------------------------------------
collection_sbm_power_full$bestModel$fittedSBM$covarParam

## ----power_sampling------------------------------------------------------
nWar <- nrow(beligerent_adjacency)
parameters_sample <- 600
sampleNet_power_miss <- missSBM::sample(
   beligerent_adjacency,
   sampling = "covar-node",
   parameters = parameters_sample, covariates = list(beligerent_power), intercept = -2
  )
boxplot(1/(1 + exp(-cbind(1,beligerent_power) %*% c(-2, parameters_sample))) ~ sampleNet_power_miss$observedNodes, ylab="mil power",xlab = "observed node")
plot(sampleNet_power_miss)

## ----fit power missing---------------------------------------------------
collection_sbm_power_miss <- estimate(sampleNet_power_miss, vBlocks = vBlocks, sampling = "covar-node", useCovariates = TRUE)

## ----estimated parameters sample-----------------------------------------
collection_sbm_power_miss$bestModel$fittedSampling$parameters

## ----estimated parameters SBM--------------------------------------------
collection_sbm_power_miss$bestModel$fittedSBM$covarParam

## ----trade---------------------------------------------------------------
trade <- beligerent_trade 
trade[is.na(trade)] <- 0
trade <- trade + t(trade)
trade <- log(trade + 1)
diag(trade) = 0

## ----samptrade-----------------------------------------------------------
parameters_sample <- 1
sampleNet_trade_miss <- missSBM::sample(beligerent_adjacency, sampling = "covar-dyad", parameters = parameters_sample, covariates = list(trade), intercept = -2)
plot(sampleNet_trade_miss)

## ----estimate trade------------------------------------------------------
collection_sbm_trade_miss <- estimate(sampleNet_trade_miss ,vBlocks = vBlocks, sampling = "covar-dyad", useCovariates  = TRUE)
collection_sbm_trade_miss$bestModel$fittedSampling$parameters
collection_sbm_trade_miss$bestModel$fittedSBM$covarParam

