#' Estimation of simple SBMs with missing data
#'
#' Variational EM inference of Stochastic Block Models indexed by block number from a partially observed network.
#'
#' @param adjacencyMatrix The N x N adjacency matrix of the network data. If \code{adjacencyMatrix} is symmetric,
#' we assume an undirected network with no loop; otherwise the network is assumed to be directed.
#' @param vBlocks The vector of number of blocks considered in the collection.
#' @param sampling The model used to described the process that originates the missing data:
#' MAR designs ("dyad", "node","covar-dyad","covar-node","snowball") and MNAR designs
#' ("double-standard", "block-dyad", "block-node" , "degree") are available. See details.
#' @param covariates An optional list with M entries (the M covariates). If the covariates are node-centered, each entry of \code{covariates}
#' must be a size-N vector;  if the covariates are dyad-centered, each entry of \code{covariates} must be N x N matrix.
#' @param control a list of parameters controlling advanced features. See details.
#'
#' @return Returns an R6 object with class \code{\link{missSBM_collection}}.
#'
#' @details Internal functions use \code{future_lapply}, so set your plan to \code{'multisession'} or
#' \code{'multicore'} to use several cores/workers.
#' The list of parameters \code{control} tunes more advanced features, such as the
#' initialization, how covariates are handled in the model, and the variational EM algorithm:
#'  * useCov logical. If \code{covariates} is not null, should they be used for the
#'         for the SBM inference (or just for the sampling)? Default is TRUE.
#'  * clusterInit Initial method for clustering: either a character ("spectral")
#'         or a list with \code{length(vBlocks)} vectors, each with size  \code{ncol(adjacencyMatrix)},
#'         providing a user-defined clustering. Default is "spectral".
#'  similarity An R x R -> R function to compute similarities between node covariates. Default is
#'         \code{l1_similarity}, that is, -abs(x-y). Only relevant when the covariates are node-centered
#'         (i.e. \code{covariates} is a list of size-N vectors).
#'  * threshold V-EM algorithm stops stop when an optimization step changes the objective function or the parameters
#'         by less than threshold. Default is 1e-2.
#'  * maxIter V-EM algorithm stops when the number of iteration exceeds maxIter.
#'        Default is 50.
#'  * fixPointIter number of fix-point iterations in the V-E step. Default is 3.
#'  * exploration character indicating the kind of exploration used among "forward", "backward", "both" or "none". Default is "both".
#'  * iterates integer for the number of iterations during exploration. Only relevant when \code{exploration} is different from "none". Default is 1.
#'  * trace logical for verbosity. Default is TRUE.
#'
#' @details The different sampling designs are split into two families in which we find dyad-centered and
#' node-centered samplings. See \doi{10.1080/01621459.2018.1562934} for a complete description.
#' * Missing at Random (MAR)
#'    * dyad parameter = p = Prob(Dyad(i,j) is observed)
#'    * node parameter = p = Prob(Node i is observed)
#'    * covar-dyad": parameter = beta in R^M, such that Prob(Dyad (i,j) is observed) = logistic(parameter' covarArray (i,j, .))
#'    * covar-node": parameter = nu in R^M such that Prob(Node i is observed)  = logistic(parameter' covarMatrix (i,)
#'    * snowball": parameter = number of waves with Prob(Node i is observed in the 1st wave)
#' * Missing Not At Random (MNAR)
#'    * double-standard parameter = (p0,p1) with p0 = Prob(Dyad (i,j) is observed | the dyad is equal to 0), p1 = Prob(Dyad (i,j) is observed | the dyad is equal to 1)
#'    * block-node parameter = c(p(1),...,p(Q)) and p(q) = Prob(Node i is observed | node i is in cluster q)
#'    * block-dyad parameter = c(p(1,1),...,p(Q,Q)) and p(q,l) = Prob(Edge (i,j) is observed | node i is in cluster q and node j is in cluster l)
#     * degree parameter = c(a,b) and logit(a+b*degree(i)) = Prob(Node i is observed | Degree(i))
#'
#' @seealso \code{\link{observeNetwork}}, \code{\link{missSBM_collection}} and \code{\link{missSBM_fit}}.
#' @examples
#' ## SBM parameters
#' N <- 100 # number of nodes
#' Q <- 3   # number of clusters
#' pi <- rep(1,Q)/Q     # block proportion
#' theta <- list(mean = diag(.45,Q) + .05 ) # connectivity matrix
#'
#' ## Sampling parameters
#' samplingParameters <- .75 # the sampling rate
#' sampling  <- "dyad"      # the sampling design
#'
#' ## generate a undirected binary SBM with no covariate
#' sbm <- sbm::sampleSimpleSBM(N, pi, theta)
#'
#' ## Uncomment to set parallel computing with future
#' ## future::plan("multicore", workers = 2)
#'
#' ## Sample some dyads data + Infer SBM with missing data
#' collection <-
#'    observeNetwork(sbm$networkData, sampling, samplingParameters) %>%
#'    estimateMissSBM(vBlocks = 1:4, sampling = sampling)
#' plot(collection, "monitoring")
#' plot(collection, "icl")
#'
#' collection$ICL
#' coef(collection$bestModel$fittedSBM, "connectivity")
#'
#' myModel <- collection$bestModel
#' plot(myModel, "expected")
#' plot(myModel, "imputed")
#' plot(myModel, "meso")
#' coef(myModel, "sampling")
#' coef(myModel, "connectivity")
#' predict(myModel)[1:5, 1:5]
#'
#' @export
estimateMissSBM <- function(adjacencyMatrix, vBlocks, sampling, covariates = list(), control = list()) {

  ## Sanity checks
  stopifnot(sampling %in% available_samplings)
  stopifnot(is.numeric(vBlocks))
  stopifnot(is.character(sampling))
  stopifnot(is.list(covariates))

  ## Default control parameters overwritten by user specification
  ctrl <- list(
    threshold = 1e-2, trace = TRUE, imputation = "median", similarity = l1_similarity, useCov = TRUE,
    maxIter = 50, fixPointIter = 3, iterates = 1, exploration = "both", clusterInit = NULL
    )
  ctrl[names(control)] <- control
  ## If no covariate is provided, you cannot ask for using them
  if (length(covariates) == 0) ctrl$useCov <- FALSE
  if (ctrl$useCov) stopifnot(sampling %in% available_samplings_covariates)

  ## Prepare network data for estimation with missing data
  partlyObservedNet <- partlyObservedNetwork$new(adjacencyMatrix, covariates, ctrl$similarity)
  clusterInit <- ctrl$clusterInit
  if (is.null(clusterInit)) clusterInit <- partlyObservedNet$clustering(vBlocks, ctrl$imputation)

  ## Instantiate the collection of missSBM_fit
  myCollection <- missSBM_collection$new(
      partlyObservedNet  = partlyObservedNet,
      sampling           = sampling,
      clusterInit        = clusterInit,
      control            = ctrl
  )

  ## Launch estimation of each missSBM_fit
  myCollection$estimate(ctrl)

  ## Looking for better models around
  myCollection$explore(ctrl)

  ## Return the collection of adjusted missSBM_fit
  myCollection
}

