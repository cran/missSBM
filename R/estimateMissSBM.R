#' Estimation of simple SBMs with missing data
#'
#' Variational EM inference of Stochastic Block Models indexed by block number from a partially observed network.
#'
#' @param adjacencyMatrix The N x N adjacency matrix of the network data. If \code{adjacencyMatrix} is symmetric,
#' we assume an undirected network with no loop; otherwise the network is assumed to be directed.
#' @param vBlocks The vector of number of blocks considered in the collection.
#' @param sampling The model used to described the process that originates the missing data:
#' MAR designs ("dyad", "node","covar-dyad","covar-node","snowball") and NMAR designs
#' ("double-standard", "block-dyad", "block-node" , "degree") are available. See details.
#' @param covariates A list with M entries (the M covariates). If the covariates are node-centered, each entry of \code{covariates}
#' must be a size-N vector;  if the covariates are dyad-centered, each entry of \code{covariates} must be N x N matrix.
#' @param control a list of parameters controlling advanced features. See details.
#'
#' @return Returns an R6 object with class \code{\link{missSBM_collection}}.
#'
#' @details The list of parameters \code{control} tunes more advanced features, such as the
#' initialization, how covariates are handled in the model, and the variational EM algorithm:
#'  \itemize{
#'  \item{"useCovSBM": }{logical. If \code{covariates} is not null, should they be used for the
#'         for the SBM inference (or just for the sampling)? Default is TRUE.}
#'  \item{"clusterInit": }{Initial method for clustering: either a character in "hierarchical", "spectral"
#'         or "kmeans", or a list with \code{length(vBlocks)} vectors, each with size
#'         \code{ncol(adjacencyMatrix)},  providing a user-defined clustering. Default is "spectral".}
#'  \item{"similarity": }{An R x R -> R function to compute similarities between node covariates. Default is
#'         \code{missSBM:::l1_similarity}, that is, -abs(x-y). Only relevant when the covariates are node-centered
#'         (i.e. \code{covariates} is a list of size-N vectors).}
#'  \item{"threshold": }{V-EM algorithm stops stop when an optimization step changes the objective function
#'         by less than threshold. Default is 1e-3.}
#'  \item{"maxIter": }{V-EM algorithm stops when the number of iteration exceeds maxIter.
#'        Default is 100 with no covariate, 50 otherwise.}
#'  \item{"fixPointIter": }{number of fix-point iterations in the V-E step.
#'        Default is 5 with no covariate, 2 otherwise.}
#'  \item{"cores": }{integer for number of cores used. Default is 1.}
#'  \item{"trace": }{integer for verbosity (0, 1, 2). Default is 1. Useless when \code{cores} > 1}
#' }
#'
#' @details The different sampling designs are split into two families in which we find dyad-centered and
#' node-centered samplings. See \doi{10.1080/01621459.2018.1562934} for a complete description.
#' \itemize{
#' \item Missing at Random (MAR)
#'   \itemize{
#'     \item{"dyad": parameter = p = Prob(Dyad(i,j) is observed)}
#'     \item{"node": parameter = p = Prob(Node i is observed)}
#'     \item{"covar-dyad": parameter = beta in R^M, such that Prob(Dyad (i,j) is observed) = logistic(parameter' covarArray (i,j, .))}
#'     \item{"covar-node": parameter = nu in R^M such that Prob(Node i is observed)  = logistic(parameter' covarMatrix (i,)}
#'     \item{"snowball": parameter = number of waves with Prob(Node i is observed in the 1st wave)}
#'   }
#' \item Not Missing At Random (NMAR)
#'   \itemize{
#'     \item{"double-standard": parameter = (p0,p1) with p0 = Prob(Dyad (i,j) is observed | the dyad is equal to 0), p1 = Prob(Dyad (i,j) is observed | the dyad is equal to 1)}
#'     \item{"block-node": parameter = c(p(1),...,p(Q)) and p(q) = Prob(Node i is observed | node i is in cluster q)}
#'     \item{"block-dyad": parameter = c(p(1,1),...,p(Q,Q)) and p(q,l) = Prob(Edge (i,j) is observed | node i is in cluster q and node j is in cluster l)}
#'     \item{"degree": parameter = c(a,b) and logit(a+b*degree(i)) = Prob(Node i is observed | Degree(i))}
#'   }
#' }
#' @seealso \code{\link{observeNetwork}}, \code{\link{missSBM_collection}} and \code{\link{missSBM_fit}}.
#' @examples
#' ## SBM parameters
#' N <- 150 # number of nodes
#' Q <- 3   # number of clusters
#' pi <- rep(1,Q)/Q     # block proportion
#' theta <- list(mean = diag(.45,Q) + .05 ) # connectivity matrix
#'
#' ## Sampling parameters
#' samplingParameters <- .5 # the sampling rate
#' sampling  <- "dyad"      # the sampling design
#'
#' ## generate a undirected binary SBM with no covariate
#' sbm <- sbm::sampleSimpleSBM(N, pi, theta)
#'
#' ## Sample some dyads data + Infer SBM with missing data
#' collection <-
#'    observeNetwork(sbm$netMatrix, sampling, samplingParameters) %>%
#'    estimateMissSBM(vBlocks = 1:5, sampling = sampling)
#' collection$ICL
#' coef(collection$bestModel$fittedSBM, "connectivity")
#'
#' myModel <- collection$bestModel
#' plot(myModel, "network")
#' coef(myModel, "sampling")
#' coef(myModel, "connectivity")
#' predict(myModel)[1:5, 1:5]
#' fitted(myModel)[1:5, 1:5]
#'
#' @import R6 parallel
#' @export
estimateMissSBM <- function(adjacencyMatrix, vBlocks, sampling, covariates = NULL, control = list()) {

  ## Sanity checks
  stopifnot(sampling %in% available_samplings)
  stopifnot(is.numeric(vBlocks))
  stopifnot(is.character(sampling))

  ## If no covariate is provided, you cannot ask for using them
  if (is.null(covariates)) control$useCovSBM <- FALSE
  ## If nothing specified by the user, use covariates by default
  else if (is.null(control$useCovSBM)) control$useCovSBM <- TRUE

  ## Defaut control parameters overwritten by user specification
  ctrl <- list(threshold = 1e-3, trace = 1, cores = 1, clusterInit = "hierarchical", similarity = l1_similarity)
  if (control$useCovSBM) {
    stopifnot(sampling %in% available_samplings_covariates)
    ctrl <- c(ctrl,list(maxIter = 50, fixPointIter = 2))
  } else {
    ctrl <- c(ctrl, list(maxIter = 100, fixPointIter = 5))
  }
  ctrl[names(control)] <- control

  ## Prepare network data for estimation with missing data
  partlyObservedNet <- partlyObservedNetwork$new(adjacencyMatrix, covariates, ctrl$similarity)

  ## Instantiate the collection of missSBM_fit
  myCollection <- missSBM_collection$new(
      partlyObservedNet  = partlyObservedNet,
      vBlocks            = vBlocks,
      sampling           = sampling,
      clusterInit        = ctrl$clusterInit,
      cores              = ctrl$cores,
      trace              = (ctrl$trace > 0),
      useCov             = ctrl$useCovSBM
  )

  ## Launch estimation of each missSBM_fit
  myCollection$estimate(ctrl)

  ## Return the collection of adjusted missSBM_fit
  myCollection
}
