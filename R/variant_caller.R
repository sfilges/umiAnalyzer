#' Beta distribution
#'
#' Calculates the negative log likelihood for the beta distribution.
#'
#' @param params Non-negative parameters of the Beta distribution.
#' @param data consensus.data table of a UMisample or UMIexperiment object.
#'
#' @importFrom stats dbeta
#' 
#' @noRd
#'
#' @return Negative log-likelihood for beta distribution.
#'
betaNLL <- function(params, data) {
  a <- params[1]
  b <- params[2]
  
  # negative log likelihood for beta
  return(-sum(stats::dbeta(data, shape1 = a, shape2 = b, log = TRUE)))
}

#' callVariants using beta binomial distribution
#'
#' Calculate variant p-values using permutation-based testing. A prior is fitted
#' to model the background error using maximum likelihood estimation of a beta
#' distribution. The maximum likelihood estimate of the beta distribution is then
#' used to define the shape of a beta-binomial distribution used to estimate
#' variant P-Values. This can be interpreted as a probability for a variant to
#' not have arisen by chance.
#'
#' @export
#'
#' @param object A UMIErrorCorrect object.
#' @param minDepth Minimum consensus depth required default is 3
#' @param minCoverage Minimum Coverage to use, default is 100 reads.
#' @param computePrior Should a new distribution be derived from data? Default is FALSE.
#'
#' @importFrom dplyr mutate progress_estimated
#' @importFrom tibble as_tibble
#' @importFrom stats nlm var p.adjust
#' @importFrom utils install.packages
#' @importFrom graphics plot
#'
#' @return Object containing raw and FDR-adjusted P-Values
#' @seealso \code{\link{filterVariants}} on how to filter variants.
#'
#'
callVariants <- function(
  object,
  minDepth = 3,
  minCoverage = 100,
  computePrior = FALSE
) {
  
  if (missing(x = object)) {
    stop("Must provide a umiExperiment object.")
  } else if(!class(object) == "UMIexperiment"){
    stop("Object is not of class UMIexperiment.")
  } else if(!is.numeric(minDepth) || minDepth < 3){
    if(minDepth < 0){
      stop("minDepth must be a positive integer.")
    } else if (minDepth < 3){
      warning("Mindepth is smaller than 3, this may affect error correction.")
    }
  } else if(!is.numeric(minCoverage) || minCoverage < 100){
    if(minCoverage < 0){
      stop("minCoverage must be a positive integer.")
    } else if(minCoverage < 100){
      warning("minCoverage is smaller than 100, this may affect error correction.")
    }
  }
  
  object <- filterUmiObject(
    object = object,
    name = "varCalls",
    minDepth = minDepth,       # Require minDepth
    minCoverage = minCoverage, # Require at least minCoverage cons reads
    minFreq = 0,               # no minimum allele freq
    minCount = 0              # no minimum variant allele count
  )
  
  cons.table <- object@filters["varCalls"][[1]]
  
  a1 <- cons.table$Coverage * cons.table$`Max Non-ref Allele Frequency` # No. of variant alleles
  b1 <- cons.table$Coverage # Total coverage
  
  m <- mean(a1 / b1) # average background count
  v <- var(a1 / b1) # variance of background counts
  
  # Calculate initial values
  a0 <- m * (m * (1 - m) / v - 1)
  b0 <- (1 - m) * (m * (1 - m) / v - 1)
  params0 <- c(a0, b0)
  
  # If computePrior == TRUE fit a new beta binomial background distribution,
  # otherwise use pre-computed values (default)
  if(computePrior){
    fit <- stats::nlm(betaNLL, params0, a0 / b0)
    a <- fit$estimate[1]
    b <- fit$estimate[2]
    pval <- NULL
    
  } else {
    a <- 2.168215069116764
    b <- 3531.588541594945
    pval <- NULL
  }
  
  for (i in 1:length(a1)) { # for each named amplicon position:
    #r1 <- VGAM::rbetabinom.ab(10000, b1[i], shape1 = a, shape2 = b) # Calculate probability of success
    r1 <- beta_binom(10000, b1[i], shape1 = a, shape2 = b)
    pval[i] <- sum(r1 > a1[i]) / 10000 # Estimate p value of variant
  }
  
  padj <- stats::p.adjust(p = pval, method = 'fdr')
  
  cons.table <- dplyr::mutate(cons.table, pval = pval)
  cons.table <- dplyr::mutate(cons.table, p.adjust = padj)
  
  object@variants <- cons.table
  
  object <- addMetaData(
    object = object,
    attributeName = "varCalls", "varCalls"
  )
  
  return(object)
}