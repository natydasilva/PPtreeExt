#' Find Optimal Projection for Class Separation
#' 
#' Finds an optimal 1D projection of multivariate data that best separates classes 
#' using Linear Discriminant Analysis (LDA) or Penalized Discriminant Analysis (PDA), 
#' then determines a cutpoint for classification based on entropy splitting.
#' @param origclass Factor or numeric vector containing the class labels for each observation.
#' @param origdata Numeric matrix or data frame containing the predictor variables. Each row represents an observation and each column represents a variable.
#' @param PPmethod Character string specifying the projection pursuit method. 
#'   Either \code{"LDA"} (Linear Discriminant Analysis, default) or \code{"PDA"} 
#'   (Penalized Discriminant Analysis).
#' @param q Integer specifying the dimension of the projected data. Default is 1 
#'   for 1D projection.
#' @param weight Logical indicating whether to use weighted LDA index calculation. 
#'   Default is \code{TRUE}.
#' @param lambda Numeric penalty parameter for the PDA method. Default is 0.1. 
#'   Only used when \code{PPmethod = "PDA"}.
#' @return A list with the following components:
#' \item{Index}{Numeric value representing the optimization criterion achieved by 
#'   the best projection. Higher values indicate better class separation.}
#' \item{Alpha}{Numeric vector of length \code{ncol(origdata)} containing the optimal 
#'   projection direction coefficients. This vector defines the linear combination 
#'   of original variables that maximizes class separation.}
#' \item{C}{Numeric scalar representing the optimal cutpoint (threshold) on the 
#'   projected data. This value is determined using entropy-based splitting and 
#'   divides observations into two groups for classification.}
#' \item{IOindexL}{Logical vector of length \code{nrow(origdata)} indicating which 
#'   observations have projected values less than or equal to the cutpoint C 
#'   (\code{projdata <= C}). These observations are assigned to the left node/class.}
#' \item{IOindexR}{Logical vector of length \code{nrow(origdata)} indicating which 
#'   observations have projected values greater than the cutpoint C (\code{projdata > C}). 
#'   These observations are assigned to the right node/class.}
#' @details
#' This function performs projection pursuit to find a one-dimensional projection 
#' that optimally separates classes in multivariate data. The process involves:
#' \enumerate{
#'   \item Finding the optimal projection direction using either LDA or PDA
#'   \item Projecting all observations onto this direction
#'   \item Determining an optimal cutpoint using entropy-based splitting
#'   \item Creating binary classification indicators based on the cutpoint
#' }
#' 
#' The cutpoint is calculated to minimize the weighted entropy of the resulting split. 
#' In edge cases where the cutpoint equals the maximum projected value, the function 
#' uses the second-largest value to ensure a valid split.
#' 
#' @note The vectors \code{IOindexL} and \code{IOindexR} are complementary 
#'   (mutually exclusive and exhaustive), meaning every observation is assigned 
#'   to exactly one group.
#'   
#' @references Lee, YD, Cook, D., Park JW, and Lee, EK (2013) 
#' PPtree: Projection Pursuit Classification Tree, 
#' Electronic Journal of Statistics, 7:1369-1386.
#' 
#' @useDynLib PPtreeExt
#' @importFrom Rcpp evalCpp
#' @export
findproj_Ext <- function(origclass, origdata, PPmethod = "LDA", q = 1, weight = TRUE, lambda = .1) {
  

  
  if(PPmethod == "LDA"){
    idx <- LDAopt_Ext(origclass, origdata, q , weight)
  }else{
    idx <- PDAopt_Ext(origclass, origdata, q , weight,lambda)
  }
  projdata = apply(origdata, 1, function(x) sum(x*idx$projbest) )
  
  #projdata = as.vector(origdata %*% idx$projbest)
  
  cp <- split_entro(origclass, projdata)

  #if ( cp == max(projdata) ) cp <- sort(projdata)[ length(projdata) - 1]
  if (abs(cp - max(projdata)) < 1e-10) {
    cp <- sort(projdata, decreasing = TRUE)[2]
  }

  
  list(
    Index  = idx$indexbest,
    Alpha = idx$projbest,
    C = cp, 
    IOindexL = projdata <= cp,
    IOindexR = projdata > cp
  )
}