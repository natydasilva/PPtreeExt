#' Projection Pursuit Optimization Using PDA Index
#' 
#' Finds the q-dimensional optimal projection using the Penalized Discriminant Analysis (PDA) projection pursuit index. This implementation follows the method described in PPtree and 
#' is particularly useful for high-dimensional data (large p, small n).
#' 
#' @title PP optimization using PDA index
#' @param origclass Factor or numeric vector containing the class labels for each observation.
#' @param origdata Numeric matrix or data frame containing the predictor variables without 
#'   class information. Each row represents an observation and each column represents a variable.
#' @param q Integer specifying the dimension of the projection space. Default is 1 for 
#'   1-dimensional projection.
#' @param weight Logical indicating whether to use weighted PDA index calculation. 
#'   Default is \code{TRUE}.
#' @param lambda Numeric penalty parameter for the PDA index. Controls the amount of 
#'   regularization applied. Default is 0.1. Higher values increase regularization, 
#'   which is useful for high-dimensional or collinear data.
#' @param ... Additional arguments to be passed to internal optimization methods.
#' @return An object of class \code{"PPoptim"}, which is a list containing:
#' \item{indexbest}{Numeric value representing the maximum PDA index achieved by the 
#'   optimal projection. Higher values indicate better class separation with appropriate 
#'   regularization.}
#' \item{projbest}{Numeric matrix of optimal projection coefficients with dimensions 
#'   \code{ncol(origdata)} by \code{q}. Each column represents an optimal projection 
#'   direction that maximizes the PDA index for class separation.}
#' \item{origclass}{The original class information vector passed as input, preserved 
#'   for reference.}
#' \item{origdata}{The original data matrix without class information, preserved 
#'   for reference.}
#' 
#' @details
#' The Penalized Discriminant Analysis (PDA) projection pursuit index extends LDA by 
#' incorporating a penalty term, making it particularly suitable for:
#' \itemize{
#'   \item High-dimensional data where the number of variables exceeds the number of observations (p > n)
#'   \item Data with multicollinearity among predictor variables
#'   \item Cases where standard LDA fails due to singular covariance matrices
#' }
#' 
#' The function performs the following steps:
#' \enumerate{
#'   \item Calls \code{PDAopt} to find the optimal q-dimensional projection directions with regularization
#'   \item Evaluates the PDA index for the optimal projection using \code{PDAindex2}
#'   \item Returns both the projection matrix and its associated index value
#' }
#' 
#' The \code{lambda} parameter controls the trade-off between maximizing class separation 
#' and regularization. When \code{weight = TRUE}, the index calculation accounts for 
#' class proportions in the optimization.
#' 
#' @references 
#' Lee, EK, Cook, D. (2010) 
#' A Projection Pursuit Index for Large p Small n Data, 
#' Statistics and Computing, 20:381-392.
#' 
#' @seealso \code{\link{LDAopt_Ext}}, \code{\link{findproj_Ext}}
#' @useDynLib PPtreeExt
#' @importFrom Rcpp evalCpp
#' @export
#' @keywords projection pursuit
PDAopt_Ext <-function(origclass, origdata, q = 1, weight = TRUE, lambda = 0.1, ...) {
    origdata <- as.matrix(origdata)
  
    optVector <- PDAopt(origclass = origclass, origdata=origdata , q=q,  weight=weight, lambda=lambda)
    optindex <-PDAindex2(origclass = origclass, origdata=origdata, proj = optVector, weight=weight, lambda=lambda )
    optobj <-
      list(
        indexbest = optindex,
        projbest = optVector,
        origclass = origclass,
        origdata = origdata
      )
    class(optobj) <- append(class(optobj), "PPoptim")
    return(optobj)

  }
