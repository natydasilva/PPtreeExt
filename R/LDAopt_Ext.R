#' Projection Pursuit Optimization Using LDA Index
#' 
#' Finds the q-dimensional optimal projection using the Linear Discriminant Analysis (LDA) 
#' projection pursuit index. This implementation follows the method described in PPtree.
#' @title PP optimization using LDA index
#' @param origclass Factor or numeric vector containing the class labels for each observation.
#' @param origdata Numeric matrix or data frame containing the predictor variables without 
#'   class information. Each row represents an observation and each column represents a variable.
#' @param q Integer specifying the dimension of the projection space. Default is 1 for 
#'   1-dimensional projection.
#' @param weight Logical indicating whether to use weighted LDA index calculation. 
#'   Default is \code{TRUE}.
#' @param ... Additional arguments to be passed to internal optimization methods.
#' 
#' @return An object of class \code{"PPoptim"}, which is a list containing:
#' \item{indexbest}{Numeric value representing the maximum LDA index achieved by the 
#'   optimal projection. Higher values indicate better class separation.}
#' \item{projbest}{Numeric matrix of optimal projection coefficients with dimensions 
#'   \code{ncol(origdata)} by \code{q}. Each column represents an optimal projection 
#'   direction that maximizes the LDA index for class separation.}
#' \item{origclass}{The original class information vector passed as input, preserved 
#'   for reference.}
#' \item{origdata}{The original data matrix without class information, preserved 
#'   for reference.}
#' @details
#' The LDA projection pursuit index measures class separation by maximizing the ratio 
#' of between-class variance to within-class variance in the projected space. This 
#' function:
#' \enumerate{
#'   \item Calls \code{LDAopt} to find the optimal q-dimensional projection directions
#'   \item Evaluates the LDA index for the optimal projection using \code{LDAindex2}
#'   \item Returns both the projection matrix and its associated index value
#' }
#' 
#' When \code{weight = TRUE}, the index calculation accounts for class proportions, 
#' giving appropriate weight to each class in the optimization.
#' 
#' @references 
#' Lee, EK., Cook, D., Klinke, S., and Lumley, T. (2005) 
#' Projection Pursuit for Exploratory Supervised Classification, 
#' Journal of Computational and Graphical Statistics, 14(4):831-846.
#' 
#' @seealso \code{\link{PDAopt_Ext}}, \code{\link{findproj_Ext}}
#' 
#' @useDynLib PPtreeExt
#' @importFrom Rcpp evalCpp
#' @export
LDAopt_Ext <- function(origclass, origdata, q = 1, weight = TRUE, ...) {
  
  origdata <- as.matrix(origdata)
  optVector <- LDAopt(origclass = origclass, origdata = origdata, q = q, weight = weight )
  optindex <- LDAindex2(origclass = origclass, origdata = origdata, proj = optVector, weight = weight)

   optobj <- list(indexbest = optindex, projbest = optVector, 
                  origclass = origclass, origdata = origdata)
  class(optobj) <- append(class(optobj), "PPoptim")
  return(optobj)
}