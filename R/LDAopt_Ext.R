#' PP optimization using LDA index same as PPtree
#' 
#' Find the q-dimensional optimal projection using LDA projectin pursuit index
#' @title PP optimization using LDA index
#' @usage LDAopt_Ext(origclass, origdata, q = 1, weight = TRUE,...) 
#' @param origclass class information vector of data
#' @param origdata data matrix without class information
#' @param q dimension of projection vector
#' @param weight weight flag in LDA index
#' @param ... arguments to be passed to methods
#' @return indexbest maximum LDA index value
#' @return projbest optimal q-dimensional projection matrix
#' @return origclass original class information vector
#' @return origdata  original data matrix  without class information
#' @references Lee, EK., Cook, D., Klinke, S., and Lumley, T.(2005) 
#' Projection Pursuit for Exploratory Supervised Classification, 
#' Journal of Computational and Graphical Statistics, 14(4):831-846.
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