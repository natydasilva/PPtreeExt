#' PP optimization using PDA index same as PPtree
#' 
#' Find the q-dimensional optimal projection using PDA projectin pursuit index
#' @title PP optimization using PDA index
#' @usage PDAopt_Ext(origclass,origdata, q = 1, weight = TRUE, lambda = 0.1,...) 
#' @param origclass class information vector of data
#' @param origdata data matrix without class information
#' @param q dimension of projection vector
#' @param weight weight flag in PDA index
#' @param lambda lambda in PDA index
#' @param ... arguments to be passed to methods
#' @return indexbest maximum PDA index value
#' @return projbest optimal q-dimensional projection matrix
#' @return origclass original class information vector
#' @return origdata  original data matrix without class information
#' @references Lee, EK, Cook, D.(2010) 
#' A Projection Pursuit Index for Large p Small n Data, 
#' Statistics and Computing, 20:381-392.
#' @useDynLib PPtreeExt
#' @importFrom Rcpp evalCpp
#' @export
#' @keywords projection pursuit
#' @examples
#' data(penguins)
#' penguins <- na.omit(penguins[, -c(2,7)])
#' penguins_pda_proj <- PDAopt_Ext(origclass = penguins[,1], origdata = penguins[,-1], 
#' weight=TRUE, q=1, lambda=0.1)
#' penguins_pda_proj$indexbest
#' penguins_pda_proj$projbest
PDAopt_Ext <-function(origclass, origdata, q = 1, weight = TRUE, lambda = 0.1, ...) {
    origdata <- as.matrix(origdata)
  
    optVector <- PDAopt(origclass = origclass, origdata=origdata , q=1,  weight=weight, lambda=lambda)
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
