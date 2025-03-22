#' Projection pursuit classification tree with random variable selection in each split
#' 
#' Find tree structure using various projection pursuit indices of classification in each split.
#' @usage PPtree_splitMOD(form, data, PPmethod='LDA', 
#' size.p=1,  lambda=0.1, entro ,entroindiv,...) 
#' @param form A character with the name of the class variable.
#' @param data Data frame with the complete data set.
#' @param PPmethod index to use for projection pursuit: 'LDA', 'PDA'
#' @param size.p proportion of variables randomly sampled in each split, default is 1, returns a PPtree.
#' @param lambda penalty parameter in PDA index and is between 0 to 1 . If \code{lambda = 0}, no penalty parameter is added and the PDA index is the same as LDA index. If \code{lambda = 1} all variables are treated as uncorrelated. The default value is \code{lambda = 0.1}.
#' @param entro TRUE, compute the entropy method
#' @param entroindiv TRUE, compute the entropy for each obs...clarify this
#' @param ... arguments to be passed to methods
#' @return An object of class \code{PPtreeclass} with components
#' \item{Tree.Struct}{Tree structure of projection pursuit classification tree}
#' \item{projbest.node}{1-dim optimal projections of each split node}
#' \item{splitCutoff.node}{cutoff values of each split node}
#' \item{origclass}{original class} 
#' \item{origdata}{original data}
#' @references Lee, YD, Cook, D., Park JW, and Lee, EK (2013) 
#' PPtree: Projection pursuit classification tree, 
#' Electronic Journal of Statistics, 7:1369-1386.
#' @useDynLib PPtreeExt
#' @importFrom Rcpp evalCpp
#' @export
#' @keywords tree
#' @examples
#' #crab data set
#' \dontrun{
#' train<- sample(1:200,150)
#' Tree.crab <- PPtree_splitMOD("Type~.", data = PPforest::crab[train, ],
#'  PPmethod = "LDA", size.p = 1, entro = TRUE,entroindiv=FALSE)
#' Tree.crab
#' 
#'  Tree.result <- PPtreeViz::PPTreeclass(Type~.,data = PPforest::crab[train,],"LDA")
#' Tree.result
#' 
#' PPtreeViz::PPclassify(Tree.result,PPforest::crab[-train,-1],1,crab[-train,1])
#' 
#' Tree.iris <- PPtree_splitMOD("Species~.", data = iris, PPmethod = "LDA", 
#' size.p = 1, entro=TRUE, entroindiv = FALSE)
#' Tree.iris}
PPtree_splitMOD <- function(form, data,  PPmethod = "LDA", size.p = 1,  lambda = 0.1, entro= TRUE, entroindiv = FALSE,...) {
  
     formula <- stats::as.formula(form)
     mf <- stats::model.frame(formula, data = data)
     origclass <- stats::model.response(mf)
    
    cls <- all.vars(formula)[[1]]
    
    origdata <- data[ , -which(colnames(data)%in%cls)]
    origdata <- as.matrix(origdata)
    pp <- ncol(origdata)
    origclass <- as.factor(origclass) #just change numeric here


    g <- table(origclass)
    G <- length(g)
    
    # Tree structure (see optim_index.cpp)
    Tree.final <- treeconstructMOD(origclass, origdata, Treestruct = cbind( 1:(2*G - 1), matrix(0, ncol = 4, nrow = 2*G-1) ), 
                  id = 0,  rep = 1, rep1 = 2, rep2 = 1, projbestnode = matrix(0, ncol = pp, nrow = 1), 
                  splitCutoffnode = matrix(0, ncol = 8, nrow = 1), PPmethod, lambda, size.p, entro = entro, entroindiv = entroindiv)
   
     # Tree.final <- treeconstructIND(origclass, origdata, Treestruct =  matrix(0, ncol = 5, nrow = 1) , 
     #                                id = 0,  rep = 1, rep1 = 2, rep2 = 1, projbestnode = matrix(0, ncol = pp, nrow = 1), 
     #                                splitCutoffnode = matrix(0, ncol = 1, nrow = 1), PPmethod, lambda, size.p, entro=FALSE, entroindiv=TRUE, tot =150)
     # 
    
    Tree.Struct <- Tree.final$Treestruct
    colnames(Tree.Struct) <- c("id", "L.node.ID", "R.F.node.ID", 
                               "Coef.ID", "Index")
    #projbest.node <- Tree.final$projbestnode[-1, ]
    
   
    if(nrow(Tree.final$splitCutoffnode) == 2){
      splitCutoff.node <- data.frame(splitCutoffnode = t(Tree.final$splitCutoffnode[-1, ]))
      colnames(splitCutoff.node) <- paste("Rule", 1:8, sep = "")
      projbest.node <- t(as.matrix(Tree.final$projbestnode[-1, ]))
      
    }else{
      splitCutoff.node <- data.frame(splitCutoffnode = Tree.final$splitCutoffnode[-1, ])
      colnames(splitCutoff.node) <- paste("Rule", 1:8, sep = "")
      projbest.node <- Tree.final$projbestnode[-1, ]
    }
    treeobj <- list(Tree.Struct = Tree.Struct, projbest.node = projbest.node, 
                    splitCutoff.node = splitCutoff.node, origclass = origclass, 
                    origdata = origdata)
    
    class(treeobj) <- append(class(treeobj), "PPtreeclass")
    
    return(treeobj)
}
