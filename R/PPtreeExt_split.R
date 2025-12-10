#' Projection Pursuit Classification Tree with Random Variable Selection
#' 
#' Constructs a projection pursuit classification tree using various projection pursuit 
#' indices. Optionally performs random variable selection at each split which can be used to include in a random forests methodology. When \code{size.p = 1}, this reduces to a  PPtree algorithm.
#' @param formula A formula of the form \code{class ~ x1 + x2 + ...} where \code{class} 
#'   is the factor variable containing class labels and \code{x1, x2, ...} are the 
#'   predictor variables.
#' @param data Data frame containing both the class variable and predictor variables.
#' @param PPmethod Character string specifying the projection pursuit index to use. 
#'   Either \code{"LDA"} (Linear Discriminant Analysis, default) or \code{"PDA"} 
#'   (Penalized Discriminant Analysis).
#' @param size.p Numeric value between 0 and 1 specifying the proportion of variables 
#'   to randomly sample at each split. Default is 1, which uses all variables at each 
#'   split (standard PPtree). Values less than 1 introduce randomness similar to 
#'   random forests, which can improve robustness and reduce overfitting.
#' @param lambda Numeric penalty parameter for the PDA index, ranging from 0 to 1. 
#'   When \code{lambda = 0}, no penalty is applied and PDA equals LDA. When 
#'   \code{lambda = 1}, all variables are treated as uncorrelated. Default is 0.1. 
#'   Only used when \code{PPmethod = "PDA"}.
#' @param entro Logical indicating whether to use entropy-based stopping rules for 
#'   tree construction. Default is \code{FALSE}.
#' @param entroindiv Logical indicating whether to compute entropy for each individual 
#'   observation in the 1D projection. Default is \code{FALSE}.
#' @param ... Additional arguments to be passed to internal tree construction methods.
#' @return An object of class \code{"PPtreeclass"}, which is a list containing:
#' \item{Tree.Struct}{A matrix defining the tree structure of the projection pursuit 
#'   classification tree. Each row represents a node with columns: node ID, left child 
#'   node ID, right child node ID (or final class if terminal), coefficient ID, and 
#'   index value.}
#' \item{projbest.node}{A matrix where each row contains the optimal 1-dimensional 
#'   projection coefficients for each split node. The number of columns equals the 
#'   number of predictor variables.}
#' \item{splitCutoff.node}{A data frame containing the cutoff values and splitting 
#'   rules for each split node. Contains 8 rule columns defining the classification 
#'   boundaries.}
#' \item{origclass}{Factor vector of the original class labels from the input data.}
#' \item{origdata}{Matrix of the original predictor variables (without the class variable).}
#' 
#' @details
#' This function extends the standard PPtree algorithm by incorporating random variable 
#' selection at each split, and define the split based on subsetting groups. The algorithm:
#' \enumerate{
#'   \item At each node, randomly samples \code{size.p * 100}\% of the predictor variables
#'   \item Finds the optimal projection using the selected variables and specified index (LDA or PDA)
#'   \item Determines a cutpoint based on entropy splitting if entropy parameters are set
#'   \item Recursively splits the data until stopping criteria are met
#' }
#' 
#' The \code{entro} parameter enables entropy-based stopping rules that halt splitting 
#' when nodes become sufficiently pure or small. The \code{entroindiv} parameter computes 
#' entropy at the individual observation level in the projected space, which can provide 
#' more refined splitting decisions.
#' 
#' When \code{size.p = 1}, all variables are used at each split and the function 
#' behaves as a standard PPtree. Values of \code{size.p < 1} introduce randomness 
#' that can improve model robustness, especially for high-dimensional data or when 
#' building ensemble models.
#' 
#' @references 
#' Lee, YD, Cook, D., Park JW, and Lee, EK (2013) 
#' PPtree: Projection pursuit classification tree, 
#' Electronic Journal of Statistics, 7:1369-1386.
#' 
#' @seealso \code{\link{TreeExt.construct}}, \code{\link{findproj_Ext}}, 
#'   \code{\link{LDAopt_Ext}}, \code{\link{PDAopt_Ext}}
#' @useDynLib PPtreeExt
#' @importFrom Rcpp evalCpp
#' @export
#' @keywords tree
#' @examples
#' data(penguins)
#' penguins <- na.omit(penguins[, -c(2,7, 8)])
#' require(rsample)
#' penguins_spl <- rsample::initial_split(penguins, strata=species)
#' penguins_train <- training(penguins_spl)
#' penguins_test <- testing(penguins_spl)
#' penguins_ppt2 <- PPtreeExt_split(species~bill_len + bill_dep +
#' flipper_len + body_mass, data = penguins_train, PPmethod = "LDA", tot=nrow
#' (penguins_train), tol =  0.5 , entro=TRUE)
PPtreeExt_split <- function(formula, data,  PPmethod = "LDA", size.p = 1,  lambda = 0.1, entro = FALSE, entroindiv = FALSE,...) {
  
     #formula <- stats::as.formula(form)
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
