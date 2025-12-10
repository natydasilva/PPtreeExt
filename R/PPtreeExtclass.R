#' Projection Pursuit Classification Tree with Extensions
#' 
#' Constructs a projection pursuit classification tree using various projection pursuit 
#' indices (LDA or PDA) at each split. This extended version includes customizable 
#' stopping rules based on entropy and node size criteria.
#' @title Projection pursuit classification tree 
#' @usage PPtreeExtclass(formula, data, PPmethod = "LDA", weight = TRUE,
#'                    lambda = 0.1,srule, tot = nrow(data), tol = 0.5,...) 
#' @param formula An object of class \code{"formula"} of the form \code{class ~ x1 + x2 + ...} 
#'   where \code{class} is the factor variable containing class labels and 
#'   \code{x1, x2, ...} are the predictor variables. Interaction terms (using \code{*}) 
#'   are not supported.
#' @param data Data frame containing both the class variable and predictor variables 
#'   specified in the formula.
#' @param PPmethod Character string specifying the projection pursuit index to use. 
#'   Either \code{"LDA"} (Linear Discriminant Analysis, default) or \code{"PDA"} 
#'   (Penalized Discriminant Analysis).
#' @param weight Logical indicating whether to use weighted index calculation in LDA 
#'   and PDA. When \code{TRUE} (default), class proportions are accounted for in the 
#'   optimization.
#' @param lambda Numeric penalty parameter for the PDA index, ranging from 0 to 1. 
#'   Default is 0.1. Only used when \code{PPmethod = "PDA"}.
#' @param srule Logical flag for stopping rule. If \code{TRUE} (default), uses entropy-based 
#'   and size-based stopping criteria. If \code{FALSE}, stops only when nodes are pure 
#'   (single class) or empty.
#' @param tot Integer specifying the total number of observations in the original dataset. 
#'   Default is \code{nrow(data)}. Used in conjunction with stopping rules to determine 
#'   minimum node sizes.
#' @param tol Numeric tolerance value for the entropy-based stopping rule. Nodes with 
#'   entropy below this threshold will not be split further. Default is 0.5. Lower 
#'   values create deeper trees.
#' @param ... Additional arguments to be passed to internal tree construction methods.
#' @return An object of class \code{c("PPtreeExtclass", "PPtreeclass")}, which is a list containing:
#' \item{Tree.Struct}{A matrix defining the tree structure. Each row represents a node 
#'   with 5 columns: node ID, left child node ID, right/final node ID (class label if 
#'   terminal node), coefficient ID (projection index), and optimization index value.}
#' \item{projbest.node}{A matrix where each row contains the optimal 1-dimensional 
#'   projection coefficients for each split node. Each row has length equal to 
#'   \code{ncol(origdata)}, defining the projection direction used at that node.}
#' \item{splitCutoff.node}{A numeric vector or matrix containing the cutoff values 
#'   (thresholds) used at each split node for classification decisions.}
#' \item{origclass}{Factor vector of the original class labels from the input data.}
#' \item{origdata}{Matrix of the original predictor variables (without the class variable).}
#' \item{terms}{The terms object from the model frame, preserving the formula structure.}
#' 
#' @details
#' This function builds a binary classification tree where each split is determined by 
#' finding an optimal projection of the data onto a one-dimensional space using either 
#' LDA or PDA indices. The algorithm works as follows:
#' 
#' \subsection{Tree Construction Process}{
#' \enumerate{
#'   \item At each node, find the optimal 1D projection that best separates classes
#'   \item Project the data onto this direction and find an optimal cutpoint
#'   \item Split observations based on the cutpoint into left and right child nodes
#'   \item Recursively repeat until stopping criteria are met
#' }
#' }
#' 
#' \subsection{Stopping Rules}{
#' When \code{srule = TRUE}, a node stops splitting if any of the following conditions hold:
#' \itemize{
#'   \item The node is pure (contains only one class)
#'   \item The node contains fewer than 5\% of the total observations (\code{n/tot <= 0.05})
#'   \item The node entropy is below the tolerance threshold (\code{entropy < tol})
#' }
#' 
#' When \code{srule = FALSE}, splitting only stops for pure or empty nodes, potentially 
#' creating deeper, more complex trees.
#' }
#' 
#' \subsection{Projection Methods}{
#' \itemize{
#'   \item \strong{LDA}: Suitable for most classification problems with moderate dimensionality
#'   \item \strong{PDA}: Recommended for high-dimensional data (p > n) or data with multicollinearity
#' }
#' }
#' 
#' The \code{tol} parameter controls tree complexity: smaller values allow more splits 
#' (deeper trees with potentially better training accuracy but higher risk of overfitting), 
#' while larger values create simpler trees (better generalization but potentially 
#' underfitting).
#' 
#' @note This function does not support interaction terms in the formula. Use only 
#' additive terms (e.g., \code{y ~ x1 + x2}) and not multiplicative terms (e.g., 
#' \code{y ~ x1 * x2}).
#' 
#' @references 
#' Lee, YD, Cook, D., Park JW, and Lee, EK (2013) 
#' PPtree: Projection Pursuit Classification Tree, 
#' Electronic Journal of Statistics, 7:1369-1386.
#' 
#' @seealso \code{\link{TreeExt.construct}}, \code{\link{PPtreeExt_split}}, 
#'   \code{\link{findproj_Ext}}, \code{\link{predict.PPtreeExtclass}}
#' @export
#' @keywords tree
#' @examples
#' set.seed(666)
#' data(penguins)
#' penguins <- na.omit(penguins[, -c(2,7, 8)])
#' require(rsample)
#' penguins_spl <- rsample::initial_split(penguins, strata=species)
#' penguins_train <- training(penguins_spl)
#' penguins_test <- testing(penguins_spl)
#' penguins_ppt <- PPtreeExtclass(species~bill_len + bill_dep +
#' flipper_len + body_mass, data = penguins_train, PPmethod = "LDA", tot=nrow
#' (penguins_train), tol =  0.2 , srule = TRUE)
PPtreeExtclass <- function(formula,
                           data,
                           PPmethod = "LDA",
                           weight = TRUE,
                           lambda = 0.1,
                           srule = TRUE,
                           tot = nrow(data),
                           tol = 0.5,
                           ...) {
  data <- data.frame(data)
  
  Call <- match.call()
  indx <- match(c("formula", "data"), names(Call), nomatch = 0L)
  if (indx[1] == 0L)
    stop("a 'formula' argument is required")
  temp <- Call[c(1L, indx)]
  temp[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(temp)
  Terms <- attr(m, "terms")
  formula <- as.character(formula)
  class.n <- formula[2]
  data.n <- strsplit(formula[3], " \\+ ")[[1]]
  int.flag <- any(strsplit(formula[3], " \\* ")[[1]] == formula[3])
  if (data.n[1] == ".") {
    tot.n <- class.n
  } else{
    tot.n <- c(class.n, data.n)
  }
  if (!int.flag) {
    stop("PPTreeclass cannot treat interaction terms")
  } else if (!sum(duplicated(c(colnames(data), tot.n))[-c(1:ncol(data))]) ==
             length(tot.n)) {
    
  } else{
    origclass <- data[, class.n]
    if (data.n[1] == ".") {
      origdata <- data[, colnames(data) != class.n]
    } else {
      origdata <- data[, data.n, drop = FALSE]
    }
  }
  TOL <- NULL
  origdata <- as.matrix(origdata)
  
  class.table <- table(origclass)
  
  g <- length(class.table)
  
  class.name <- names(class.table)
  
  rm(class.table)
  rm(g)
  rm(class.name)
  
  splitCutoff.node <- NULL
  projbest.node <- NULL
  Tree.Struct <- NULL
  id <- 1
  rep1 <- 2
  rep2 <- 1
  rep <- 1

  Tree.final <- TreeExt.construct(origclass = origclass, origdata = origdata, Tree.Struct = Tree.Struct,
    id = id, rep = rep, rep1 = rep1, rep2 = rep2, projbest.node = projbest.node, splitCutoff.node = splitCutoff.node, PPmethod = PPmethod, lambda = lambda, q = 1, weight = weight, srule = srule,  tot = tot, tol = tol, ...)
  
  
  
  Tree.Struct <- Tree.final$Tree.Struct
  colnames(Tree.Struct) <- c("id", "L.node.ID", "R.F.node.ID", "Coef.ID", "Index")
  projbest.node <- Tree.final$projbest.node
  splitCutoff.node <- Tree.final$splitCutoff.node
  #colnames(splitCutoff.node)<-paste("Rule",1,sep="")
  treeobj <- list(
    Tree.Struct = Tree.Struct,
    projbest.node = projbest.node,
    splitCutoff.node = splitCutoff.node,
    origclass = origclass,
    origdata = origdata,
    terms = Terms
  )
  class(treeobj) <- append(class(treeobj), c("PPtreeExtclass", "PPtreeclass"))
  return(treeobj)
}
