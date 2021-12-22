#' 1D projection for each node partition using entropy (NEW)
#' 
#' @usage findproj_MOD(origclass, origdata, PPmethod="LDA", q=1,weight=TRUE,lambda =.1) 
#' @param origclass original class 
#' @param origdata original data
#' @param PPmethod method for projection pursuit; "LDA", "PDA"
#' @param q  1D proj
#' @param weight weight flag in LDA, PDA and Lr index
#' @param lambda lambda in PDA index
#' @references Lee, YD, Cook, D., Park JW, and Lee, EK(2013) 
#' PPtree: Projection Pursuit Classification Tree, 
#' Electronic Journal of Statistics, 7:1369-1386.
#' @useDynLib PPtreeExt
#' @importFrom Rcpp evalCpp
#' @export
#' @export

findproj_MOD <- function(origclass, origdata, PPmethod="LDA", q=1,weight=TRUE,lambda =.1) {
  
  
  
  cat("\n\n ---- Corrida findproj_MOD PRE PDA_opt---- \n\n")
  cat("dimensión de origdata: ", dim(origdata), "\n")
  class.table <- table(origclass)
  cat("class.table: ", class.table, "\n")
  g <- length(class.table)
  cat("g = ", g, "\n")
  class.name <- names(class.table)
  cat("class.name = ", class.name, "\n")
  p <- ncol(origdata)
  n <- nrow(origdata)
  # comentarios
  cat("p = ", p, "\n")
  cat("n = ", n, "\n")
  cat("dimensión origdata: ", dim(origdata), "\n")
  cat("origclass = ", origclass, "\n")
  rm(class.table)
  rm(g)
  rm(class.name)
  rm(p)
  rm(n)
  cat("\n\n ------------- \n\n")
  
  
  
  if(PPmethod=="LDA"){
    idx <- LDAopt_MOD(origclass, origdata)
  }else{
    idx <- PDAopt_MOD(origclass, origdata, q,weight,lambda)
  }
  projdata = apply(origdata, 1, function(x) sum(x*idx$projbest) )
  cp <- split_entro(origclass, projdata)
  #pm <- mean(projdata)
  if ( cp == max(projdata) ) cp <- sort(projdata)[ length(projdata) - 1]
  
  
  
  cat("\n\n ---- Corrida findproj_MOD POST PDA_opt ---- \n\n")
  cat("dimensión de origdata: ", dim(origdata), "\n")
  class.table <- table(origclass)
  cat("class.table: ", class.table, "\n")
  g <- length(class.table)
  cat("g = ", g, "\n")
  class.name <- names(class.table)
  cat("class.name = ", class.name, "\n")
  p <- ncol(origdata)
  n <- nrow(origdata)
  # comentarios
  cat("p = ", p, "\n")
  cat("n = ", n, "\n")
  cat("dimensión origdata: ", dim(origdata), "\n")
  cat("origclass = ", origclass, "\n")
  rm(class.table)
  rm(g)
  rm(class.name)
  rm(p)
  rm(n)
  cat("\n\n ------------- \n\n")
  
  
  
  
  list(
    Index  = idx$indexbest, 
    Alpha = idx$projbest, 
    C = cp, 
    IOindexL = projdata <= cp,
    IOindexR = projdata > cp
  )
}