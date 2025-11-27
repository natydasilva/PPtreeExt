#' Construct the projection pursuit classification tree with the extensions 
#' 
#' Find tree structure using various projection pursuit indices of 
#' classification in each split.
#' @title Projection pursuit classification tree MOD
#' @usage TreeExt.construct(origclass, origdata, Tree.Struct, id, rep, rep1, rep2,
#' projbest.node, splitCutoff.node, PPmethod, 
#' lambda = NULL, q = 1, weight = TRUE, srule=TRUE, tot=NULL, tol = .5,...) 
#' @param origclass original class 
#' @param origdata original data
#' @param Tree.Struct tree structure of projection pursuit classification tree
#' @param id tree node id
#' @param rep internal counter for nodes
#' @param rep1 internal counter for nodes
#' @param rep2 internal counter for nodes
#' @param projbest.node bests projection node
#' @param splitCutoff.node cutof node
#' @param PPmethod method for projection pursuit; "LDA", "PDA"
#' @param lambda lambda in PDA index
#' @param q numeric value with dimension of the projected data, if it is 1 then 1D projection is used
#' @param weight weight flag in LDA, PDA 
#' @param srule stopping rule flag; if TRUE use stopping rule, if FALSE stop only for pure or empty nodes
#' @param tot total number of observations
#' @param tol tolerance value for entropy stopping rule for splitting a node
#' @param ... additional arguments to pass trough
#' @useDynLib PPtreeExt
#' @importFrom Rcpp evalCpp
#' @export
TreeExt.construct <- function(origclass, origdata, Tree.Struct, id, rep, rep1,
                              rep2, projbest.node, splitCutoff.node,
                              PPmethod, lambda = NULL,
                               q = 1, weight = TRUE,
                              srule=TRUE, tot = NULL, tol = .5,...) {
  
  origclass <- as.integer(origclass)
  origdata <- as.matrix(origdata)
  n <- nrow(origdata)
  g <- table(origclass)
  G <- length(g)
  if (length(Tree.Struct) == 0) {
    Tree.Struct <- matrix(0, nrow = 1, ncol = 5)
  }
  if (id > nrow(Tree.Struct)) {
    Tree.Struct <- rbind(Tree.Struct, matrix(0, nrow = id - nrow(Tree.Struct), ncol =
                                               5))
  }
  

  
  if(srule){
   if( is.na(entropy(origclass))) ent <- 0 
   if (length(entropy(origclass)) == 0 ) ent <- 0
    
  cnd <- (G == 1) |
          (length(origclass)/tot <= 0.05) |
         (entropy(origclass) < tol)
  }else{
  cnd <- (G==1)
  }
  
  if (cnd) {
    Tree.Struct[id, 3] <- as.integer(names(g)[which.max(g)])
    Tree.Struct[, 1] <- 1:nrow(Tree.Struct)
    list(
      Tree.Struct = Tree.Struct,
      projbest.node = projbest.node,
      splitCutoff.node = splitCutoff.node,
      rep = rep,
      rep1 = rep1,
      rep2 = rep2
    )
    
  } else {
    Tree.Struct.row <- numeric(5)
    Tree.Struct.row[1] <- id
    Tree.Struct.row[2] <- rep1
    rep1 <- rep1 + 1
    Tree.Struct.row[3] <- rep1
    rep1 <- rep1 + 1
    Tree.Struct.row[4] <- rep2
    rep2 <- rep2 + 1
    a <- findproj_Ext(origclass, origdata, PPmethod, q = 1,  weight = TRUE, lambda)
    Tree.Struct.row[5] <- a$Index
    
    Tree.Struct[id, ] <- Tree.Struct.row
    
    splitCutoff.node <- rbind(splitCutoff.node, a$C)
    projbest.node <- rbind(projbest.node, matrix(a$Alpha, ncol = length(a$Alpha)))
    t.class <- origclass
    t.data <- origdata
    t.class <- t.class * a$IOindexL
    t.n <- length(t.class[t.class == 0])
    t.index <- sort.list(t.class)
    t.index <- sort(t.index[-(1:t.n)])
    t.class <- t.class[t.index]
    t.data <- origdata[t.index, ]
    
    #cat("Node size:", length(origclass), "  tot:", tot, "  G:", G, "\n")
    
    b <- TreeExt.construct(origclass= t.class, origdata = t.data, Tree.Struct =Tree.Struct, id = Tree.Struct[id, 2],
                           rep = rep, rep1=rep1, rep2=rep2, projbest.node =projbest.node, splitCutoff.node = splitCutoff.node,
                           PPmethod = PPmethod,  lambda=lambda,srule=srule, tot=tot, tol=tol, ...)
    
    Tree.Struct <- b$Tree.Struct
    projbest.node <- b$projbest.node
    splitCutoff.node <- b$splitCutoff.node
    rep <- b$rep
    rep1 <- b$rep1
    rep2 <- b$rep2
    t.class <- origclass
    t.data <- origdata
    t.class <- (t.class * a$IOindexR)
    t.n <- length(t.class[t.class == 0])
    t.index <- sort.list(t.class)
    t.index <- sort(t.index[-(1:t.n)])
    t.class <- t.class[t.index]
    t.data <- origdata[t.index, ]
    n <- nrow(t.data)
    G <- length(table(t.class))
    
    b <- TreeExt.construct(origclass=t.class, origdata=t.data, Tree.Struct=Tree.Struct, id =Tree.Struct[id, 3], rep=rep, rep1=rep1,
                           rep2=rep2, projbest.node=projbest.node, splitCutoff.node=splitCutoff.node, PPmethod = PPmethod, 
                           lambda =lambda,  srule=srule, tot=tot, tol=tol, ... )
    Tree.Struct <- b$Tree.Struct
    projbest.node <- b$projbest.node
    splitCutoff.node <- b$splitCutoff.node
    rep <- b$rep
    rep1 <- b$rep1
    rep2 <- b$rep2
  }
  list(Tree.Struct = Tree.Struct,
    projbest.node = projbest.node,
    splitCutoff.node = splitCutoff.node,
    rep = rep,
    rep1 = rep1,
    rep2 = rep2
  )
}