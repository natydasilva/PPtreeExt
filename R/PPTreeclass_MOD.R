#' Construct the projection pursuit classification tree MOD
#' 
#' Find tree structure using various projection pursuit indices of 
#' classification in each split.
#' @title Projection pursuit classification tree 
#' @usage PPTreeclass_MOD(formula,data, PPmethod="LDA",weight=TRUE,
#'                      r=1,lambda=0.1,energy=0,maxiter=50000,strule = 1,tot,...) 
#' @param formula an object of class "formula"
#' @param data data frame
#' @param PPmethod method for projection pursuit; "LDA", "PDA"
#' @param weight weight flag in LDA, PDA and Lr index
#' @param r r in Lr index
#' @param lambda lambda in PDA index
#' @param energy parameter for the probability to take new projection
#' @param maxiter maximum iteration number
#' @param strule select the stoping rule rule based in G=1 pure node
#' @param tot total obs original class
#' @param ... arguments to be passed to methods
#' @return Tree.Struct tree structure of projection pursuit classification tree
#' @return projbest.node 1 dimensional optimal projections of each node split
#' @return splitCutoff.node cutoff values of each node split 
#' @return origclass original class 
#' @return origdata original data
#' @references Lee, YD, Cook, D., Park JW, and Lee, EK(2013) 
#' PPtree: Projection Pursuit Classification Tree, 
#' Electronic Journal of Statistics, 7:1369-1386.
#' @export
#' @keywords tree
#' @examples
#' data(iris)
#' Tree.result <- PPTreeclass_MOD(Species~.,data = iris,"LDA")
#' Tree.result
PPTreeclass_MOD <- function(formula, data, PPmethod = "LDA", weight = TRUE, r = 1,
                         lambda = 0.1, energy = 0, maxiter = 50000, strule = 1, tot, ...){
  Call <- match.call()
  indx<-match(c("formula", "data"), names(Call), nomatch=0L)
  if(indx[1] == 0L) 
    stop("a 'formula' argument is required")
  temp<-Call[c(1L,indx)]
  temp[[1L]]<-quote(stats::model.frame)
  m<-eval.parent(temp)
  Terms<-attr(m,"terms")
  formula<-as.character(formula)
  class.n<-formula[2]
  data.n<-strsplit(formula[3]," \\+ ")[[1]]
  int.flag<-any(strsplit(formula[3]," \\* ")[[1]] == formula[3])
  if(data.n[1]=="."){
    tot.n<-class.n
  } else{
    tot.n<-c(class.n,data.n)
  }     
  if(!int.flag){
    stop("PPTreeclass cannot treat interaction terms")
  }else if(!sum(duplicated(c(colnames(data),tot.n))[-c(1:ncol(data))])==length(tot.n)){
    
  }else{
    origclass<-data[,class.n]
    if(data.n[1]=="."){
      origdata<-data[,colnames(data)!=class.n]
    }else {
      origdata<-data[,data.n,drop=FALSE]
    }
  }
  TOL<-NULL
  origdata<-as.matrix(origdata)
  
  class.table <- table(origclass)

  g <- length(class.table)
  
  class.name <- names(class.table)
 
  rm(class.table)
  rm(g)
  rm(class.name)
  
  
  splitCutoff.node<-NULL
  projbest.node<-NULL
  Tree.Struct<-NULL
  id<-1
  rep1<-2
  rep2<-1
  rep<-1
  Tree.final <- Tree.construct_MOD(origclass, origdata, Tree.Struct, id, rep, rep1, 
                                      rep2, projbest.node, splitCutoff.node,
                                      PPmethod, r, lambda, TOL, maxiter,strule = strule, tot=tot, ...)                            
  Tree.Struct <- Tree.final$Tree.Struct
  colnames(Tree.Struct)<-c("id","L.node.ID","R.F.node.ID","Coef.ID","Index")
  projbest.node<-Tree.final$projbest.node
  splitCutoff.node<-Tree.final$splitCutoff.node
  #colnames(splitCutoff.node)<-paste("Rule",1,sep="")
  treeobj<-list(Tree.Struct=Tree.Struct,projbest.node=projbest.node, 
                splitCutoff.node=splitCutoff.node,origclass=origclass,
                origdata=origdata,terms=Terms)
  class(treeobj)<-append(class(treeobj),"PPtreeclassMOD")
  return(treeobj)
}
