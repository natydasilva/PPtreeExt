#' Construct the projection pursuit classification tree MOD (NEW)
#' 
#' Find tree structure using various projection pursuit indices of 
#' classification in each split.
#' @title Projection pursuit classification tree MOD
#' @usage Tree.construct_MOD(origclass,origdata,Tree.Struct, id,rep,rep1,rep2,
#' projbest.node,splitCutoff.node,PPmethod,r = NULL, 
#'lambda=NULL,TOL,maxiter=50000,q=1,weight=TRUE,tol = .5,...) 
#' @param origclass original class 
#' @param origdata original data
#' @param Tree.Struct tree structure of projection pursuit classification tree
#' @param id something
#' @param rep something
#' @param rep1 something
#' @param rep2 something
#' @param projbest.node somenthing
#' @param splitCutoff.node something
#' @param PPmethod method for projection pursuit; "LDA", "PDA"
#' @param r r in Lr index
#' @param lambda something
#' @param TOL dfasd
#' @param maxiter something
#' @param q something
#' @param weight weight flag in LDA, PDA 
#' @param tol something
#' @param ... something
#' @useDynLib PPtreeExt
#' @importFrom Rcpp evalCpp
#' @export
#' @examples
#' data(iris)
#' Tree.result <- PPtreeViz::PPTreeclass(Species~.,data = iris,"LDA")
#' Tree.result
Tree.construct_MOD <- 
  function(origclass,origdata,Tree.Struct, id,rep,rep1,rep2,projbest.node,splitCutoff.node,PPmethod,
           r = NULL, lambda=NULL,TOL,maxiter=50000,q=1,weight=TRUE,tol = .5,...) {
    
    origclass <- as.integer(origclass)
    origdata <- as.matrix(origdata)
    n <- nrow(origdata)
    g <- table(origclass)
    G <- length(g)
    if(length(Tree.Struct) == 0){
      #Tree.Struct<-matrix(1:(2*G-1),ncol=1)
      # Tree.Struct<- cbind(Tree.Struct,0,0,0,0)
      Tree.Struct <- matrix(0, nrow=1, ncol=5)
    }
    if (id > nrow(Tree.Struct) ) {
      Tree.Struct <- rbind(Tree.Struct, matrix(0, nrow=id - nrow(Tree.Struct), ncol=5) )
    }
    
    end.node = (G==1 | length(origclass) <= 30 | entropy(origclass)<tol)
    if( end.node ){
      Tree.Struct[id,3] <- as.integer( names(g)[which.max(g)] )
      Tree.Struct[,1] <- 1:nrow(Tree.Struct)
      list(Tree.Struct=Tree.Struct,projbest.node=projbest.node, 
           splitCutoff.node=splitCutoff.node,rep=rep,rep1=rep1,rep2=rep2)
    } else {
      
      Tree.Struct.row <- numeric(5)
      Tree.Struct.row[1] <- id
      Tree.Struct.row[2]<-rep1
      rep1<-rep1+1
      Tree.Struct.row[3]<-rep1
      rep1<-rep1+1
      Tree.Struct.row[4]<-rep2
      rep2<-rep2+1
      a<-findproj_MOD(origclass,origdata,PPmethod,q=1,weight=TRUE,lambda)
      #a<-findproj_modLDA(origclass,origdata )
      Tree.Struct.row[5]<-a$Index
      
      
      Tree.Struct[id, ] <- Tree.Struct.row
      
      splitCutoff.node<-rbind(splitCutoff.node, a$C)
      projbest.node<-rbind(projbest.node,matrix(a$Alpha, ncol=length(a$Alpha)) )
      t.class<-origclass
      t.data<-origdata
      t.class<-t.class*a$IOindexL
      t.n<-length(t.class[t.class==0])
      t.index<-sort.list(t.class)
      t.index<-sort(t.index[-(1:t.n)])
      t.class<-t.class[t.index]
      t.data<-origdata[t.index,]
      
      #Tree.Struct<- Tree.Struct
      
      b<-Tree.construct_MOD(t.class,t.data,Tree.Struct, 
                            Tree.Struct[id, 2],rep,rep1,rep2,projbest.node, 
                            splitCutoff.node,PPmethod,r,lambda,TOL,maxiter,...)
      Tree.Struct<-b$Tree.Struct
      projbest.node<-b$projbest.node
      splitCutoff.node<-b$splitCutoff.node
      rep<-b$rep
      rep1<-b$rep1
      rep2<-b$rep2
      t.class<-origclass
      t.data<-origdata
      t.class<-(t.class*a$IOindexR)
      t.n<-length(t.class[t.class==0])
      t.index<-sort.list(t.class)
      t.index<-sort(t.index[-(1:t.n)])
      t.class<-t.class[t.index]
      t.data<-origdata[t.index,]
      n<-nrow(t.data)
      G<-length(table(t.class))
      b<-Tree.construct_MOD(t.class,t.data,Tree.Struct, 
                            Tree.Struct[id,3],rep,rep1,rep2,projbest.node, 
                            splitCutoff.node,PPmethod,r,lambda,TOL,maxiter,...)
      Tree.Struct<-b$Tree.Struct
      projbest.node<-b$projbest.node
      splitCutoff.node<-b$splitCutoff.node
      rep<-b$rep
      rep1<-b$rep1
      rep2<-b$rep2
    }
    list(Tree.Struct=Tree.Struct,projbest.node=projbest.node, 
         splitCutoff.node=splitCutoff.node,rep=rep,rep1=rep1,rep2=rep2)
  }