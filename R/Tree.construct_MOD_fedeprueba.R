#' Construct the projection pursuit classification tree MOD (NEW)
#' 
#' Find tree structure using various projection pursuit indices of 
#' classification in each split.
#' @title Projection pursuit classification tree MOD
#' @usage Tree.construct_MOD_fedeprueba(origclass,origdata,Tree.Struct, id,rep,rep1,rep2,
#' projbest.node,splitCutoff.node,PPmethod,r = NULL, 
#'lambda=NULL,TOL,maxiter=50000,q=1,weight=TRUE,tol = .5,strule ,tot,...) 
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
#' @param strule select the stoping rule rule based in G=1 pure node
#' @param tot total obs original class
#' @param ... something
#' @useDynLib PPtreeExt
#' @importFrom Rcpp evalCpp
#' @export

Tree.construct_MOD_fedeprueba <- 
  function(origclass,origdata,Tree.Struct, id,rep,rep1,rep2,projbest.node,splitCutoff.node,PPmethod,
           r = NULL, lambda=NULL,TOL,maxiter=50000,q=1,weight=TRUE,tol = .5,strule, tot,...) {

    class.table <- table(origclass)
    class.name <- names(class.table)
 
    rm(class.table)
    rm(class.name)
    
    flag = FALSE
    if(!is.matrix(origdata)) {
      flag <- TRUE
    }
    
    
    origclass <- as.integer(origclass)
    origdata <- as.matrix(origdata)
    n <- nrow(origdata)
    g <- table(origclass)
    G <- length(g)
    
    class.table <- table(origclass)

    class.name <- names(class.table)
   
    rm(class.table)
    rm(class.name)
 
    if(length(Tree.Struct) == 0){
      #Tree.Struct<-matrix(1:(2*G-1),ncol=1)
      # Tree.Struct<- cbind(Tree.Struct,0,0,0,0)
      Tree.Struct <- matrix(0, nrow=1, ncol=5)
    }
    if (id > nrow(Tree.Struct) ) {
      Tree.Struct <- rbind(Tree.Struct, matrix(0, nrow=id - nrow(Tree.Struct), ncol=5) )
    }
    ##To see the effect of diferent rules
    #end.node = (G==1 | length(origclass)/tot <= .5| entropy(origclass)<tol)
    end.node <- 0
    if(strule==1) {
        end.node <- 1*(G == 1)
    }else if(strule==2) {
       end.node <- 1*(length(origclass)/tot <= .05)
    }else{
     end.node <- 1*(entropy(origclass) < tol)
    }
    #end.node = (G==1 | length(origclass) <= 30 | entropy(origclass)<tol)
    #,pure=TRUE,nodesize=FALSE,entronode=FALSE,tot,
     cnd <- (end.node == 1) | (NROW(origdata) < 10) | !is.matrix(origdata) | flag
    
    
    if( cnd > 0 ){
    
      #if (nrow(origdata) > 0 ) {
      Tree.Struct[id,3] <- as.integer( names(g)[which.max(g)] )
      Tree.Struct[,1] <- 1:nrow(Tree.Struct)
      list(Tree.Struct=Tree.Struct,projbest.node=projbest.node, 
           splitCutoff.node=splitCutoff.node,rep=rep,rep1=rep1,rep2=rep2)
      # } else {
      #   list(Tree.Struct=Tree.Struct,projbest.node=projbest.node, 
      #        splitCutoff.node=splitCutoff.node,rep=rep,rep1=rep1,rep2=rep2)
      # }
    } else {
      
      Tree.Struct.row <- numeric(5)
      Tree.Struct.row[1] <- id
      Tree.Struct.row[2]<-rep1
      rep1<-rep1+1
      Tree.Struct.row[3]<-rep1
      rep1<-rep1+1
      Tree.Struct.row[4]<-rep2
      rep2<-rep2+1
      
      
      
      class.table <- table(origclass)
      class.name <- names(class.table)
      rm(class.table)
      rm(class.name)
 
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
      
    
      class.table <- table(t.class)
    
      class.name <- names(class.table)
     
      rm(class.table)
      rm(class.name)
 
      b<-Tree.construct_MOD(t.class,t.data,Tree.Struct, 
                            Tree.Struct[id, 2],rep,rep1,rep2,projbest.node, 
                            splitCutoff.node,PPmethod,r,lambda,TOL,maxiter,strule=strule,tot=tot,...)
      
    
      class.table <- table(t.class)
      class.name <- names(class.table)
      
      rm(class.table)
      rm(class.name)
      
      
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
      
      class.table <- table(t.class)
      class.name <- names(class.table)
      rm(class.table)
      rm(class.name)
     
    
      b<-Tree.construct_MOD(t.class,t.data,Tree.Struct, 
                            Tree.Struct[id,3],rep,rep1,rep2,projbest.node, 
                            splitCutoff.node,PPmethod,r,lambda,TOL,maxiter,strule=strule,tot=tot,...)
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
