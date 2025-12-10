#' Print Method for PPtreeExtclass Objects
#' 
#' Prints a summary of a fitted projection pursuit classification tree, including 
#' the tree structure, optionally the projection coefficients and cutoff values, 
#' and the training error rate.
#' 
#' @param x An object of class \code{"PPtreeExtclass"} from 
#'   \code{\link{PPtreeExtclass}} or \code{\link{PPtreeExt_split}}.
#' @param coef.print Logical indicating whether to print the projection coefficients 
#'   for each split node. Default is \code{FALSE}.
#' @param cutoff.print Logical indicating whether to print the cutoff values for 
#'   each split node. Default is \code{FALSE}.
#' @param verbose Logical indicating whether to print the tree structure and error 
#'   rate. If \code{FALSE}, the function returns the tree structure invisibly without 
#'   printing. Default is \code{TRUE}.
#' @param ... Additional arguments (currently not used).
#' @details
#' The function traverses the tree structure stored in \code{x$Tree.Struct} and 
#' creates a hierarchical text representation. When \code{coef.print = TRUE}, 
#' the projection coefficients (linear combinations of features) used at each 
#' split are displayed. When \code{cutoff.print = TRUE}, the threshold values 
#' used to determine left/right splits are shown.
#' 
#' The training error rate is computed by applying the fitted tree to the 
#' original training data.
#' 
#' @seealso \code{\link{PPtreeExtclass}}, \code{\link{PPtreeExt_split}}, 
#'   \code{\link{predict.PPtreeExtclass}}
#' @export
print.PPtreeExtclass <- function(x, coef.print  = FALSE, cutoff.print = FALSE,
                               verbose = TRUE,...){
  PPtreeOBJ <- x
  TS <- PPtreeOBJ$Tree.Struct
  Alpha <- PPtreeOBJ$projbest.node
  cut.off <- PPtreeOBJ$splitCutoff.node
  gName <- names(table(PPtreeOBJ$origclass))
  pastemake <- function(k,arg,sep.arg=""){
    temp <- ""
    for(i in 1:k)
      temp <- paste(temp, arg, sep = sep.arg)
    return(temp)
  }
  TreePrint <- "1) root"
  i <- 1  
  flag.L <- rep(FALSE,nrow(TS))
  keep.track <- 1
  depth.track <- 0
  depth<-0
  while(sum(flag.L)!=nrow(TS)){
    if(!flag.L[i]){                    
      if(TS[i,2] == 0) {
        flag.L[i] <-TRUE
        n.temp <-length(TreePrint)
        tempp <-strsplit(TreePrint[n.temp],") ")[[1]]
        temp.L <-paste(tempp[1],")*",tempp[2],sep="")
        temp.L <- paste(temp.L,"  ->  ","\"",gName[TS[i,3]],"\"",sep="")
        TreePrint <- TreePrint[-n.temp]
        id.l <-length(keep.track)-1
        i <- keep.track[id.l]
        depth <-depth -1
      } else if(!flag.L[TS[i,2]]){
        depth <- depth+1
        emptyspace <- pastemake(depth, "   ")
        temp.L <- paste(emptyspace, TS[i,2], ")  proj",
                      TS[i,4],"*X < cut",TS[i,4], sep = "")
        i <- TS[TS[i,2],1]   
      } else{
        depth <- depth +1
        emptyspace <- pastemake(depth,"   ")          
        temp.L <- paste(emptyspace,TS[i,3],")  proj",
                       TS[i,4],"*X >= cut",TS[i,4],sep="")
        flag.L[i] <- TRUE
        i <- TS[TS[i,3],1]
      } 
      keep.track <- c(keep.track,i)
      depth.track <- c(depth.track,depth)
      TreePrint <-c(TreePrint,temp.L)
    } else{
      id.l<-id.l-1
      i<-keep.track[id.l]
      depth<-depth.track[id.l]
    }
  }
  colnames(Alpha)<-colnames(PPtreeOBJ$origdata)
  rownames(Alpha)<-paste("proj",1:nrow(Alpha),sep="")  
  TreePrint.output<-
    paste("=============================================================",                          
          "\nProjection Pursuit Classification Tree Extension result",                           
          "\n=============================================================\n")
  for(i in 1:length(TreePrint))
    TreePrint.output<-paste(TreePrint.output,TreePrint[i],sep="\n")
  TreePrint.output<-paste(TreePrint.output,"\n",sep="")
  sample.data.X<-PPtreeOBJ$origdata
  sample.data.class<-PPtreeOBJ$origclass
  #error.rate<- PPtreeViz::PPclassify(PPtreeOBJ,sample.data.X)$predict.error
  error.rate.aux <- predict.PPtreeExtclass(object = PPtreeOBJ, newdata = sample.data.X, 
                                           true.class = sample.data.class)
  error.rate <- error.rate.aux$predict.error
  colnames(Alpha)<-paste(1:ncol(Alpha),":\"",colnames(Alpha),"\"",sep="")
  if(verbose){
    cat(TreePrint.output)
    if(coef.print){
      cat("\nProjection Coefficient in each node",
          "\n-------------------------------------------------------------\n")
      print(round(Alpha,4))
    }
    if(cutoff.print){
      cat("\nCutoff values of each node",
          "\n-------------------------------------------------------------\n")
      print(round(cut.off,4))
    }    
    cat("\nError rates",
        "\n-------------------------------------------------------------\n")
    print(round(error.rate,4))
    
  }    
  return(invisible(TreePrint)) 
}

