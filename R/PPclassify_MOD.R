#' Predict class for the test set with the fitted projection pursuit 
#' classification tree and calculate prediction error.
#' @title predict PPtree
#' @usage PPclassify_MOD(Tree.result,test.data,true.class=NULL,...)  
#' @param Tree.result PPtreeclass object 
#' @param test.data  the test dataset
#' @param true.class true class of test dataset if available
#' @param ... arguments to be passed to methods
#' @return predict.class predicted class
#' @return predict.error number of the prediction errors
#' @references Lee, YD, Cook, D., Park JW, and Lee, EK(2013) 
#' PPtree: Projection Pursuit Classification Tree, 
#' Electronic Journal of Statistics, 7:1369-1386.
#' @export
#' @keywords tree
#' @examples
#' data(penguins)
#' penguins <- na.omit(penguins[, -c(2,7, 8)])
#' require(rsample)
#' penguins_spl <- rsample::initial_split(penguins, strata=species)
#' penguins_train <- training(penguins_spl)
#' penguins_test <- testing(penguins_spl)
#' penguins_ppt <- PPTreeclass_MOD(species~bill_len + bill_dep +
#'   flipper_len + body_mass, data = penguins_train, PPmethod = "LDA")
#' PPclassify_MOD(penguins_ppt, test.data=penguins_test[,-1], true.class = penguins_test$species)

PPclassify_MOD <- function(Tree.result, test.data = NULL, true.class = NULL,...) {
  
  if(is.null(test.data))
    test.data<-Tree.result$origdata
  test.data<-as.matrix(test.data)
 
   if(!is.null(true.class)){  
    true.class<-as.matrix(true.class); 
    if(nrow(true.class)==1) 
      true.class<-t(true.class)
    if(!is.numeric(true.class)) {
      class.name<-names(table(true.class))
      temp<-rep(0,nrow(true.class))
      for(i in 1:length(class.name))
        temp<-temp+(true.class==class.name[i])*i
      true.class<-temp
    }
  }   
  
  PP.Classification<-function(Tree.Struct,test.class.index,IOindex,
                              test.class,id,rep){
    if(Tree.Struct[id,4]==0){
      i.class<-test.class
      i.class[i.class>0]<-1
      i.class<-1-i.class
      test.class<-test.class+IOindex*i.class*Tree.Struct[id, 3]
      return(list(test.class=test.class,rep=rep))
    } else{  
      IOindexL<-IOindex*test.class.index[rep,]
      IOindexR<-IOindex*(1-test.class.index[rep,])
      rep<-rep+1
      a<-PP.Classification(Tree.Struct,test.class.index,IOindexL,
                           test.class,Tree.Struct[id,2],rep)
      test.class<-a$test.class
      rep<-a$rep;
      a<-PP.Classification(Tree.Struct,test.class.index,IOindexR,
                           test.class,Tree.Struct[id,3],rep)
      test.class<-a$test.class
      rep<-a$rep
    }
    list(test.class=test.class,rep=rep)
  }
  
  PP.Class.index<-function(class.temp,test.class.index,test.data,
                           Tree.Struct,Alpha.Keep,C.Keep,id){
    class.temp<-as.integer(class.temp)
    if(Tree.Struct[id,2]==0){
      return(list(test.class.index=test.class.index,class.temp=class.temp))
    } else{
      t.class<-class.temp 
      t.n<-length(t.class[t.class==0])
      t.index<-sort.list(t.class)
      if(t.n)
        t.index<-sort(t.index[-(1:t.n)])
      t.data<-test.data[t.index,]
      id.proj<-Tree.Struct[id,4]
      
      proj.test<-as.matrix(test.data)%*%as.matrix(Alpha.Keep[id.proj,])
      proj.test<-as.double(proj.test)
      class.temp<-t(proj.test<C.Keep[id.proj]) 
      test.class.index<-rbind(test.class.index,class.temp)
      a<-PP.Class.index(class.temp,test.class.index,test.data,
                        Tree.Struct,Alpha.Keep,C.Keep,
                        Tree.Struct[id,2])
      test.class.index<-a$test.class.index
      a<-PP.Class.index(1-class.temp,test.class.index,test.data,
                        Tree.Struct,Alpha.Keep,C.Keep,
                        Tree.Struct[id,3])
      test.class.index<-a$test.class.index;
    }
    list(test.class.index=test.class.index,class.temp=class.temp)
  }
  
  n<-nrow(test.data)
  class.temp<-rep(1,n)
  test.class.index<-NULL
  temp <- PP.Class.index(class.temp,test.class.index,test.data,
                       Tree.result$Tree.Struct,Tree.result$projbest.node,
                       Tree.result$splitCutoff.node,1)
  test.class<-rep(0,n)
  IOindex<-rep(1,n)
  temp<-PP.Classification(Tree.result$Tree.Struct,temp$test.class.index,
                          IOindex,test.class,1,1)
  if(!is.null(true.class)){
    predict.error<-sum(true.class!=temp$test.class)
  } else {
    predict.error<-NA
  }  
  class.name<-names(table(Tree.result$origclass))
  predict.class <- class.name[temp$test.class]
  list(predict.error= predict.error, predict.class=predict.class)

}