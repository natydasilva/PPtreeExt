#' Shiny app to explore PPtree algirith partitions
#' 
#' @usage explorapp(ui, server) 
#' @param ui A 
#' @param server aaa
#' @useDynLib PPtreeExt
#' @importFrom Rcpp evalCpp
#' @export
#' @examples
#' \dontrun{
#' explorapp(ui,server)
#' }
explorapp <-function(ui,server){
  
  X1 <- NULL
  X2<- NULL
  ppred <- NULL
  Sim <- NULL
  pred <- NULL
  predict <- NULL
  #column <- NULL
simu3 <- function(mux1, mux2, muy1, muy2, muz1, muz2,  cor1, cor2, cor3, n1 = 100, n2 = 100, n3 = 100) {
  set.seed(666)
  bivn <- MASS::mvrnorm(n1, mu = c(mux1, mux2), Sigma = matrix(c(1, cor1, cor1, 1), 2))
  bivn2 <- MASS::mvrnorm(n2, mu = c(muy1, muy2), Sigma = matrix(c(1, cor2, cor2, 1), 2))
  bivn3 <- MASS::mvrnorm(n3, mu = c(muz1, muz2), Sigma = matrix(c(1, cor3, cor3, 1), 2))
  
  d1 <- data.frame(Sim = "sim1", bivn)
  d2 <- data.frame(Sim = "sim2", bivn2)
  d3 <- data.frame(Sim = "sim3", bivn3)
  rbind(d1, d2, d3)
}




ppbound <- function(ru, data , meth, entro ,title,simM=FALSE){

  grilla <- base::expand.grid(X1 = seq((min(data$X1) + sign( min(data$X1))*.5) , (max(data$X1) + sign(max(data$X1))*.5), , 100),
                              X2 = seq((min(data$X2 ) + sign( min(data$X2))*.5), (max(data$X2) + sign(max(data$X2))*.5), , 100))

   if(meth == "Original"){
    pptree <- PPtreeViz::PPTreeclass(Sim~., data = data, "LDA")
    ppred.sim <- PPtreeViz::PPclassify(pptree, test.data = grilla, Rule = ru)
    grilla$pred <- ppred.sim[[2]]
    err <- round(PPtreeViz::PPclassify(pptree, test.data=data[,-1], true.class = data[,1], Rule=ru)[[1]]/nrow(data[,-1]),3)*100
  }
  if(meth == "Rpart"){
    rpart.mod<- rpart::rpart(Sim ~., data = data)
    grilla$pred <- predict(rpart.mod, newdata = grilla, type = "class")
    err <- round(1-sum(diag(table(predict(rpart.mod, newdata=data[,-1],type = "class") ,data[,1])))/nrow(data[,-1]),3)*100
    
    }
  
  if(entro) {mod = 2
  }else{
    mod = 1
  }
  
  if(meth == "Modified"){
    pptree <- PPtree_splitMOD(Sim~., data = data, "LDA", entro = entro)
    ppred.sim <- PPtreeViz::PPclassify(pptree, test.data = grilla, Rule = ru)
    grilla$pred <- paste("sim", ppred.sim[[2]], sep = "")
    err <- round(PPtreeViz::PPclassify(pptree, test.data=data[,-1], true.class = data[,1], Rule=ru)[[1]]/nrow(data[,-1]),3)*100
    
  }
  
  #ruleid <- pptree$splitCutoff.node[,ru]
  if(simM){
    pl.pp  <- ggplot2::ggplot(data = grilla ) + ggplot2::geom_point( ggplot2::aes(x = X1, y = X2, color = as.factor(pred) ), alpha = .20) +
      ggplot2::scale_colour_brewer(name = "Class",type = "qual", palette = "Dark2" ) + ggplot2::theme_bw() +
     ggplot2::geom_point(data = data, ggplot2::aes(x = X1 , y = X2, group = Sim,  color =Sim), size = I(3)  ) + 
      ggplot2::theme(legend.position = "none",aspect.ratio = 1) + 
      ggplot2::scale_y_continuous(expand = c(0, 0) ) + ggplot2::scale_x_continuous(expand = c(0, 0)) + 
      ggplot2::labs(x = " ", y = "", title = paste(title, "error", err,"%" ))
  }else{
    pl.pp  <- ggplot2::ggplot(data = grilla ) + ggplot2::geom_point( ggplot2::aes(x = X1, y = X2, color = as.factor(pred), shape = as.factor(pred) ), alpha = .20) +
       ggplot2::scale_colour_brewer(name = "Class",type = "qual", palette = "Dark2" ) + ggplot2::theme_bw() +
    ggplot2::scale_shape_discrete(name='Class') + ggplot2::geom_point(data = data, ggplot2::aes(x = X1 , y = X2, group = Sim, shape = Sim, color = Sim), size = I(3)  ) + 
    ggplot2::theme(legend.position = "none",aspect.ratio = 1) + 
    ggplot2::scale_y_continuous(expand = c(0, 0) ) + ggplot2::scale_x_continuous(expand = c(0, 0)) + 
    ggplot2::labs(x = " ", y = "", title = paste(title, "error", err,"%" ))
  }

  pl.pp
}


ppboundMOD <- function( data , meth = "MOD", entro = FALSE, entroindiv = TRUE, title, simM=FALSE ){
  
   grilla <- base::expand.grid(X1 = seq((min(data$X1) + sign( min(data$X1))*.5) , (max(data$X1) + sign(max(data$X1))*.5), , 100),
                              X2 = seq((min(data$X2 ) + sign( min(data$X2))*.5), (max(data$X2) + sign(max(data$X2))*.5), , 100))

  
  
  pptree <- PPTreeclass_MOD(Sim~. ,  data = data, PPmethod = 'LDA')
  
  ppred.sim <- PPclassify_MOD(pptree, test.data = grilla)
  grilla$ppred <-ppred.sim[[2]]
  err <- round(PPclassify_MOD(pptree, test.data=data[,-1], true.class = data[,1])[[1]]/nrow(data[,-1]),3)*100
  
if(simM){
  pl.pp <- ggplot2::ggplot(data = grilla ) + ggplot2::geom_point( ggplot2::aes(x = X1, y = X2, color = ppred), alpha = .20) +
  ggplot2::scale_colour_brewer(name = "Class",type = "qual", palette = "Dark2" ) + ggplot2::theme_bw() +
 ggplot2::scale_shape_discrete(name='Class') + ggplot2::geom_point(data = data, ggplot2::aes(x = X1 , y = X2, group = Sim, color = Sim), size = I(3)  ) + 
    ggplot2::theme(legend.position = "none",aspect.ratio = 1) + 
    ggplot2::scale_y_continuous(expand = c(0, 0) ) + ggplot2::scale_x_continuous(expand = c(0, 0)) + 
    ggplot2::labs(x = " ", y = "", title = paste(title,"error",err,"%" ))
  
}else{
  pl.pp <- 
    ggplot2::ggplot(data = grilla ) + ggplot2::geom_point( ggplot2::aes(x = X1, y = X2, color = ppred, shape = ppred ), alpha = .20) +
    ggplot2::scale_colour_brewer(name = "Class",type = "qual", palette = "Dark2" ) + ggplot2::theme_bw() +
    ggplot2::scale_shape_discrete(name='Class') + ggplot2::geom_point(data = data, ggplot2::aes(x = X1 , y = X2, group = Sim, shape = Sim, color = Sim), size = I(3)  ) + 
    ggplot2::theme(legend.position = "none",aspect.ratio = 1) + 
    ggplot2::scale_y_continuous(expand = c(0, 0) ) + ggplot2::scale_x_continuous(expand = c(0, 0)) + 
    ggplot2::labs(x = " ", y = "", title = paste(title,"error",err,"%" ))
  
}
  pl.pp
}

ui <- shiny::fluidPage(
    shiny::mainPanel(
    shiny::tabsetPanel(
      shiny::tabPanel(
        "SIM", 
        shiny::fluidRow(shiny::column(3, shiny::selectInput(inputId = "rule", label ="Rule", choices = 1:8, selected = 1 ) ), shiny::column(3, shiny::selectInput(inputId = "modi", label ="Modification", choices = 1:3, selected = 1 )  )),
        shiny::fluidRow( shiny::column(4,
      shiny::textInput( inputId = 'mean', label = 'Group means ', value = 
                   "-1, 0.6, 0, -0.6, 2,-1" ) ),
      shiny::column(4, shiny::textInput(inputId = "cor", label = "Correlations",
                value = "0.95, 0.95, 0.95")) ,
      shiny::column(4, shiny::textInput(inputId = "sample", label = "Group sample",
        value = "100, 100, 100") )), shiny::fluidRow(shiny::actionButton("do", label = "OK")), 
    shiny::fluidRow(
      shiny::plotOutput("distPlot" ) ) ),
    shiny::tabPanel("SIM outliers",
                    shiny::fluidRow(shiny::column(4, shiny::selectInput(inputId = "rule2",label ="Rule", choices = 1:8, selected = 1 ) ), shiny::column(3, shiny::selectInput(inputId = "modi2", label ="Modification", choices = 1:3, selected = 1 ) )),
                    shiny::fluidRow( shiny::column(4, shiny::textInput( inputId = 'mean2', label = 'Group means ', value = 
                                           "-1, 0.6, 0, -0.6, 2,-1" )),
                      shiny:: column(4, shiny::textInput(inputId = "cor2", label = "Correlations",
                                           value = "0.95, 0.95, 0.95")) ,
                       shiny::column(4, shiny::textInput(inputId = "sample2", label = "Group sample",
                                           value = "100, 100, 100") )),
                    shiny::fluidRow(shiny::column(4, shiny::selectInput(inputId = "group",label ="Add outliers to class", choices = 1:3, selected = 1 ))),
                    shiny::fluidRow(shiny::column(4, shiny::textInput( inputId = 'meanout', label = 'Out. X1, X2 means ', value = 
                                             "-3, 3" )),
                    shiny::column(4, shiny::textInput( inputId = 'sdout', label = 'Out. X1, X2 sd ', value =".2,.2" 
                                             )),
                    shiny::column(4, shiny::textInput(inputId = "sampleout", label = "Out. sample size",
                                          value = "10") ),
                      shiny::fluidRow(shiny::actionButton("do2", label = "OK"))), 
             
                    shiny::fluidRow(
                      shiny::plotOutput("distPlot2") )
             ), 
    ##
    shiny::tabPanel(
      "MixSim",
      shiny::fluidRow(shiny::column(4,
          shiny::selectInput(inputId = "rule3", label ="Rule", choices = 1:8, selected = 1 )),
          shiny::column(4,shiny::selectInput(inputId = "modi3", label ="Modification", choices = 1:3, selected = 1 ))),
      shiny::fluidRow(shiny::column(4, shiny::numericInput("size", label = "Sample size", value = 100)),
            shiny::column(4, shiny::numericInput("BarOmega", label = "BarOmega desired average overlap", value = 0.05))),
            shiny::fluidRow(shiny::column(4,shiny::numericInput("MaxOmega", label = "MaxOmega desired maximum overlap", value = 0.15)),
      shiny::column(4,shiny::numericInput("K", label = "K number of components", value = 3))),
          #shiny::numericInput("p", label = "number of dimensions", value = 5),
          shiny::fluidRow(shiny::actionButton("simmaitra","simula")),
        shiny::fluidRow(
        
          shiny::plotOutput("plotsmaitra")
        )
    ##
    
    )
  )))


server <- function(input, output) {

  output$distPlot <- shiny::renderPlot({
    if(input$do){
   
    x1 <- shiny::isolate(as.numeric(unlist(strsplit(input$mean, ",") ) ))
    x2 <- shiny::isolate(as.numeric(unlist(strsplit(input$cor, ",") ) ))
    x3 <- shiny::isolate(as.numeric(unlist(strsplit(input$sample, ",") ) ))
    
    dat.pl2 <-  shiny::isolate(simu3(x1[1], x1[2], x1[3],x1[4],x1[5],x1[6],
                     x2[1], x2[2], x2[3], x3[1], x3[2], x3[3]))
    
   
    if(  input$modi==1){
     modpl <- ppbound(ru = as.numeric(input$rule),  data = dat.pl2, meth = "Modified" , entro = FALSE, title = "Modified subset")
    }
    if(  input$modi==2){
      modpl <- ppbound(ru =  as.numeric(input$rule),  data = dat.pl2, meth = "Modified" , entro = TRUE, title = "Modified entropy mp groups")
    }
    if(  input$modi==3){
      modpl <-  ppboundMOD(data = dat.pl2, meth = "MOD", entro = FALSE, entroindiv = TRUE, title = "Modified entropy ind partitions")
    }
    
     gridExtra::grid.arrange(ppbound(ru =  as.numeric(input$rule),data = dat.pl2, meth = "Rpart", entro = TRUE ,title ="Rpart"),
                 ppbound(ru =  as.numeric(input$rule), data = dat.pl2, meth = "Original" , entro = FALSE, title ="Original"), 
                 #ppbound(ru =  as.numeric(input$rule),  data = dat.pl2, meth = "Modified" , entro = TRUE),
                modpl, ncol = 3)
    }
  })

  
  

    output$distPlot2 <- shiny::renderPlot({
      if(input$do2){
      x1 <- shiny::isolate(as.numeric(unlist(strsplit(input$mean2, ",") ) ) )
      x2 <- shiny::isolate(as.numeric(unlist(strsplit(input$cor2, ",") ) ) )
      x3 <- shiny::isolate(as.numeric(unlist(strsplit(input$sample2, ",") ) ) )
      x4 <- shiny::isolate(as.numeric(unlist(strsplit(input$meanout, ",") ) ) )
      x5 <- shiny::isolate(as.numeric(unlist(strsplit(input$sdout, ",") ) ) )
      x6 <- shiny::isolate(as.numeric(unlist(strsplit(input$sampleout, ",") ) ) )
      dat.pl2 <- simu3(x1[1], x1[2], x1[3],x1[4],x1[5],x1[6],
                       x2[1], x2[2], x2[3], x3[1], x3[2], x3[3])
      
      
       set.seed(123)
       aux <- data.frame(Sim = rep(paste("sim",as.numeric(input$group),sep=""), x6), X1 = stats::rnorm( n = x6, mean = x4[1], sd = x5[1]), X2 = stats::rnorm( n = x6, mean = x4[2], sd =x5[2]))
      dat.pl2 <- rbind(dat.pl2, aux) 
      
      if(  input$modi2==1){
        modpl <- ppbound(ru =  as.numeric(input$rule),  data = dat.pl2, meth = "Modified" , entro = FALSE, title="Modified subset")
      }
      if(  input$modi2==2){
        modpl <- ppbound(ru =  as.numeric(input$rule),  data = dat.pl2, meth = "Modified" , entro = TRUE, title="Modified entropy mp groups")
      }
      if(  input$modi2==3){
        modpl <-  ppboundMOD(data = dat.pl2, meth = "MOD", entro = FALSE, entroindiv = TRUE, title = "Modified entropy ind partitions")
      }
      
      gridExtra::grid.arrange( ppbound(ru =  as.numeric(input$rule2),  data = dat.pl2, meth = "Rpart" , entro = FALSE, title ="Rpart"),
                               ppbound(ru =  as.numeric(input$rule2), data = dat.pl2, meth = "Original", entro = TRUE, title= "Original"),
                   
                  # ppbound(ru =  as.numeric(input$rule2), FALSE, data = dat.pl2, meth = "Modified" , entro = TRUE),
                  modpl, ncol = 3)
    
      
      }
      })
  
  
    
    output$plotsmaitra <- shiny::renderPlot({
      if(input$simmaitra){
      
        #Q <- MixSim(BarOmega = 0.01, K = 4, p = 2)
        set.seed(123)
        repeat{
        Q <- MixSim::MixSim(BarOmega = shiny::isolate(as.numeric(input$BarOmega)), MaxOmega = shiny::isolate(as.numeric(input$MaxOmega)), K = shiny::isolate(as.numeric(input$K)),
                    p=2, sph = FALSE, hom = TRUE,
               ecc = 0.90, PiLow = 1.0, int = c(0.0, 1.0), resN = 100,
               eps = 1e-06, lim = 1e06)
        if (Q$fail == 0) break
        }
        A <- MixSim::simdataset(n = shiny::isolate(as.numeric(input$size)), Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
        dat.pl2 <- data.frame(Sim=paste("sim",A[[2]],sep=""),X1=scale(A[[1]][,1]), X2=scale(A[[1]][,2]))
        
        
        if(  input$modi3==1){
          modpl <- ppbound(ru =  as.numeric(input$rule3),  data = dat.pl2, meth = "Modified" , entro = FALSE, title = "Modified subset", simM = TRUE )
        }
        if( input$modi3==2){
          modpl <- ppbound(ru =  as.numeric(input$rule3),  data = dat.pl2, meth = "Modified" , entro = TRUE, title="Modified entropy mp groups", simM = TRUE)
        }
        if(input$modi3 == 3){
          modpl <-  ppboundMOD(data = dat.pl2, meth = "MOD", entro = FALSE, entroindiv = TRUE, title = "Modified entropy ind partitions", simM = TRUE)
        }
        
        gridExtra::grid.arrange(ppbound(ru =  as.numeric(input$rule3),data = dat.pl2, meth = "Rpart", entro = TRUE ,title ="Rpart",simM=TRUE),
                                ppbound(ru =  as.numeric(input$rule3), data = dat.pl2, meth = "Original" , entro = FALSE, title ="Original", simM=TRUE), 
                                #ppbound(ru =  as.numeric(input$rule),  data = dat.pl2, meth = "Modified" , entro = TRUE),
                                modpl, ncol = 3)
      }
    })  
  
   
}

shiny::shinyApp(ui = ui, server = server)
}




