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
  #column <- NULL
simu3 <- function(mux1, mux2, muy1, muy2, muz1, muz2,  cor1,cor2,cor3, n1 = 100, n2 = 100, n3 = 100) {
  set.seed(666)
  bivn <- MASS::mvrnorm(n1, mu = c(mux1, mux2), Sigma = matrix(c(1, cor1, cor1, 1), 2))
  bivn2 <- MASS::mvrnorm(n2, mu = c(muy1, muy2), Sigma = matrix(c(1, cor2, cor2, 1), 2))
  bivn3 <- MASS::mvrnorm(n3, mu = c(muz1, muz2), Sigma = matrix(c(1, cor3, cor3, 1), 2))
  
  d1 <- data.frame(Sim = "sim1", bivn)
  d2 <- data.frame(Sim = "sim2", bivn2)
  d3 <- data.frame(Sim = "sim3", bivn3)
  rbind(d1, d2, d3)
}


ppbound <- function(ru, leg = TRUE, data , meth, entro ){
  grilla <- base::expand.grid(X1 = seq((min(data$X1) + sign( min(data$X1))*.5) , (max(data$X1) + sign(max(data$X1))*.5), , 100),
                              X2 = seq((min(data$X2 ) + sign( min(data$X2))*.5), (max(data$X2) + sign(max(data$X2))*.5), , 100))
  if(meth == "Original"){
    pptree <- PPtreeViz::PPTreeclass(Sim~., data = data, "LDA")
    ppred.sim <- PPtreeViz::PPclassify(pptree, test.data = grilla, Rule = ru)
    grilla$ppred <- ppred.sim[[2]]
  }
  
  if(entro) {mod = 2
  }else{
    mod = 1
  }
  
  if(meth == "Modified"){
    pptree <- PPtree_splitMOD(Sim~., data = data, "LDA", entro = entro)
    ppred.sim <- PPtreeViz::PPclassify(pptree, test.data = grilla, Rule = ru)
    grilla$ppred <- paste("sim", ppred.sim[[2]], sep = "")
    
  }
  
  ruleid <- pptree$splitCutoff.node[,ru]
  
  p <- ggplot2::ggplot(data = grilla ) + ggplot2::geom_point( ggplot2::aes(x = X1, y = X2, color = as.factor(ppred), shape = as.factor(ppred) ), alpha = .20) +
    ggplot2::geom_abline(intercept = ruleid[[1]] / pptree$projbest.node[[3]], slope = -pptree$projbest.node[[1]]/pptree$projbest.node[[3]], size = 1 )+ ggplot2::scale_colour_brewer(name = "Class",type = "qual", palette = "Dark2" ) + ggplot2::theme_bw() +
    ggplot2::geom_abline(intercept = ruleid[[2]] / pptree$projbest.node[[4]], slope = -pptree$projbest.node[[2]] / pptree$projbest.node[[4]], size = 1) + ggplot2::scale_shape_discrete(name='Class')
  
  if(leg){
    pl.pp <- p + ggplot2::geom_point(data = data, ggplot2::aes(x = X1 , y = X2, group = Sim, shape = Sim, color = Sim), size = I(3)  ) + 
      ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 6)) + 
      ggplot2::scale_y_continuous(expand = c(0, 0) ) + ggplot2::scale_x_continuous(expand = c(0, 0) ) + 
      ggplot2::labs(title = paste("Rule", ru, paste(meth, mod, sep = "" ) ) )
  }else{
    pl.pp <- p + ggplot2::geom_point(data = data, ggplot2::aes(x = X1 , y = X2, group = Sim, shape = Sim, color = Sim), size = I(3)  ) + 
      ggplot2::theme(legend.position = "none",aspect.ratio = 1) + 
      ggplot2::scale_y_continuous(expand = c(0, 0) ) + ggplot2::scale_x_continuous(expand = c(0, 0)) + 
      ggplot2::labs(x = " ", y = "", title = paste("Rule", ru, paste(meth, mod, sep = "") ) )
  }
  
  pl.pp
}


ui <- shiny::fluidPage(

    shiny::mainPanel(
    shiny::tabsetPanel(
      shiny::tabPanel(
        "SIM 1", 
        shiny::fluidRow(shiny::column(3, shiny::selectInput(inputId = "rule", label ="Rule", choices = 1:8, selected = 1 ) ) ),
        shiny::fluidRow( shiny::column(4,
      shiny::textInput( inputId = 'mean', label = 'Group means ', value = 
                   "-1, 0.6, 0, -0.6, 2,-1" )),
      shiny::column(4, shiny::textInput(inputId = "cor", label = "Correlations",
                value = "0.95, 0.95, 0.95")) ,
      shiny::column(4, shiny::textInput(inputId = "sample", label = "Group sample",
        value = "100, 100, 100") ), shiny::fluidRow(shiny::actionButton("do", label = "OK"))), 
     
   
    shiny::fluidRow(
      shiny::plotOutput("distPlot", width = "150%", height = "400px" ) ) ),
    shiny::tabPanel("SIM 2",
                    shiny::fluidRow(shiny::column(4, shiny::selectInput(inputId = "rule2",label ="Rule", choices = 1:8, selected = 1 ) ) ),
                    shiny::fluidRow( shiny::column(4, shiny::textInput( inputId = 'mean2', label = 'Group means ', value = 
                                           "-1, 0.6, 0, -0.6, 2,-1" )),
                      shiny::column(4, shiny::textInput(inputId = "cor2", label = "Correlations",
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
                      shiny::plotOutput("distPlot2", width = "150%", height = "400px" ) )
             )
  )))


server <- function(input, output) {

  output$distPlot <- shiny::renderPlot({
    input$do
   
    x1 <- shiny::isolate(as.numeric(unlist(strsplit(input$mean, ",") ) ))
    x2 <- shiny::isolate(as.numeric(unlist(strsplit(input$cor, ",") ) ))
    x3 <- shiny::isolate(as.numeric(unlist(strsplit(input$sample, ",") ) ))
    
    dat.pl2 <- simu3(x1[1], x1[2], x1[3],x1[4],x1[5],x1[6],
                     x2[1], x2[2], x2[3], x3[1], x3[2], x3[3])
    

    
   gridExtra::grid.arrange(ppbound(ru =  as.numeric(input$rule), leg = FALSE, data = dat.pl2, meth = "Original", entro = TRUE ),
                 ppbound(ru =  as.numeric(input$rule), FALSE, data = dat.pl2, meth = "Modified" , entro = FALSE), 
                 ppbound(ru =  as.numeric(input$rule), FALSE, data = dat.pl2, meth = "Modified" , entro = TRUE), ncol = 3)
  })

  
  

    output$distPlot2 <- shiny::renderPlot({
      input$do2
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
      
      
      
      gridExtra::grid.arrange(ppbound(ru =  as.numeric(input$rule2), leg = FALSE, data = dat.pl2, meth = "Original", entro = TRUE ),
                   ppbound(ru =  as.numeric(input$rule2), FALSE, data = dat.pl2, meth = "Modified" , entro = FALSE), 
                   ppbound(ru =  as.numeric(input$rule2), FALSE, data = dat.pl2, meth = "Modified" , entro = TRUE), ncol = 3)
    })
  
  
  
   
}

shiny::shinyApp(ui = ui, server = server)
}




