
utils::globalVariables(c("Sim", "X1", "X2", "predict", "pred", "ppred"))

# Description: Function to simulate multivariate (bivariate) normal distributions
# Input: 
# mean vector, 
# covariance matrix of the variables
# number of samples
# Output: 
# data.frame with dimension (n1+n2+n3)x(3)

simu3 <- function(mux1, mux2, muy1, muy2, muz1, muz2,
           cor1, cor2, cor3,
           n1 = 100, n2 = 100, n3 = 100
  ) {
    
    set.seed(666)
    bivn <- MASS::mvrnorm(n1, mu = c(mux1, mux2), Sigma = matrix(c(1, cor1, cor1, 1), 2))
    bivn2 <- MASS::mvrnorm(n2, mu = c(muy1, muy2), Sigma = matrix(c(1, cor2, cor2, 1), 2))
    bivn3 <- MASS::mvrnorm(n3, mu = c(muz1, muz2), Sigma = matrix(c(1, cor3, cor3, 1), 2))
    
    d1 <- data.frame(Sim = "sim1", bivn)
    d2 <- data.frame(Sim = "sim2", bivn2)
    d3 <- data.frame(Sim = "sim3", bivn3)
    rbind(d1, d2, d3)
  }

# Description:
# Input:
# ru = split rule = {1, 2, 3, 4, 5, 6, 7}
# data = dataset of simulate values
# meth = 
# entro = 
# Output:
# Generate grid values to use as test set.

ppbound <- function(ru, data , meth, entro , title, simM = FALSE) {

  n_grid <- 100
  
  grilla <-
    base::expand.grid(
      X1 = seq((min(data$X1) + sign(min(data$X1)) * .5), (max(data$X1) + sign(max(data$X1)) * .5), length.out = n_grid),
      X2 = seq((min(data$X2) + sign(min(data$X2)) * .5), (max(data$X2) + sign(max(data$X2)) * .5), length.out = n_grid)
    )
  
  if (meth == "Original") {
    # Si Sim no es factor se genera error. Agregar un chequeo de clases y que salte error/warning y se cambien.
    
    # Dependent variable must be a factor
    data$Sim <- as.factor(data$Sim)
    # Construct the projection pursuit classification tree (using LDA Index)
    pptree <- PPtreeViz::PPTreeclass(Sim ~ ., data = data, PPmethod = "LDA")
    # Predict projection pursuit classification tree ...NOT USED!!!
    ppred.sim <- PPtreeViz::PPclassify(pptree, test.data = grilla, Rule = ru)
    #PROBLEMA: POR QUE ppred.sim$predict.error es NA?
    grilla$pred <- ppred.sim[[2]]
    err <-
      round(
        PPtreeViz::PPclassify(
          pptree,
          test.data = data[, -1],
          true.class = data[, 1],
          Rule = ru
        )[[1]] / nrow(data),
        3
      ) * 100
  }
  if (meth == "Rpart") {
    rpart.mod <- rpart::rpart(Sim ~ ., data = data)
    grilla$pred <- predict(rpart.mod, newdata = grilla, type = "class")
    err <-
      round(1 - sum(diag(table(
        predict(rpart.mod, newdata = data[, -1], type = "class") , data[, 1]
      ))) / nrow(data), 3) * 100
  }
  
  # Why this???
  if (entro) {
    mod = 2
  } else{
    mod = 1
  }
  
  if (meth == "Modified") {
    
    # Projection pursuit classification tree with random variable selection in each split
    pptree <- PPtree_splitMOD(Sim ~ ., data = data, "LDA", entro = entro)
    # 
    ppred.sim <- PPtreeViz::PPclassify(pptree, test.data = grilla, Rule = ru)
    grilla$pred <- paste("sim", ppred.sim[[2]], sep = "")
    err <-
      round(
        PPtreeViz::PPclassify(
          pptree,
          test.data = data[, -1],
          true.class = data[, 1],
          Rule = ru
        )[[1]] / nrow(data),
        3
      ) * 100
    
  }
  
  #ruleid <- pptree$splitCutoff.node[,ru]
  if (simM) {
    pl.pp  <-
      ggplot2::ggplot(data = grilla) + 
      ggplot2::geom_point(ggplot2::aes(
        x = X1,
        y = X2,
        color = as.factor(pred)
      ), alpha = .20) +
      ggplot2::scale_colour_brewer(name = "Class",
                                   type = "qual",
                                   palette = "Dark2") + ggplot2::theme_bw() +
      ggplot2::geom_point(
        data = data,
        ggplot2::aes(
          x = X1 ,
          y = X2,
          group = Sim,
          color = Sim
        ),
        size = I(3)
      ) +
      ggplot2::theme(legend.position = "none", aspect.ratio = 1) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) + 
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::labs(
        x = " ",
        y = "",
        title = paste(title, "error", err, "%")
      )
  } else {
    pl.pp  <-
      ggplot2::ggplot(data = grilla) + 
      ggplot2::geom_point(ggplot2::aes(
        x = X1,
        y = X2,
        color = as.factor(pred),
        shape = as.factor(pred)
      ),
      alpha = .20) +
      ggplot2::scale_colour_brewer(name = "Class",
                                   type = "qual",
                                   palette = "Dark2") + 
      ggplot2::theme_bw() +
      ggplot2::scale_shape_discrete(name = 'Class') + 
      ggplot2::geom_point(
        data = data,
        ggplot2::aes(
          x = X1 ,
          y = X2,
          group = Sim,
          shape = Sim,
          color = Sim
        ),
        size = I(3)
      ) +
      ggplot2::theme(legend.position = "none", aspect.ratio = 1) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) + 
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::labs(
        x = " ",
        y = "",
        title = paste(title, "error", err, "%")
      )
  }
  
  pl.pp
}

ppboundMOD <-
  function(data ,
           meth = "MOD",
           entro = FALSE,
           entroindiv = TRUE,
           title,
           simM = FALSE,
           strule,
           tot) {
    
    # Construct grid values
    grilla <-
      base::expand.grid(
        X1 = seq((min(data$X1) + sign(min(data$X1)) * .5), (max(data$X1) + sign(max(data$X1)) * .5), length.out = 100), 
        X2 = seq((min(data$X2) + sign(min(data$X2)) * .5), (max(data$X2) + sign(max(data$X2)) * .5), length.out = 100)
      )
    # Responde variable must be a factor
    data$Sim <- as.factor(data$Sim)
    
    # Construct the projection pursuit classification tree using LDA index.
    pptree <-
      PPTreeclass_MOD(
        Sim ~ . ,
        data = data,
        PPmethod = 'LDA',
        strule = strule,
        tot = tot
      )
    
    ppred.sim <- PPclassify_MOD(pptree, test.data = grilla)
    
    grilla$ppred <- ppred.sim[[2]]
    
    err <-
      round(PPclassify_MOD(pptree, test.data = data[, -1], true.class = data[, 1])[[1]] /
              nrow(data),
            3) * 100
    
    if (simM) {
      pl.pp <-
        ggplot2::ggplot(data = grilla) + 
        ggplot2::geom_point(ggplot2::aes(x = X1, y = X2, color = ppred), alpha = .20) +
        ggplot2::scale_colour_brewer(name = "Class",
                                     type = "qual",
                                     palette = "Dark2") + 
        ggplot2::theme_bw() +
        ggplot2::scale_shape_discrete(name = 'Class') + 
        ggplot2::geom_point(
          data = data,
          ggplot2::aes(
            x = X1 ,
            y = X2,
            group = Sim,
            color = Sim
          ),
          size = I(3)
        ) +
        ggplot2::theme(legend.position = "none", aspect.ratio = 1) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) + 
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::labs(
          x = " ",
          y = "",
          title = paste(title, "error", err, "%")
        )
      
    } else {
      pl.pp <-
        ggplot2::ggplot(data = grilla) + 
        ggplot2::geom_point(ggplot2::aes(
          x = X1,
          y = X2,
          color = ppred,
          shape = ppred
        ), 
        alpha = .20) +
        ggplot2::scale_colour_brewer(name = "Class",
                                     type = "qual",
                                     palette = "Dark2") + 
        ggplot2::theme_bw() +
        ggplot2::scale_shape_discrete(name = 'Class') + 
        ggplot2::geom_point(
          data = data,
          ggplot2::aes(
            x = X1 ,
            y = X2,
            group = Sim,
            shape = Sim,
            color = Sim
          ),
          size = I(3)
        ) +
        ggplot2::theme(legend.position = "none", aspect.ratio = 1) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) + 
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::labs(
          x = " ",
          y = "",
          title = paste(title, "error", err, "%")
        )
      
    }
    pl.pp
  }

# Comentó esto y lo agrego en el server.
# Mejor dicho...Para qué esta esto si no se asigna a ningún objeto ???
# ppboundMOD(
#   data = dat.pl2, # Datos simulados en server. 
#   meth = "MOD",
#   entro = FALSE,
#   entroindiv = TRUE,
#   title = "Modified multi_sp",
#   strule = x4,
#   tot = sum(x3)
# )