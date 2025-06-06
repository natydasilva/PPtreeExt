#' PP optimization using PDA index same as PPtree
#' 
#' Find the q-dimensional optimal projection using PDA projectin pursuit index
#' @title PP optimization using PDA index
#' @usage PDAopt_MOD(origclass,origdata, q = 1, weight = TRUE, lambda = 0.1,...) 
#' @param origclass class information vector of data
#' @param origdata data matrix without class information
#' @param q dimension of projection vector
#' @param weight weight flag in PDA index
#' @param lambda lambda in PDA index
#' @param ... arguments to be passed to methods
#' @return indexbest maximum PDA index value
#' @return projbest optimal q-dimensional projection matrix
#' @return origclass original class information vector
#' @return origdata  original data matrix without class information
#' @references Lee, EK, Cook, D.(2010) 
#' A Projection Pursuit Index for Large p Small n Data, 
#' Statistics and Computing, 20:381-392.
#' @useDynLib PPtreeExt
#' @importFrom Rcpp evalCpp
#' @export
#' @keywords projection pursuit
#' @examples
#' data(penguins)
#' penguins <- na.omit(penguins[, -c(2,7)])
#' 
#' penguins_pda_proj <- PDAopt_MOD(penguins[,1],penguins[,-1], weight=TRUE, q=2, lambda=0.1)
#' penguins_pda_proj$indexbest
#' penguins_pda_proj$projbest
PDAopt_MOD <-
  function(origclass,
           origdata,
           q = 1,
           weight = TRUE,
           lambda = 0.1,
           ...) {
    data.std <- as.matrix(origdata)
    class.table <- table(origclass)
    g <- length(class.table)
    class.name <- names(class.table)
    p <- ncol(data.std)
    n <- nrow(data.std)

    mean.g <- matrix(apply(data.std, 2, function(x) {
      tapply(x, origclass, mean, na.rm = TRUE)
    }),
    ncol = p)
   
    mean.all <- matrix(apply(data.std, 2, mean), ncol = p)
    B <- matrix(0, ncol = p, nrow = p)
    W <- matrix(0, ncol = p, nrow = p)
    for (i in 1:length(class.table)) {
      sel.id <- which(origclass == class.name[i])
      temp.m1 <- mean.g[i, ] - mean.all
      temp.m2 <- data.std[sel.id, ] -
        matrix(1, length(sel.id), ncol = 1) %*% mean.g[i, , drop = FALSE]
      gn1 <- ifelse(weight, length(sel.id), n / g)
      B <- B + gn1 * t(temp.m1) %*% temp.m1
      W <- W + gn1 * t(temp.m2) %*% temp.m2 / length(sel.id)
    }
    
    W.t <- (1 - lambda) * W
    diag(W.t) <- diag(W)
    WB.t <- W.t + B
    opt <- eigen(MASS::ginv(WB.t) %*% B)
    
    optVector <- matrix(as.numeric(opt$vectors[, 1:q]), ncol = q)
    proj.data <- data.std %*% optVector
    optindex <-
      PDAindex2(origclass, proj.data, optVector, weight, lambda = lambda)
    optobj <-
      list(
        indexbest = optindex,
        projbest = optVector,
        origclass = origclass,
        origdata = origdata,
        data.std = data.std
      )
    class(optobj) <- append(class(optobj), "PPoptim")
    return(optobj)
  }
