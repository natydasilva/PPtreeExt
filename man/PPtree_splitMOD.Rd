% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PPtree_splitMOD.R
\name{PPtree_splitMOD}
\alias{PPtree_splitMOD}
\title{Projection pursuit classification tree with random variable selection in each split}
\usage{
PPtree_splitMOD(form, data, PPmethod='LDA', 
size.p=1,  lambda=0.1, entro ,entroindiv,...)
}
\arguments{
\item{form}{A character with the name of the class variable.}

\item{data}{Data frame with the complete data set.}

\item{PPmethod}{index to use for projection pursuit: 'LDA', 'PDA'}

\item{size.p}{proportion of variables randomly sampled in each split, default is 1, returns a PPtree.}

\item{lambda}{penalty parameter in PDA index and is between 0 to 1 . If \code{lambda = 0}, no penalty parameter is added and the PDA index is the same as LDA index. If \code{lambda = 1} all variables are treated as uncorrelated. The default value is \code{lambda = 0.1}.}

\item{entro}{TRUE, compute the entropy method}

\item{entroindiv}{TRUE, compute the entropy for each obs...clarify this}

\item{...}{arguments to be passed to methods}
}
\value{
An object of class \code{PPtreeclass} with components
\item{Tree.Struct}{Tree structure of projection pursuit classification tree}
\item{projbest.node}{1-dim optimal projections of each split node}
\item{splitCutoff.node}{cutoff values of each split node}
\item{origclass}{original class} 
\item{origdata}{original data}
}
\description{
Find tree structure using various projection pursuit indices of classification in each split.
}
\examples{
#crab data set
\dontrun{
train<- sample(1:200,150)
Tree.crab <- PPtree_splitMOD("Type~.", data = PPforest::crab[train, ],
 PPmethod = "LDA", size.p = 1, entro = TRUE,entroindiv=FALSE)
Tree.crab

 Tree.result <- PPtreeViz::PPTreeclass(Type~.,data = PPforest::crab[train,],"LDA")
Tree.result

PPtreeViz::PPclassify(Tree.result,PPforest::crab[-train,-1],1,crab[-train,1])
Tree.iris <- PPtree_splitMOD("Species~.", data = iris, PPmethod = "PDA", 
size.p = 1, entro = TURE, entroindiv = FALSE)
Tree.iris}
}
\references{
Lee, YD, Cook, D., Park JW, and Lee, EK (2013) 
PPtree: Projection pursuit classification tree, 
Electronic Journal of Statistics, 7:1369-1386.
}
\keyword{tree}
