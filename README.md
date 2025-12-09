
<!-- README.md is generated from README.Rmd. Please edit that file -->
# PPtreeExt

PPtreeExt is an R package with extensons to the Projection Pursuit Tree (PPtree) algorithm to improve its performance in multi-class settings and under nonlinear separations.
The PPtree classifier finds separations between classes  based on linear combinations of variables by optimizing a projection pursuit index. One of its limitations is a rigid structure: the depth of a PPtree object is at most
$G$-1, where $G$ is the number of classes, with each class assigned to a single terminal node.
The proposed modifications enhance predictive performance in multi-class contexts, particularly in situations involving outliers or asymmetries. The objective is to increase the classifier's flexibility to handle more complex scenarios, while retaining interpretability.
The package includes an interactive web application to explore the behavior of the original and modified PPtree classifiers under a variety of scenarios.
This interactive tool played a key role in identifying limitations of the original algorithm and informing the design of the proposed modifications.


## Installation

Install version from CRAN (not available yet!): 

```r
install.packages("PPTreeExt")
```

Install the development version from GitHub:
```r
# install.packages("devtools")
devtools::install_github("natydasilva/PPtreeExt")
```
---

## Simple Example

```r
set.seed(249)
n <- nrow(iris)
tot <- c(1:n)
n.train <- round(n*0.9)
train <- sample(tot,n.train)
test <- tot[-train]
Tree.result <- PPtreeExtclass(formula = Species~., data = iris[train,],
PPmethod = "LDA",  srule = TRUE, tol = 0.1)
Tree.result
============================================================= 
Projection Pursuit Classification Tree Extension result 
=============================================================

1) root
   2)  proj1*X < cut1
      4)  proj2*X < cut2
         6)* proj3*X < cut3  ->  "virginica"
         7)* proj3*X >= cut3  ->  "virginica"
      5)  proj2*X >= cut2
         8)* proj4*X < cut4  ->  "versicolor"
         9)* proj4*X >= cut4  ->  "virginica"
   3)* proj1*X >= cut1  ->  "setosa"

Error rates 
-------------------------------------------------------------
[1] 1

 pred <- predict(object = Tree.result, newdata = iris[test,1:4], true.class = iris[test,5])
 
 
pred$predict.error
[1] 0

pred$predict.class
 [1] "setosa"     "setosa"     "setosa"     "setosa"     "setosa"     "versicolor"
 [7] "versicolor" "versicolor" "versicolor" "versicolor" "virginica"  "virginica" 
[13] "virginica"  "virginica"  "virginica" 

```

