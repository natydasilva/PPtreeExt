
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis-CI Build Status](https://travis-ci.org/natydasilva/PPtreeExt.svg?branch=master)](https://travis-ci.org/natydasilva/PPtreeExt)

# PPtreeExt
PPtreeExt is an R package that extends the PPtree algorithm  for classification problems. The extentions allows 
to handle multiple splits per group and more complex structures than the original algorithm. 


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
n <- nrow(iris)
tot <- c(1:n)
n.train <- round(n*0.9)
train <- sample(tot,n.train)
test <- tot[-train]
Tree.result <- PPTreeclass_MOD(formula = Species~.,data = iris[train,],PPmethod = "LDA")
Tree.result
============================================================= 
Projection Pursuit Classification Tree result 
=============================================================

1) root
   2)  proj1*X < cut1
      4)  proj2*X < cut2
         6)* proj3*X < cut3  ->  "virginica"
         7)  proj3*X >= cut3
            8)  proj4*X < cut4
               10)  proj5*X < cut5
                  12)* proj6*X < cut6  ->  "versicolor"
                  13)  proj6*X >= cut6
                     14)* proj7*X < cut7  ->  "versicolor"
                     15)  proj7*X >= cut7
                        16)  proj8*X < cut8
                           18)* proj9*X < cut9  ->  "virginica"
                           19)  proj9*X >= cut9
                              20)  proj10*X < cut10
                                 22)* proj11*X < cut11  ->  "versicolor"
                                 23)  proj11*X >= cut11
                                    24)  proj12*X < cut12
                                       26)* proj13*X < cut13  ->  "versicolor"
                                       27)* proj13*X >= cut13  ->  "versicolor"
                                    25)  proj12*X >= cut12
                                       28)* proj14*X < cut14  ->  "virginica"
                                       29)  proj14*X >= cut14
                                          30)* proj15*X < cut15  ->  "virginica"
                                          31)* proj15*X >= cut15  ->  "virginica"
                              21)* proj10*X >= cut10  ->  "versicolor"
                        17)* proj8*X >= cut8  ->  "versicolor"
               11)* proj5*X >= cut5  ->  "versicolor"
            9)  proj4*X >= cut4
               32)* proj16*X < cut16  ->  "setosa"
               33)  proj16*X >= cut16
                  34)  proj17*X < cut17
                     36)* proj18*X < cut18  ->  "versicolor"
                     37)  proj18*X >= cut18
                        38)* proj19*X < cut19  ->  "setosa"
                        39)* proj19*X >= cut19  ->  "setosa"
                  35)* proj17*X >= cut17  ->  "setosa"
      5)* proj2*X >= cut2  ->  "virginica"
   3)* proj1*X >= cut1  ->  "virginica"

Error rates of various cutoff values 
-------------------------------------------------------------
[1] NA

PPclassify_MOD(Tree.result,test.data = iris[test,1:4], true.class = iris[test,5])
$predict.error
[1] 2

$predict.class
 [1] "setosa"     "setosa"     "setosa"     "setosa"     "setosa"     "setosa"    
 [7] "setosa"     "setosa"     "virginica"  "versicolor" "versicolor" "setosa"    
[13] "virginica"  "virginica"  "virginica" 

```

