
# compdag

<!-- badges: start -->
<!-- badges: end -->

The goal of compdag is to provide a score-based causal discovery algorithm to learn causality in compositional data.

## Installation

You can install the development version of cyclicvb from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("wangzhuofan/compdag")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(compdag)
## basic example code
# library(gtools)
generate_data<- function(n,px,py,signal,seed){
  set.seed(1)
  # n = 100; px = 30; py = 40;signal = 5
  m = matrix(0,py,px)
  for (colind in 1:px) {
    randNum = sample(1:3,1,replace = FALSE)
    randLoc = sample(1:py,randNum,replace = FALSE)
    m[randLoc,colind] = 1/randNum
  }
  set.seed(seed)
  x = t(rdirichlet(n,rep(0.1,px)))
  y <- matrix(0,py,n)
  temp =signal * m %*% x+0.1
  for (i in 1:n) {
    y[,i] = rdirichlet(1,temp[,i])
  }
  return(list(x=x,y=y))
}
data <- generate_data(50,5,5,5,1)
x = data$x
y = data$y
compdag(x,y,0.1,0.1)
```

