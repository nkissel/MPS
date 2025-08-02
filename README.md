# Model Path Selection
The "kisselmentch2021" folder contains R code to reproduce the simulations, applications, and motivation examples in paper [Forward Stability and Model Path Selection](https://link.springer.com/article/10.1007/s11222-024-10395-8).

The "ModelPath" folder can be installed as an R package. It provides tools for performing MPS.
 
 ## R Package Install
`devtools::install_github(repo='nkissel/MPS', subdir='ModelPath')`

# Getting started
For an example of how to use MPS
```
library(MASS)
n <- 1000
p <- 10
x <- mvrnorm(n, rep(0, p), diag(p))
signal <- rowSums(x[,1:3])
noise <- rnorm(n, 0, 1/4)
y <- signal + noise
mydata <- data.frame(x, y)
mps1 <- full.select.gen(myframe = mydata, resp.name = 'y', depth = 3,
   r = 100, model = 'lm', condense = FALSE)
mps1
build.tree(mps1) #graph

#merged paths
mps2 <- full.select.gen(myframe = mydata, resp.name = 'y', depth = 3,
   r = 100, model = 'lm', condense = TRUE)
mps2
build.tree(mps2) #graph
```
