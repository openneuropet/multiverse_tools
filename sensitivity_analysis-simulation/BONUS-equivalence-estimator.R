### BONUS-equivalence-estimator.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 13 2022 (10:51) 
## Version: 
## Last-Updated: sep 14 2022 (17:34) 
##           By: Brice Ozenne
##     Update #: 23
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Packages
library(mvtnorm)
library(ggplot2)
library(data.table)
library(lava)
library(LMMstar) ## require development version >= 0.8.0 which can be installed using remotes::install_github("bozenne/LMMstar")

## * Simulate data
## ** settings
n <- 25
J <- 5
rho <- 0.5
mu <- 1:J
Sigma <- rho + diag(1-rho,J,J)
dmu <- rep(0.2,J)

## Sigma <- tcrossprod(1:5) * rbind(c(1,0.5,0.25,0.125,0),
##                                  c(0.5,1,0.5,0.25,0.125),
##                                  c(0.25,0.5,1,0.5,0.25),
##                                  c(0.125,0.25,0.5,1,0.5),
##                                  c(0,0.125,0.25,0.5,1))
## dmu <- 1:J 

## ** simulation
set.seed(10)
Y_G1 <- rmvnorm(n, mean = mu, sigma = Sigma)
Y_G2 <- rmvnorm(n, mean = mu + dmu, sigma = Sigma)

## ** reshape
dtW <- rbind(data.table(id = paste0("H",formatC(1:n, width = 3, format = "d", flag = "0")),
                         group = "G1",
                         Y_G1),
              data.table(id = paste0("C",formatC(1:n, width = 3, format = "d", flag = "0")),
                         group = "G2",
                         Y_G2))
dtL <- melt(dtW, id.vars = c("id","group"), variable.name = "pipeline")
dtL
setkeyv(dtL,c("id","pipeline"))
## dtL <- dtL[-c(1,10,50,89,150,249)]

## * PCA formulation
e.mlmm <- mlmm(value~group, repetition = ~1|id, data = dtL, df = FALSE, robust = TRUE,
                by = "pipeline", effects = "groupG2=0")

## ** GLS estimator
oneGLS <- matrix(1,J,1)
SigmaGLS <- vcov(e.mlmm)
betaGLS <- coef(e.mlmm)
solve(t(oneGLS) %*% solve(SigmaGLS) %*% oneGLS) %*% t(oneGLS) %*% solve(SigmaGLS) %*% betaGLS
##           [,1]
## [1,] 0.2620123

## ** PCA estimator (automatic)
model.tables(e.mlmm, method = "pool.pca")
##                                     estimate        se  df      lower     upper   p.value
## pipeline=<V1,V2,V3,V4,V5>: groupG2 0.2620123 0.2451637 Inf -0.2184997 0.7425243 0.2851941

## ** PCA estimator (by hand)
iEigen <- eigen(SigmaGLS)
iLambda <- iEigen$values
iP <- iEigen$vectors
iPsum <- colSums(iEigen$vector)
iPstar <- sweep(iP, FUN = "/", MARGIN = 2, STATS = iPsum)
iW <- iPsum^2/iLambda
weighted.mean(betaGLS %*% iPstar, w = iW)
## [1] 0.2620123

## ** equivalence left term
iX <- matrix(1,J,1)
sum(colSums(iP)^2/iEigen$values)
sum( (t(iX) %*% iP)^2 /iEigen$values)
sum( (t(iX) %*% iP)^2 /iEigen$values)
sum(1/iOmega)

t(oneGLS) %*% solve(SigmaGLS) %*% oneGLS

## ** equivalence right term
iPsum %*% solve(diag(iEigen$values)) %*% t(iP)
sweep(iP, FUN = "*", MARGIN = 2, STATS = iPsum/iEigen$values) %*% matrix(1,nrow=5,ncol = 1)

## ** overall
nu <- 1/sum( colSums(iEigen$vectors)^2 / iEigen$values)
nu * rowSums(sweep(iEigen$vectors, FUN = "*", MARGIN = 2, STATS = colSums(iEigen$vectors)/iEigen$values)) %*% betaGLS

## * Equivalence stratified analysis + pooling with GLS and joint model

## ** overview
e.mlmm <- mlmm(value~group, repetition = ~1|id, data = dtL, df = FALSE, robust = TRUE,
                by = "pipeline", effects = "groupG2=0")
e.glsUN <- lmm(value ~ pipeline + group, repetition = ~pipeline|id, data = dtL,
               structure = "UN")


#### joint model
coef(e.glsUN)["groupG2"]
##   groupG2 
## 0.2620113 

#### gls estimator
oneGLS <- matrix(1,J,1)
SigmaGLS <- vcov(e.mlmm)
betaGLS <- coef(e.mlmm)
solve(t(oneGLS) %*% solve(SigmaGLS) %*% oneGLS) %*% t(oneGLS) %*% solve(SigmaGLS) %*% betaGLS
##           [,1]
## [1,] 0.2620123

#### pca estimator
model.tables(e.mlmm, method = "pool.pca")
##                                     estimate        se  df      lower     upper   p.value
## pipeline=<V1,V2,V3,V4,V5>: groupG2 0.2620123 0.2451637 Inf -0.2184997 0.7425243 0.2851941

## ** detail

## precompute functions of X
X <- model.matrix(e.mlmm$model[[1]])
XXm1 <- solve(crossprod(X)) ## solve(n, n/2, n/2, n/2)
                            ##  1/(n^2/4)  (n/2, -n/2, -n/2, n) = (2/n, -2/n, -2/n, 4/n)
XXXm1 <- X %*% XXm1 ## second element -2/n or 2/n
table(XXXm1[,2])

XXXm12 <- XXXm1^2 %*% c(0,1) ##  -2/n or 2/n squared, i.e. 4/n^2
table(XXXm12)

## betaGLS
Y <- do.call(cbind,tapply(dtL$value,dtL$pipeline,function(x){x}, simplify = FALSE))

solve(crossprod(X)) %*% t(X) %*% Y
colMeans(Y[X[,"groupG2"]==1,,drop=FALSE]) - colMeans(Y[X[,"groupG2"]==0,,drop=FALSE])
coef(e.mlmm)

## SigmaGLS
epsilon <- do.call(cbind,lapply(e.mlmm$model, residuals))
vcov(e.mlmm) - crossprod(epsilon) * XXXm12[1]



##----------------------------------------------------------------------
### BONUS-equivalence-estimator.R ends here
