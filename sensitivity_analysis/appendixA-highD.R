### appendixA-highD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 12 2022 (16:13) 
## Version: 
## Last-Updated: okt 10 2022 (12:29) 
##           By: Brice Ozenne
##     Update #: 18
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
J <- 60
rho <- 0.5
mu <- 1:J
Sigma <- rho + diag(1-rho,J,J)
dmu <- rep(0.2,J)

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

## * Univariate analyses
## ** fit linear regression
ls.lm <- lapply(paste0("V",1:J), function(iPipe){lm(value~group, data = dtL[pipeline==iPipe])})

## ** extract summary statistics
e.mlmm <- mlmm(value~group, repetition = ~1|id, data = dtL, df = FALSE, robust = TRUE,
                by = "pipeline", effects = "groupG2=0")
plot(e.mlmm)

## * Pooling
psi <- coef(e.mlmm)
Sigma_psi <- vcov(e.mlmm)

## average
model.tables(e.mlmm, method = "average")
##                              estimate        se  df     lower     upper   p.value
## pipeline=<V1:V60>: groupG2 0.06388668 0.2289877 Inf -0.384921 0.5126944 0.7802478
mean(psi)
## [1] 0.06388668

## pool according to uncertainty
model.tables(e.mlmm, method = "pool.se")
##                              estimate        se  df      lower     upper   p.value
## [1] 0.07305267
weighted.mean(psi, w = 1/diag(Sigma_psi)) ## robust se
## [1] 0.07305267

## pool according to uncertainty and correlation (GLS)
model.tables(e.mlmm, method = "pool.gls")
##                               estimate         se  df      lower     upper   p.value
## pipeline=<V1:V60>: groupG2 -0.08840572 0.09647198 Inf -0.2774873 0.1006759 0.3594637
eigen_psi <- eigen(Sigma_psi)
weight <- colSums(eigen_psi$vectors)^2/eigen_psi$values
eigenvector.norm <- sweep(eigen_psi$vectors, FUN = "/", MARGIN = 2, STATS = colSums(eigen_psi$vectors))
weighted.mean(psi %*% eigenvector.norm[,eigen_psi$values>1e-12], weight[eigen_psi$values>1e-12])

##----------------------------------------------------------------------
### appendixA-highD.R ends here
