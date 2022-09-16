### simulation-scenario1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 15 2022 (16:48) 
## Version: 
## Last-Updated: sep 16 2022 (15:39) 
##           By: Brice Ozenne
##     Update #: 36
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## cd /projects/biostat01/people/hpl802/pipeline/
## source("simulation-scenario1.R")

## * Dependencies
library(mvtnorm)
library(data.table)
library(Matrix)
library(LMMstar)
library(pbapply)
source("simulation-FCT.R")

## * Parameters
seqN <- c(10, 25, 50, 100, 200)
beta <- 0.5
n.sim <- c(10000,250) ## number of simulations, number of times proportion is estimated
n.cpus <- 25
Sigma <- as.matrix(bdiag(diag(0.05,15,15) + 0.95,matrix(1),matrix(1),matrix(1),matrix(1),matrix(1)))
Sigmatot <- Sigma+diag(1,NROW(Sigma),NROW(Sigma))

## * GLS weights
## solve(matrix(1, nrow = 1, ncol = NCOL(Sigma)) %*% solve(Sigma) %*% matrix(1, nrow = NCOL(Sigma), ncol = 1)) %*% matrix(1, nrow = 1, ncol = NCOL(Sigma)) %*% solve(Sigma)
table(round(100*colSums(solve(Sigmatot))/sum(solve(Sigmatot)),2))
## 1.88 14.37 
##   15     5 

## * Simulation
##  test
## dt.ex <- simData(n.obs = 250, sigma.pipe = Sigma, beta = 0)
## analyzeData(dt.ex)

## dt.ex <- simData(n.obs = 100, sigma.pipe = Sigma, beta = beta)
## system.time(
##     analyzeData(dt.ex, proportion = TRUE)
## )
## system.time(
##     analyzeData(dt.ex, proportion = FALSE)
## )

ls.sim1 <- pblapply(1:n.sim[1], function(iSim){ ## iSim <- 1
    iOut <- NULL
    for(iN in seqN){ ## iN <- 200
        iOut <- rbind(iOut,
                      analyzeData(simData(n.obs = iN, sigma.pipe = Sigma, beta = 0), proportion = iSim<=n.sim[2]),
                      analyzeData(simData(n.obs = iN, sigma.pipe = Sigma, beta = beta), proportion = iSim<=n.sim[2])
                      )
    }
    return(cbind(sim=iSim,iOut))
}, cl = n.cpus)

## * Export
out1 <- do.call(rbind,ls.sim1)
rownames(out1) <- NULL
saveRDS(out1, file = "results/sim1.rds")
sessionInfo()


##----------------------------------------------------------------------
### simulation-scenario1.R ends here
