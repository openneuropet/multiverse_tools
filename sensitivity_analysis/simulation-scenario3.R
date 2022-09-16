### simulation-scenario1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 15 2022 (16:48) 
## Version: 
## Last-Updated: sep 16 2022 (16:14) 
##           By: Brice Ozenne
##     Update #: 35
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## cd /projects/biostat01/people/hpl802/pipeline/
## source("simulation-scenario3.R")

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
n.sim <- c(250,100) ## number of simulations, number of times proportion is esitmated
n.cpus <- 50
Sigma <- as.matrix(bdiag((diag(0.05,15,15) + 0.95)*2.5, matrix(0.25), matrix(5), matrix(7.5), matrix(10), matrix(15)))
Sigmatot <- Sigma+diag(1,NROW(Sigma),NROW(Sigma))
## cov2cor(Sigmatot)

## * GLS weights
## solve(matrix(1, nrow = 1, ncol = NCOL(Sigma)) %*% solve(Sigma) %*% matrix(1, nrow = NCOL(Sigma), ncol = 1)) %*% matrix(1, nrow = 1, ncol = NCOL(Sigma)) %*% solve(Sigma)
sqrt(diag(Sigma))/min(sqrt(diag(Sigma)))

unique(round(100 * colSums(solve(Sigma))/sum(solve(Sigma)),2))
unique(round(100 * (1/diag(Sigma))/sum((1/diag(Sigma))),2))

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

ls.sim3 <- pblapply(1:n.sim[1], function(iSim){ ## iSim <- 1
    iOut <- NULL
    for(iN in seqN){ ## iN <- 10
        iOut <- rbind(iOut,
                      analyzeData(simData(n.obs = iN, sigma.pipe = Sigma, beta = 0), proportion = iSim<=n.sim[2]),
                      analyzeData(simData(n.obs = iN, sigma.pipe = Sigma, beta = beta), proportion = iSim<=n.sim[2])
                      )
    }
    return(cbind(sim=iSim,iOut))
}, cl = n.cpus)

## * Export
out3 <- do.call(rbind,ls.sim3)
rownames(out3) <- NULL
saveRDS(out3, file = "results/sim3.rds")
sessionInfo()


##----------------------------------------------------------------------
### simulation-scenario1.R ends here
