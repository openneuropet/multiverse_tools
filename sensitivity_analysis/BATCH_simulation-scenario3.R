### BATCH_simulation-scenario3.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 15 2022 (16:48) 
## Version: 
## Last-Updated: okt  5 2022 (13:46) 
##           By: Brice Ozenne
##     Update #: 49
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## cd /projects/biostat01/people/hpl802/pipeline/
## source("BATCH_simulation-scenario3.R")

## * SLURM
iter_sim <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
n.iter_sim <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
if(is.na(iter_sim)){iter_sim <- 1}
if(is.na(n.iter_sim)){n.iter_sim <- 10}

set.seed(1)
seqSeed <- sample(1:max(1e5,n.iter_sim),size=n.iter_sim,replace=FALSE)
iSeed <- seqSeed[iter_sim]
set.seed(iSeed)

cat("iteration ",iter_sim," over ",n.iter_sim," (seed ",iSeed,")\n",sep="")

## * Prepare export
path <- "."
path.res <- file.path(path,"Results","simulation-scenario3")
if(dir.exists(path.res)==FALSE){
    if(dir.exists(file.path(path,"Results"))==FALSE){
    dir.create(file.path(path,"Results"))
    }
    dir.create(path.res)
}
path.output <- file.path(path,"output","simulation-scenario3")
if(dir.exists(path.output)==FALSE){
    if(dir.exists(file.path(path,"output"))==FALSE){
    dir.create(file.path(path,"output"))
    }
    dir.create(path.output)
}

## * Dependencies
data.table::setDTthreads(1)
library(data.table)
library(mvtnorm)
library(Matrix)
library(LMMstar)
library(lava)
library(pbapply)
source("FCT_simulation.R")

## * Parameters
seqN <- c(10, 25, 50, 100, 200, 500)
beta <- 0.5
n.sim <- c(1000,25) ## number of simulations, number of times proportion is esitmated
n.cpus <- 50
Sigma <- as.matrix(bdiag((diag(0.05,15,15) + 0.95)*2.5, matrix(0.25), matrix(5), matrix(7.5), matrix(10), matrix(15)))
Sigmatot <- Sigma+diag(1,NROW(Sigma),NROW(Sigma))
## cov2cor(Sigmatot)

## * GLS weights
## solve(matrix(1, nrow = 1, ncol = NCOL(Sigma)) %*% solve(Sigma) %*% matrix(1, nrow = NCOL(Sigma), ncol = 1)) %*% matrix(1, nrow = 1, ncol = NCOL(Sigma)) %*% solve(Sigma)
cat("\n Relative variance: ")
relativeVar <- sqrt(diag(Sigmatot))/min(sqrt(diag(Sigmatot)))
print(table(round(relativeVar,4)))
     ## 1 1.6733 2.1909 2.6077 2.9665 3.5777 
     ## 1     15      1      1      1      1 

cat("\n Pool weights: ")
Sigma2weight <- round(100*(1/diag(Sigmatot))/sum((1/diag(Sigmatot))),2)
print(table(unname(Sigma2weight)))
## 1.13  1.65  2.13  3.02  5.17 14.48 
##    1     1     1     1    15     1 

cat("\n GLS weights: ")
Sigma2weight <- round(100*colSums(solve(Sigmatot))/sum(solve(Sigmatot)),2)
print(table(unname(Sigma2weight)))
## 1.65   3.8  5.52  7.15 10.13 48.61 
##   15     1     1     1     1     1 

## * Simulation
res <- NULL
for(iSim in 1:n.sim[1]){
    cat(iSim,": ")
    for(iN in seqN){ ## iN <- 200
        cat(iN," ")
        res <- rbind(res,
                     cbind(sim = iSim, seed = iSeed, analyzeData(simData(n.obs = iN, sigma.pipe = Sigma, beta = 0), proportion = iSim<=n.sim[2])),
                     cbind(sim = iSim, seed = iSeed, analyzeData(simData(n.obs = iN, sigma.pipe = Sigma, beta = beta), proportion = iSim<=n.sim[2]))
                     )
    }
    cat("\n")
    saveRDS(res, file = file.path(path.res,paste0("scenario3_",iter_sim,"(tempo).rds")))
}
cat("\n")

## * Export
rownames(res) <- NULL
saveRDS(res, file = file.path(path.res,paste0("scenario3_",iter_sim,".rds")))
sessionInfo()


##----------------------------------------------------------------------
### BATCH_simulation-scenario3.R ends here
