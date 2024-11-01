### BATCH_simulation-scenario4.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 15 2022 (16:48) 
## Version: 
## Last-Updated: nov  1 2024 (11:35) 
##           By: Brice Ozenne
##     Update #: 77
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## cd /projects/biostat01/people/hpl802/pipeline/
## ## INTERACTIVE
## source("BATCH_simulation-scenario4.R")
## ## BATCH
## for ITER in `seq 1 10`;
## do
## eval 'R CMD BATCH --vanilla "--args iter_sim='$ITER' n.iter_sim=10" BATCH_simulation-scenario4.R output/simulation-scenario4/simulation-scenario4-'$ITER'.Rout &'
## done
## 421816-421825

## * SLURM
args <- commandArgs(TRUE) ## BATCH MODE
if(length(args)>0){
    for (arg in args){
        eval(parse(text=arg))
    }
}else{ ## SLUMR
    iter_sim <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
    n.iter_sim <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
}
if(is.na(iter_sim)){iter_sim <- 1}
if(is.na(n.iter_sim)){n.iter_sim <- 10}

set.seed(1)
seqSeed <- sample(1:max(1e5,n.iter_sim),size=n.iter_sim,replace=FALSE)
iSeed <- seqSeed[iter_sim]
set.seed(iSeed)

cat("iteration ",iter_sim," over ",n.iter_sim," (seed ",iSeed,")\n",sep="")

## * Prepare export
path <- "."
path.res <- file.path(path,"Results","simulation-scenario4")
if(dir.exists(path.res)==FALSE){
    if(dir.exists(file.path(path,"Results"))==FALSE){
    dir.create(file.path(path,"Results"))
    }
    dir.create(path.res)
}
path.output <- file.path(path,"output","simulation-scenario4")
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
n.sim <- c(100,10) ## number of simulations, number of times proportion is estimated
n.cpus <- 25
Sigma <- as.matrix(bdiag((diag(0.05,15,15) + 0.95)*2.5, matrix(0.25), matrix(5), matrix(7.5), matrix(10), matrix(15)))
Sigmatot <- Sigma+diag(1,NROW(Sigma),NROW(Sigma))

## * Weights
## solve(matrix(1, nrow = 1, ncol = NCOL(Sigma)) %*% solve(Sigma) %*% matrix(1, nrow = NCOL(Sigma), ncol = 1)) %*% matrix(1, nrow = 1, ncol = NCOL(Sigma)) %*% solve(Sigma)
cat("\n Relative variance: ")
relativeVar <- sqrt(diag(Sigmatot))/min(sqrt(diag(Sigmatot)))
print(table(relativeVar))
##  1 
## 20 

cat("\n Pool weights: ")
Sigma2weight <- round(100*(1/diag(Sigmatot))/sum((1/diag(Sigmatot))),2)
print(table(unname(Sigma2weight)))
##  5 
## 20 

cat("\n GLS weights: ")
Sigma2weight <- round(100*colSums(solve(Sigmatot))/sum(solve(Sigmatot)),2)
print(table(unname(Sigma2weight)))
## 1.88 14.37 
##   15     5 
cat("\n")
## analyzeData(simData(n.obs = 1000, sigma.pipe = Sigma, beta = 0), proportion = FALSE, print.weight = TRUE)

## * Simulation
res <- NULL
for(iSim in 1:n.sim[1]){ ## 
    cat(iSim,": ")
    for(iN in seqN){ ## iN <- 200
        cat(iN," ")
        try(res <- rbind(res,
                         cbind(sim = iSim, seed = iSeed, analyzeData(simData(n.obs = iN, sigma.pipe = Sigma, beta = 0, df = 3), proportion = iSim<=n.sim[2]))))
        try(res <- rbind(res,
                         cbind(sim = iSim, seed = iSeed, analyzeData(simData(n.obs = iN, sigma.pipe = Sigma, beta = beta, df = 3), proportion = iSim<=n.sim[2]))))
    }
    cat("\n")
    saveRDS(res, file = file.path(path.res,paste0("scenario4_",iter_sim,"(tempo).rds")))
}
cat("\n")
    
## * Export
rownames(res) <- NULL
saveRDS(res, file = file.path(path.res,paste0("scenario4_",iter_sim,".rds")))
sessionInfo()


##----------------------------------------------------------------------
### BATCH_simulation-scenario4.R ends here
