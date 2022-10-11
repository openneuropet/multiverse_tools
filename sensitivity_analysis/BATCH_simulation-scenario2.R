### BATCH_simulation-scenario2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 15 2022 (16:48) 
## Version: 
## Last-Updated: okt 10 2022 (15:01) 
##           By: Brice Ozenne
##     Update #: 64
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
## source("BATCH_simulation-scenario2.R")
## ## BATCH
## for ITER in `seq 1 10`;
## do
## eval 'R CMD BATCH --vanilla "--args iter_sim='$ITER' n.iter_sim=10" BATCH_simulation-scenario2.R output/simulation-scenario2/simulation-scenario2-'$ITER'.Rout &'
## done
## 346982-346991

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
path.res <- file.path(path,"Results","simulation-scenario2")
if(dir.exists(path.res)==FALSE){
    if(dir.exists(file.path(path,"Results"))==FALSE){
    dir.create(file.path(path,"Results"))
    }
    dir.create(path.res)
}
path.output <- file.path(path,"output","simulation-scenario2")
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
n.sim <- c(1000,25) ## number of simulations, number of times proportion is estimated
n.cpus <- 25
Sigma <- diag(c(2.5,0.25,5,7.5,10,15))
Sigmatot <- Sigma+diag(1,NROW(Sigma),NROW(Sigma))

## * GLS weights!
## solve(matrix(1, nrow = 1, ncol = NCOL(Sigma)) %*% solve(Sigma) %*% matrix(1, nrow = NCOL(Sigma), ncol = 1)) %*% matrix(1, nrow = 1, ncol = NCOL(Sigma)) %*% solve(Sigma)
cat("\n Relative variance: ")
relativeVar <- sqrt(diag(Sigmatot))/min(sqrt(diag(Sigmatot)))
print(relativeVar)
## [1] 1.673320 1.000000 2.190890 2.607681 2.966479 3.577709

cat("\n Pool weights: ")
Sigma2weight <- round(100*(1/diag(Sigmatot))/sum((1/diag(Sigmatot))),2)
print(table(unname(Sigma2weight)))
## 4.1  5.97  7.72 10.94 18.75 52.51 
##   1     1     1     1     1     1 

cat("\n GLS weights: ")
Sigma2weight <- round(100*colSums(solve(Sigmatot))/sum(solve(Sigmatot)),2)
print(table(unname(Sigma2weight)))
## 4.1  5.97  7.72 10.94 18.75 52.51 
##   1     1     1     1     1     1 

cat("\n")

## * Simulation
res <- NULL
for(iSim in 1:n.sim[1]){
    cat(iSim,": ")
    for(iN in seqN){ ## iN <- 200
        cat(iN," ")
        try(res <- rbind(res,
                         cbind(sim = iSim, seed = iSeed, analyzeData(simData(n.obs = iN, sigma.pipe = Sigma, beta = 0), proportion = iSim<=n.sim[2]))))
        try(res <- rbind(res,
                         cbind(sim = iSim, seed = iSeed, analyzeData(simData(n.obs = iN, sigma.pipe = Sigma, beta = beta), proportion = iSim<=n.sim[2]))))
    }
    cat("\n")
    saveRDS(res, file = file.path(path.res,paste0("scenario2_",iter_sim,"(tempo).rds")))
}
cat("\n")

## * Export
rownames(res) <- NULL
saveRDS(res, file = file.path(path.res,paste0("scenario2_",iter_sim,".rds")))
sessionInfo()


##----------------------------------------------------------------------
### BATCH_simulation-scenario2.R ends here
