### simulation-FCT.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 15 2022 (16:48) 
## Version: 
## Last-Updated: sep 16 2022 (11:07) 
##           By: Brice Ozenne
##     Update #: 21
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * simData
simData <- function(n.obs, sigma.pipe, beta){
    require(mvtnorm)
    require(data.table)
    
    n.pipe <- NCOL(sigma.pipe)
    
    ## design matrix
    X <- factor(c(rep("Male",n.obs),rep("Female",n.obs)))

    ## generate Y
    dtW.data <- data.table(id = 1:(2*n.obs),
                          X = X,
                          X.num = as.numeric(X)-1,
                          Y = rnorm(2*n.obs, mean = beta*(as.numeric(X)-1), sd = 1))

    ## generate observed Y
    noise.pipeline <- mvtnorm::rmvnorm(n=2*n.obs, sigma = sigma.pipe)

    for(iP in 1:n.pipe){
        dtW.data[, c(paste0("Ypip.",iP)) := Y + noise.pipeline[,iP]]
    }

    ## move to long format
    dtL.data <- melt(dtW.data, id.vars = c("id","X","X.num"),
                     measure.vars = paste0("Ypip.",1:n.pipe),
                     value.name = "Ypip", variable.name = "pipeline")

    ## export
    return(list(wide = dtW.data,
                long = dtL.data,
                beta = beta))
}

## * analyzeData
analyzeData <- function(data, proportion){
    require(LMMstar)

    e.mlmm <- do.call(mlmm, list(Ypip ~ X, repetition = ~1|id, data = data$long, by = "pipeline", effects = "XMale=0", df = TRUE))

    out <- rbind(average = model.tables(e.mlmm, method = "average"),
                 fixse = model.tables(e.mlmm, method = "pool.fixse"),
                 gls = model.tables(e.mlmm, method = "pool.gls"))
    if(proportion){
        out <- rbind(out, proportion = proportion(e.mlmm, method = "single-step", n.sample = 0))
    }
    
    return(cbind(beta = data$beta, n = NROW(data$wide), type = rownames(out), out))
}

##----------------------------------------------------------------------
### simulation-FCT.R ends here
