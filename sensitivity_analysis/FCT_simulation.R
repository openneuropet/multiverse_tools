### FCT_simulation.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 15 2022 (16:48) 
## Version: 
## Last-Updated: feb 24 2023 (16:36) 
##           By: Brice Ozenne
##     Update #: 68
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * simData
simData <- function(n.obs, sigma.pipe, beta, df = Inf){
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
    if(is.infinite(df)){
        noise.pipeline <- mvtnorm::rmvnorm(n=2*n.obs, sigma = sigma.pipe)
    }else{
        noise.pipeline <- mvtnorm::rmvt(n=2*n.obs, sigma = sigma.pipe, df = df)
    }
    
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
analyzeData <- function(data, proportion, print.weight = FALSE, df = TRUE){
    require(LMMstar)
    require(lme4)
    require(lmerTest)
    
    ## table(data$long[, .N, by = "id"][[2]])
    e.lmerCS <- lmer(Ypip ~ X + (1|id), data = data$long) 
    pool.ranef <- data.frame(estimate = unname(fixef(e.lmerCS)["XMale"]),
                             se = summary(e.lmerCS)$coef["XMale","Std. Error"],
                             df = summary(e.lmerCS)$coef["XMale","df"],
                             lower = confint(e.lmerCS, parm = "XMale", method = "Wald", quiet = TRUE)[1,1],
                             upper = confint(e.lmerCS, parm = "XMale", method = "Wald", quiet = TRUE)[1,2],
                             p.value = summary(e.lmerCS)$coef["XMale","Pr(>|t|)"])
    ##e.lmerUN <- lmer(Ypip ~ X + (0+pipeline|id), data = data$long, control = lmerControl(check.nobs.vs.nRE = "ignore"))

    e.mlmm <- do.call(mlmm, list(Ypip ~ X, repetition = ~1|id, data = data$long, by = "pipeline", effects = "XMale=0", df = df, trace = FALSE))
    pool.average <- confint(e.mlmm, method = "average", columns = c("estimate","se","df","lower","upper","p.value"))
    pool.fixse <- confint(e.mlmm, method = "pool.fixse", columns = c("estimate","se","df","lower","upper","p.value"))
    pool.gls <- confint(e.mlmm, method = "pool.gls", columns = c("estimate","se","df","lower","upper","p.value"))
    pool.gls1 <- confint(e.mlmm, method = "pool.gls1", columns = c("estimate","se","df","lower","upper","p.value"))
    M.iid <- do.call(cbind,lapply(e.mlmm$model, function(iM){iid(iM)[,"XMale",drop=FALSE]}))
    M.rSigmaM1 <- robustInverse(M.iid)
    pool.robust <- data.frame(estimate = sum(M.rSigmaM1 %*% coef(e.mlmm)) / sum(M.rSigmaM1), se = NA, df = NA, lower = NA, upper = NA, p.value = NA)

    if(print.weight){
        weightTable.average <- table(round(attr(pool.average,"contrast"),4))
        names(weightTable.average) <- as.character(round(as.numeric(names(weightTable.average)),3))
        weightTable.fixse <- table(round(attr(pool.fixse,"contrast"),4))
        names(weightTable.fixse) <- as.character(round(as.numeric(names(weightTable.fixse)),3))
        weightTable.gls <- table(round(attr(pool.gls,"contrast"),4))
        names(weightTable.gls) <- as.character(round(as.numeric(names(weightTable.gls)),3))
        weightTable.gls1 <- table(round(attr(pool.gls1,"contrast"),4))
        names(weightTable.gls1) <- as.character(round(as.numeric(names(weightTable.gls1)),3))
        
        cat("Weights:\n",
            "- average: ", paste(names(weightTable.average), collapse = "/"), " (",paste(weightTable.average, collapse = "/"),")\n",
            "- fixse: ", paste(names(weightTable.fixse), collapse = "/"), " (",paste(weightTable.fixse, collapse = "/"),")\n", 
            "- gls: ", paste(names(weightTable.gls), collapse = "/"), " (",paste(weightTable.gls, collapse = "/"),")\n",
            "- gls1: ", paste(names(weightTable.gls1), collapse = "/"), " (",paste(weightTable.gls1, collapse = "/"),")\n",
            sep = "")
    }

    out <- rbind(ranef = pool.ranef,
                 average = pool.average,
                 fixse = pool.fixse,
                 gls = pool.gls,
                 gls1 = pool.gls1,
                 robust = pool.robust)
    if(proportion){
        out <- rbind(out, proportion.np = c(estimate = mean(model.tables(e.mlmm, method = "single-step")$p.value<=0.05),
                                            se = NA, df = NA, lower = NA, upper = NA, p.value = NA))
        out <- rbind(out, proportion.p = proportion(e.mlmm, method = "single-step", n.sample = 0))
    }
    
    return(cbind(beta = data$beta, n = NROW(data$wide), type = rownames(out), out))
}

## * robustInverse
## Input: n times p matrix
robustInverse <- function(object, robust = TRUE){

    ## range(svd(O)$v %*% diag(svd(O)$d^2) %*% t(svd(O)$v) - crossprod(O))
    ## [1] -9.094947e-13  2.728484e-12
    if(robust){fct.center <- median}else{fct.center <- mean}
    if(robust){fct.scale <- mad}else{fct.scale <- sd}

    ## normalize
    O <- as.matrix(object)
    n.O <- NROW(O)
    p.O <- NCOL(O)
    ## substract median (column-wise)
    O.centerscale <- scale(O, center = apply(O,2,fct.center), scale = apply(O,2,fct.scale))
    ## substract mean  (column-wise)
    O.centerscale2 <- scale(O.centerscale, center = apply(O.centerscale,2,mean), scale = FALSE)
    ## compute SVD
    O.svd <- svd(O.centerscale2)
    ## subset
    svd.subset <- which(abs(O.svd$d) > 1e-10)
    ## inverse
    OM1 <- O.svd$v[,svd.subset,drop=FALSE] %*% diag(1/O.svd$d[svd.subset]^2) %*% t(O.svd$v[,svd.subset,drop=FALSE])
    return(OM1)
}
## SANITY CHECK
## range(robustInverse(iid55, robust = FALSE) - solve(cov2cor(Sigma55))/9)
## [1] -4.218847e-15  3.552714e-15

##----------------------------------------------------------------------
### FCT_simulation.R ends here
