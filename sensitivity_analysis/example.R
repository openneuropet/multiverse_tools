### example.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  9 2022 (13:24) 
## Version: 
## Last-Updated: sep  9 2022 (13:24) 
##           By: Brice Ozenne
##     Update #: 1
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
library(LMMstar) ## require development version >= 0.8.0
## can be installed using remotes::install_github("bozenne/LMMstar")

## * Simulate data
## ** settings
n <- 100
J <- 5
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

## ** extract key quantities
## design matrix
X <- model.matrix(ls.lm[[1]])
as.data.table(X)

## model parameters
ls.Theta <- lapply(ls.lm, function(iLM){c(coef(iLM),sigma = sigma(iLM))})
M.Theta <- do.call(rbind,ls.Theta)
M.Theta

## target parameter
psi <- (M.Theta %*% c(0,1,0))[,1]
psi

## influence function
M.iid <- do.call(cbind,lapply(ls.lm, function(iLM){iid(iLM)[,2]}))*sqrt(2*n)
as.data.table(M.iid)

## (manually)
score <- X[1,,drop=FALSE] * residuals(ls.lm[[1]])[1] / M.Theta[1,3]
information_indiv <- solve(crossprod(X)) * (M.Theta[1,3]) * sqrt(2*n)
score %*% information_indiv %*% c(0,1)

## ** combine across models
## variance-covariance
Sigma_psi <- 1/(2*n)*crossprod(M.iid)
Sigma_psi
## (manually)
sum(sapply(1:(2*n), function(iObs){residuals(ls.lm[[1]])[iObs]*residuals(ls.lm[[2]])[iObs] * (X[iObs,,drop=FALSE] %*% solve(crossprod(X)))^2 %*% c(0,1)}))

## correlation
psi.se <- sqrt(diag(Sigma_psi))
Sigma_t <- Sigma_psi / tcrossprod(psi.se)
Sigma_t

## Waæd statistic
waldStat <- psi/psi.se
waldStat

## * Testing hypotheses
## ** No effects across all pipelines
## max statistic
waldStatMax <- max(abs(waldStat))
waldStatMax

## adjusted p-value
1 - pmvnorm(lower = -rep(waldStatMax,J), upper = rep(waldStatMax,J), mean = rep(0,J), corr = Sigma_t)

## critical quantile
t_c <- qmvnorm(0.95, mu = rep(0,J), sigma = Sigma_t, tail = "both.tails")$quantile
t_c

## confidence intervals
cbind(lower = psi - psi.se * t_c, estimate = psi, upper = psi + psi.se * t_c)

## LMMstar version
e.mlmm <- mlmm(value~group, repetition = ~1|id, data = dtL, df = FALSE, robust = TRUE,
                by = "pipeline", effects = "groupG2=0")
summary(e.mlmm, method = "single-step")


## ** At least one pipeline with no effect
max(2*(1-pnorm(abs(waldStat))))

## ** Percentage of significant results
## non parametric
mean(waldStat > t_c)

## parametric
1-mean(sapply(waldStat, function(iStat){pnorm(t_c - iStat)-pnorm(-t_c - iStat)}))

## bootstrap ci
Upipe <- paste0("V",1:J)
calcPC <- function(data, statistic, ...){ ## mydata <- dtL
    if(is.null(statistic)){newdata <- data}else{newdata <- data[statistic]}
    iLM <- lapply(Upipe, function(iPipe){lm(value~group, data = newdata[pipeline==iPipe])})
    iBeta <- sapply(iLM, function(iiLM){coef(iiLM)[2]})
    iSigma <- crossprod(do.call(cbind,lapply(iLM, function(iiLM){iid(iiLM)[,2]})))
    iBeta.se <- sqrt(diag(iSigma))

    iT_c <- qmvnorm(0.95, mu = rep(0,J), sigma = cov2cor(iSigma),  tail = "both.tails")$quantile
    iWald <- iBeta/iBeta.se
    out <- 1-mean(sapply(iWald, function(iiStat){pnorm(iT_c - iiStat)-pnorm(-iT_c - iiStat)}))
    return(out)
}
## calcPC(dtL, statistic  = NULL)
library(boot)
e.boot <- boot(data = dtL, statistic = calcPC, R = 500)
e.boot
##      original        bias     std. error
## t1* 0.4129262 -1.022495e-05 0.0001583698
## plot(e.boot)

## delta-method
e.mlmmML <- mlmm(value~group, repetition = ~1|id, data = dtL, df = FALSE, 
                  by = "pipeline", effects = "all", method = "ML", transform.sigma = "none")

vec.coefML <- coef(e.mlmmML)
M.vcovML <- vcov(e.mlmmML)

index.alpha <- 1 + (1:5-1)*3
index.psi <- 2 + (1:5-1)*3
index.sigma <- 3 + (1:5-1)*3

## range(vec.coefML[index.psi] - psi)

vec.sigmaML <- vec.coefML[index.sigma]
## range(vec.sigmaML - sapply(e.mlmm$model, coef, effects = "variance"))


de_dtheta <- 0*vec.coefML
de_dtheta[index.sigma] <- -psi/(J*psi.se*vec.sigmaML) * sapply(waldStat, function(iStat){dnorm(t_c - iStat)-dnorm(-t_c - iStat)})
de_dtheta[index.psi] <- 1/(J*psi.se) * sapply(waldStat, function(iStat){dnorm(t_c - iStat)-dnorm(-t_c - iStat)})

sqrt(rbind(de_dtheta) %*% M.vcovML %*% cbind(de_dtheta))
sqrt(rbind(de_dtheta[index.psi]) %*% M.vcovML[index.psi,index.psi] %*% cbind(de_dtheta[index.psi]))




calcPC2 <- function(p, ...){ ## p <- vec.coef 

    return(pnorm(t_c-p/psi.se))
    ## return(1-mean(sapply(p/psi.se, function(iiStat){pnorm(t_c - iiStat)-pnorm(-t_c - iiStat)})))

}
## calcPC2(psi)
sqrt(diag(numDeriv::jacobian(calcPC2,psi) %*% Sigma_psi %*% t(numDeriv::jacobian(calcPC2,psi))))
 ## [1] 0.017489802 0.017311625 0.017682413 0.017766174 0.017606586 0.001731572
 ## [7] 0.002936276 0.025185324 0.142186331 0.013247053

numDeriv::jacobian(calcPC2,psi)
dnorm(t_c - psi/psi.se)
## sqrt(dnorm(t_c - psi/psi.se)^2)

calcPC <- function(data, statistic, ...){ ## mydata <- dtL
    if(is.null(statistic)){newdata <- data}else{newdata <- data[statistic]}
    iLM <- lapply(Upipe, function(iPipe){lm(value~group, data = newdata[pipeline==iPipe])})
    iBeta <- sapply(iLM, function(iiLM){coef(iiLM)[2]})
    ## out <- 1-mean(sapply(iBeta/psi.se, function(iiStat){pnorm(t_c - iiStat)-pnorm(-t_c - iiStat)}))
    out <- pnorm(t_c-iBeta/psi.se)
    ## out <- pnorm(iBeta)
    return(out)
}
library(boot)
e.boot <- boot(data = dtL, statistic = calcPC, R = 500)
e.boot

## Bootstrap Statistics :
##         original      bias    std. error
## t1* 0.0004861530 0.009740197  0.03819484
## t2* 0.0008615514 0.010796368  0.04021281
## t3* 0.0093727089 0.040218325  0.09889667
## t4* 0.0754387454 0.074525983  0.18652505
## t5* 0.0045322270 0.025820022  0.06240041

## * Visualization
gg.heatmap <- ggplot(reshape2::melt(Sigma_t)) + geom_tile(aes(x=Var1,y=Var2,fill=value))
gg.heatmap <- gg.heatmap + scale_fill_gradient(limits = c(0,1.00001), low = "white", high = "red")
gg.heatmap <- gg.heatmap + labs(x="Pipeline",y="Pipeline", fill = "correlation")
gg.heatmap

plot(e.mlmm)

## * Pooling
## ** naive
mean(psi)
summary(e.mlmm, method = "average")

## ** variance pooling
psi.w <- 1/psi.se^2
weighted.mean(psi, psi.w)
sum(psi*psi.w)/sum(psi.w)

e.mlmm <- mlmm(value~group, repetition = ~1|id, data = dtL, df = FALSE, robust = TRUE,
                by = "pipeline", effects = "groupG2=0")

summary(e.mlmm, method = "pool.fixse")
summary(e.mlmm, method = "pool.se")

e.mlmmML <- mlmm(value~group, repetition = ~1|id, data = dtL, df = FALSE, robust = TRUE,
                 by = "pipeline", effects = "groupG2=0", method.fit = "ML")
summary(e.mlmmML, method = "pool.fixse")
summary(e.mlmmML, method = "pool.se")

## ** GLS
e.glsUN <- lmm(value ~ pipeline + group, repetition = ~pipeline|id, data = dtL,
               structure = "UN")
model.tables(e.glsUN)["groupG2",]

e.glsCS <- lmm(value ~ pipeline + group, repetition = ~pipeline|id, data = dtL,
               structure = "CS")
model.tables(e.glsCS)["groupG2",]

## ** PCA
eigen_psi <- eigen(Sigma_psi)
kappa <- sweep(eigen_psi$vectors, FUN = "/", MARGIN = 2, STATS = colSums(eigen_psi$vectors))
omega <- eigen_psi$values/colSums(eigen_psi$vectors)^2
weighted.mean(psi %*% kappa, w = 1/omega)

summary(e.mlmmML, method = "pool.pca")

## * Systematic differences
psi[-1] - psi[1]

eD.mlmm <- mlmm(value~group, repetition = ~1|id, data = dtL, df = FALSE,
                by = "pipeline", effects = "groupG2=0", contrast.rbind = "Dunnett")
summary(eD.mlmm)
summary(eD.mlmm, method = "single-step")

eT.mlmm <- mlmm(value~group, repetition = ~1|id, data = dtL, df = FALSE,
                by = "pipeline", effects = "groupG2=0", contrast.rbind = "Tukey")
summary(eT.mlmm)
summary(eT.mlmm, method = "single-step")
summary(eT.mlmm, method = "Westfall")

##----------------------------------------------------------------------
### example.R ends here
