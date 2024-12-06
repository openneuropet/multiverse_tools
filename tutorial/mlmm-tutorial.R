### mlmm-tutorial.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 25 2024 (16:23) 
## Version: 
## Last-Updated: okt 25 2024 (17:00) 
##           By: Brice Ozenne
##     Update #: 14
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(mvtnorm)
library(Matrix)
library(LMMstar)
library(data.table)
library(ggplot2)


## * Simulate data

## ** parametrisation
n.obs <- 500 ## number of subjects in each group
beta <- 0.5 ## group difference
noise.Sigma <- as.matrix(bdiag((diag(0.05,15,15) + 0.95)*2.5, matrix(0.25), matrix(5), matrix(7.5), matrix(10), matrix(15))) ## noise distribution
n.pipe <- NCOL(noise.Sigma)

## ** simulate
set.seed(1)
X <- factor(c(rep("Male",n.obs),rep("Female",n.obs)))
dtW.data <- data.table(id = 1:(2*n.obs),
                       X = X,
                       X.num = as.numeric(X)-1,
                       Y = rnorm(2*n.obs, mean = beta*(as.numeric(X)-1), sd = 1)) ## Y is not observed
noise.pipeline <- rmvnorm(n=2*n.obs, sigma = noise.Sigma)
dtW.data <- cbind(dtW.data, matrix(dtW.data$Y, ncol = n.pipe, nrow = 2*n.obs, byrow = FALSE) + noise.pipeline)
names(dtW.data) <- c("id","X","X.num","Y",paste0("Ypip.",1:n.pipe))

## move to the long format (observed data)
dtL.data <- melt(dtW.data, id.vars = c("id","X","X.num"),
                 measure.vars = paste0("Ypip.",1:n.pipe),
                 value.name = "Ypip", variable.name = "pipeline")

## * Analysis

## ** fit pipeline specific models & get pipeline specific estimates
## run multiple linear model, one per pipeline (argument by), and extract the test w.r.t. the gender effect (argument effects).
## the argument repetition enables to keep track of the structure of the data to be able to perform statistical inference accounting for within-subject correlation.
e.mlmm <- mlmm(Ypip ~ X, data = dtL.data,
               repetition = ~1|id,
               by = "pipeline",
               effects = "XMale=0")

summary(e.mlmm)
## 	Linear Mixed Models stratified according to "pipeline"  [TOFIX: mixed should be removed]

## Statistical inference for XMale  [TOFIX: wrong se is displayed]

##            estimate    se  df  lower upper p.value    
##    Ypip.1     0.147 0.261 198 -0.368 0.663 0.57352    
##    Ypip.2     0.237 0.263 198 -0.281 0.755 0.36763    
##    Ypip.3     0.198 0.265 198 -0.324 0.719 0.45610    
##    Ypip.4     0.238 0.264 198 -0.284 0.759 0.36971    
##    Ypip.5     0.205 0.268 198 -0.323 0.732 0.44507    
##    Ypip.6     0.211 0.269 198 -0.319  0.74 0.43379    
##    Ypip.7     0.097  0.26 198 -0.416 0.611 0.70846    
##    Ypip.8     0.104  0.27 198 -0.429 0.637 0.70109    
##    Ypip.9     0.185 0.265 198 -0.338 0.709 0.48578    
##    Ypip.10    0.293  0.26 198 -0.219 0.805 0.26075    
##    Ypip.11    0.165  0.27 198 -0.367 0.696 0.54203    
##    Ypip.12    0.214 0.272 198 -0.322  0.75 0.43274    
##    Ypip.13    0.068 0.265 198 -0.455 0.592 0.79750    
##    Ypip.14     0.16 0.264 198 -0.361 0.681 0.54626    
##    Ypip.15    0.218 0.266 198 -0.306 0.743 0.41271    
##    Ypip.16    0.685 0.151 198  0.387 0.983 < 1e-04 ***
##    Ypip.17    1.203  0.37 198  0.474 1.933 0.00134  **
##    Ypip.18     0.69 0.433 198 -0.164 1.544 0.11258    
##    Ypip.19    0.251 0.457 198 -0.649 1.152 0.58248    
##    Ypip.20    1.302  0.59 198  0.138 2.466 0.02856   *
##    -------------------------------------------------- 
##   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1.
##   Model-based standard errors are derived from the observed information (column se). 
##   Degrees of freedom were computed using a Satterthwaite approximation (column df). 

summary(e.mlmm$model[[1]])

## ** [optional] evaluate each pipeline in term of s.e. and how they correlate
e.mlmm_vcov <- vcov(e.mlmm) 

sqrt(diag(e.mlmm_vcov)) ## standard error w.r.t. each pipeline
##     Ypip.1     Ypip.2     Ypip.3     Ypip.4     Ypip.5     Ypip.6     Ypip.7 
## 0.11743723 0.11875811 0.11801756 0.11859325 0.11821299 0.11747522 0.11771105 
##     Ypip.8     Ypip.9    Ypip.10    Ypip.11    Ypip.12    Ypip.13    Ypip.14 
## 0.11730527 0.11969145 0.11711779 0.11763687 0.11889764 0.11666183 0.11840204 
##    Ypip.15    Ypip.16    Ypip.17    Ypip.18    Ypip.19    Ypip.20 
## 0.11706961 0.07212763 0.15997801 0.18662571 0.20905938 0.25548422 

range(cov2cor(e.mlmm_vcov)) ## correlation between pipeline estimates
plot(e.mlmm, "heat") ## automatic graphical display

## ** global effect estimator
confint(e.mlmm, method = "average", columns = c("estimate","se","df","lower","upper","p.value"))
confint(e.mlmm, method = "pool.fixse", columns = c("estimate","se","df","lower","upper","p.value"))
confint(e.mlmm, method = "pool.gls", columns = c("estimate","se","df","lower","upper","p.value"))
confint(e.mlmm, method = "pool.gls1", columns = c("estimate","se","df","lower","upper","p.value"))

## ** forest plot
plot(e.mlmm) ## without global effect estimators

plot(e.mlmm, method = c("none","average","pool.se","pool.gls","pool.gls1"))

## make it more pretty
shape.forest <- c(rep(19,20),c(15,17,19,8))
color.forest <- c(rep("black",20), palette.colors()[2:5])
plot(e.mlmm, method = c("none","average","pool.se","pool.gls","pool.gls1"), 
     add.args = list(size.estimate = 2.5, size.ci = 1, width.ci = 0.5, shape = shape.forest, color = color.forest))

## ** proportion estimator
confint(e.mlmm, method = "proportion", columns = c("estimate","se","df","lower","upper","p.value"))
##----------------------------------------------------------------------
### mlmm-tutorial.R ends here
