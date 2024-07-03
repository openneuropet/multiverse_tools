### TEXT-simulation.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  3 2024 (11:55) 
## Version: 
## Last-Updated: jul  3 2024 (15:34) 
##           By: Brice Ozenne
##     Update #: 24
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(Matrix)
dtS.sim <- readRDS(file.path("Results","simulation-summary-scenario.rds"))

## * Weights

## ** scenario 1
Sigma1 <- as.matrix(bdiag(diag(0.05,15,15) + 0.95,matrix(1),matrix(1),matrix(1),matrix(1),matrix(1)))
Sigmatot1 <- Sigma+diag(1,NROW(Sigma1),NROW(Sigma1))

cov2cor(Sigmatot1)[1,2]
## [1] 0.475

## average
weight1.average <- 1/NROW(Sigmatot1)
## pool-se
weight1.poolse <- (1/diag(Sigmatot1))/sum((1/diag(Sigmatot1)))
## gls
weight1.gls <- rowSums(solve(Sigmatot1))/sum(solve(Sigmatot1))

100*weight1.average
## [1] 5
100*weight1.poolse
## [1] 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5
round(100*weight1.gls[1:15],2)
## [1] 1.88 1.88 1.88 1.88 1.88 1.88 1.88 1.88 1.88 1.88 1.88 1.88 1.88 1.88 1.88
round(100*weight1.gls[16:20],2)
## [1] 14.37 14.37 14.37 14.37 14.37

## ** scenario 2
Sigma2 <- diag(c(2.5,0.25,5,7.5,10,15))
Sigmatot2 <- Sigma+diag(1,NROW(Sigma2),NROW(Sigma2))

sqrt(diag(Sigmatot2))/min(sqrt(diag(Sigmatot2)))
## [1] 1.673320 1.000000 2.190890 2.607681 2.966479 3.577709

## average
weight2.average <- 1/NROW(Sigmatot2)
## pool-se
weight2.poolse <- (1/diag(Sigmatot2))/sum((1/diag(Sigmatot2)))
## gls
weight2.gls <- rowSums(solve(Sigmatot2))/sum(solve(Sigmatot2))

round(100*weight2.average,2)
## [1] 16.67
round(100*weight2.poolse,2)
## [1] 18.75 52.51 10.94  7.72  5.97  4.10
weight2.gls - weight2.poolse
## [1] 0 0 0 0 0 0

## ** scenario 3
Sigma3 <- as.matrix(bdiag((diag(0.05,15,15) + 0.95)*2.5, matrix(0.25), matrix(5), matrix(7.5), matrix(10), matrix(15)))
cov2cor(Sigma3)[1,2]
## [1] 0.95
Sigmatot3 <- Sigma3+diag(1,NROW(Sigma3),NROW(Sigma3))
cov2cor(Sigmatot3)[1,2]
## [1] 0.6785714

## average
weight3.average <- 1/NROW(Sigmatot3)
## pool-se
weight3.poolse <- (1/diag(Sigmatot3))/sum((1/diag(Sigmatot3)))
## gls
weight3.gls <- rowSums(solve(Sigmatot3))/sum(solve(Sigmatot3))

round(100*weight3.average,2)
## [1] 5
round(100*weight3.poolse,2)
##  [1]  5.17  5.17  5.17  5.17  5.17  5.17  5.17  5.17  5.17  5.17  5.17  5.17
## [13]  5.17  5.17  5.17 14.48  3.02  2.13  1.65  1.13
round(100*weight3.gls,2)
##  [1]  1.65  1.65  1.65  1.65  1.65  1.65  1.65  1.65  1.65  1.65  1.65  1.65
## [13]  1.65  1.65  1.65 48.61 10.13  7.15  5.52  3.80


## * Estimate and bias

## ** pool
ggPool.bias <- ggplot(dtS.sim[target=="pool" & type != "robust"], aes(x = n.char, y = average, linetype = target, group = type.char)) 
ggPool.bias <- ggPool.bias + geom_hline(aes(yintercept = beta), color = "brown") + facet_grid(beta.char~scenario, scales = "free") 
ggPool.bias <- ggPool.bias + geom_point(size = 2, aes(shape = type.char, color = type.char)) + geom_line(linewidth = 1, aes(color = type.char))
ggPool.bias <- ggPool.bias + labs(x = "Sample size (per group)", y = "Estimate", shape = "", color = "") + guides(linetype = "none")
ggPool.bias

dtS.sim[beta == 0 & target == "pool", .(medianBias = median(bias), maxBias = max(abs(bias))),by="type"]
dtS.sim[target == "pool" & type != "robust"][max(abs(bias)) == abs(bias),.(scenario,beta,n,type,GS,average,bias)]
##      scenario  beta     n   type    GS   average        bias
##        <char> <num> <int> <char> <num>     <num>       <num>
## 1: scenario 3   0.5    20    gls   0.5 0.4850872 -0.01491278

unique(dtS.sim$n)
unique(dtS.sim$n.char)

## ** proportion
dtS.sim[beta==0 & target == "proportion", .(min = round(100*min(average),2), max = round(100*max(average),2)), by = "type"]
##             type   min   max
##           <char> <num> <num>
## 1: proportion.np  0.50  1.71
## 2:  proportion.p  5.35  6.79

dcast(dtS.sim[beta==0 & target == "proportion", .(round(100*average,2)), by = c("type","n","scenario")], type+scenario~n)
##             type   scenario    20    50   100   200   400  1000
##           <char>     <char> <num> <num> <num> <num> <num> <num>
## 1: proportion.np scenario 1  1.71  0.94  0.64  0.60  1.16  0.77
## 2: proportion.np scenario 2  1.23  0.93  0.97  0.88  1.02  0.82
## 3: proportion.np scenario 3  1.41  0.76  0.52  0.50  0.86  0.64
## 4:  proportion.p scenario 1  6.79  6.48  6.31  5.82  6.74  6.14
## 5:  proportion.p scenario 2  6.05  6.09  6.32  6.17  6.37  6.44
## 6:  proportion.p scenario 3  5.90  5.75  5.76  5.35  5.96  5.58

dcast(dtS.sim[beta>0 & target == "proportion", .(round(100*average,2)), by = c("type","n","scenario")], type+scenario~n)
##             type   scenario    20    50   100   200   400  1000
##           <char>     <char> <num> <num> <num> <num> <num> <num>
## 1: proportion.np scenario 1  4.18  8.70 19.10 44.00 81.14 99.70
## 2: proportion.np scenario 2  2.82  5.62 10.30 20.78 36.53 64.92
## 3: proportion.np scenario 3  2.43  4.47  8.82 19.58 44.25 85.08
## 4:  proportion.p scenario 1 10.17 16.52 27.49 45.56 73.99 98.12
## 5:  proportion.p scenario 2  8.31 12.24 16.87 26.29 40.41 63.57
## 6:  proportion.p scenario 3  7.58 10.88 16.93 26.70 46.22 79.19

## * Efficiency

## ** pool
dtS.sim[beta == 0 & scenario == "scenario 1" & type == "average",sd]
## [1] 0.56031872 0.35681882 0.24965752 0.17576386 0.12422384 0.07879121
dtS.sim[beta == 0 & scenario == "scenario 1" & type == "fixse",sd]
## [1] 0.5583603 0.3563523 0.2496001 0.1756760 0.1242233 0.0787807
round(100*(dtS.sim[beta == 0 & scenario == "scenario 1" & type == "gls",sd]/dtS.sim[beta == 0 & scenario == "scenario 1" & type == "average",sd]-1),2)
## [1] 155.69  10.97  -3.66  -8.52 -10.86 -11.29

dtS.sim[beta == 0 & scenario == "scenario 2" & type == "fixse",sd]
## [1] 0.51719850 0.32415765 0.23323722 0.16340775 0.11530477 0.07273707
dtS.sim[beta == 0 & scenario == "scenario 2" & type == "gls",sd]
## [1] 0.58035628 0.32740744 0.22736689 0.15726154 0.11053911 0.06947867
round(100*(dtS.sim[beta == 0 & scenario == "scenario 2" & type == "gls",sd]/dtS.sim[beta == 0 & scenario == "scenario 2" & type == "average",sd]-1),2)
## [1]  -9.93 -19.60 -22.74 -23.93 -24.13 -24.86

dtS.sim[beta == 0 & scenario == "scenario 3" & type == "average",sd]
## [1] 0.70034058 0.44515007 0.31235500 0.22121723 0.15596532 0.09860849
dtS.sim[beta == 0 & scenario == "scenario 3" & type == "fixse",sd]
## [1] 0.69930202 0.44702021 0.31377783 0.22215794 0.15633837 0.09883643
round(100*(dtS.sim[beta == 0 & scenario == "scenario 3" & type == "gls",sd]/dtS.sim[beta == 0 & scenario == "scenario 3" & type == "average",sd]-1),2)
## [1] 105.00  -9.53 -21.60 -26.09 -27.69 -27.75
round(100*(dtS.sim[beta == 0 & scenario == "scenario 3" & type == "gls1",sd]/dtS.sim[beta == 0 & scenario == "scenario 3" & type == "average",sd]-1),2)
## [1]  -0.05 -14.09 -21.75 -26.09 -27.69 -27.75

## ** proportion
eff.prop <- dtS.sim[beta >0 & type == "proportion.p",sd]/dtS.sim[beta >0 & type == "proportion.np",sd]-1
round(100*eff.prop,2)
##  [1] -11.98 -19.54 -29.93 -36.75 -31.03  22.96  -9.89 -17.00 -22.04 -29.30
## [11] -28.10 -35.98  -5.67 -10.15 -21.51 -30.67 -36.94 -15.29
round(100*range(eff.prop),2)
## [1] -36.94  22.96
round(100*range(eff.prop[eff.prop<0]),2)
## [1] -36.94  -5.67

## * Type 1 error
dcast(dtS.sim[beta == 0 & target == "pool",], scenario+type~n , value.var = "power")
##       scenario    type     20     50    100    200    400   1000
##         <char>  <char>  <num>  <num>  <num>  <num>  <num>  <num>
##  1: scenario 1 average 0.0631 0.0581 0.0524 0.0464 0.0518 0.0512
##  2: scenario 1   fixse 0.0689 0.0592 0.0531 0.0471 0.0521 0.0516
##  3: scenario 1     gls 0.7702 0.2484 0.1164 0.0788 0.0608 0.0586
##  4: scenario 1    gls1 0.1613 0.1471 0.1070 0.0785 0.0608 0.0586

##  6: scenario 2 average 0.0574 0.0523 0.0563 0.0536 0.0542 0.0536
##  7: scenario 2   fixse 0.0736 0.0584 0.0593 0.0541 0.0535 0.0495
##  8: scenario 2     gls 0.1720 0.0833 0.0703 0.0565 0.0555 0.0512
##  9: scenario 2    gls1 0.1620 0.0819 0.0703 0.0565 0.0555 0.0512

## 11: scenario 3 average 0.0625 0.0555 0.0524 0.0497 0.0522 0.0505
## 12: scenario 3   fixse 0.0656 0.0572 0.0540 0.0496 0.0506 0.0496
## 13: scenario 3     gls 0.7660 0.2477 0.1203 0.0790 0.0618 0.0614
## 14: scenario 3    gls1 0.2202 0.2048 0.1199 0.0790 0.0618 0.0614

max(dtS.sim[beta == 0 & type == "average",power])
## [1] 0.0631
max(dtS.sim[beta == 0 & type == "fixse",power])
## [1] 0.0736
max(dtS.sim[beta == 0 & type == "gls",power])
## [1] 0.7702
max(dtS.sim[beta == 0 & type == "gls1",power])
## [1] 0.2202
##----------------------------------------------------------------------
### TEXT-simulation.R ends here
