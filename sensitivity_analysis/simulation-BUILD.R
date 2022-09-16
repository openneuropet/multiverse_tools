### simulation-BUILD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 16 2022 (09:44) 
## Version: 
## Last-Updated: sep 16 2022 (16:43) 
##           By: Brice Ozenne
##     Update #: 23
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
library(ggplot2)
library(gridExtra)
library(mvtnorm)
library(data.table)
library(Matrix)
library(LMMstar)
library(pbapply)

if(system("whoami",intern=TRUE) == c("unicph\\hpl802")){
    path <- "x:/pipeline/"
}else if(system("whoami",intern=TRUE) == c("hpl802")){
    path <- ""
}

options(width = 120)

## * Load results
dt.sim1 <- as.data.table(readRDS(file.path(path,"results","sim1.rds")))
dt.sim2 <- as.data.table(readRDS(file.path(path,"results","sim2.rds")))
dt.sim3 <- as.data.table(readRDS(file.path(path,"results","sim3.rds")))

dt.sim <- rbind(cbind(scenario = "scenario 1", dt.sim1),
                cbind(scenario = "scenario 2", dt.sim2),
                cbind(scenario = "scenario 3", dt.sim3))


## * Process results
dt.sim[, n.char := factor(n/2, unique(n/2))]
dt.sim[, beta.char := ifelse(beta==0,paste0("Null hypothesis (\u03B2=",beta,")"),paste0("Alternative hypothesis (\u03B2=",beta,")"))]
dt.sim[, beta.char := factor(beta.char,levels=unique(beta.char))]
dt.sim[, type.char := factor(type,levels=c("average","fixse","gls","proportion"),c("pool (average)","pool (se)","pool (gls)","proportion"))]
dt.sim[, target := factor(type,levels=c("average","fixse","gls","proportion"),c("pool","pool","pool","proportion"))]
dt.sim[, GS := beta]
dt.sim[type=="proportion", GS := as.numeric(beta>0)]
dtS.sim <- dt.sim[,.(rep = .N, average = mean(estimate), sd = sd(estimate), average.se = mean(se), bias = mean(estimate)-GS, power = mean(p.value<=0.05)),
                    by = c("scenario","beta","n","type","n.char","beta.char","type.char","target","GS")]



## * display results
## bias
gg.bias <- ggplot(dtS.sim, aes(x = n.char, y = average, linetype = target, color = type.char, group = type.char)) 
gg.bias <- gg.bias + geom_hline(aes(yintercept = beta), color = "brown")
gg.bias <- gg.bias + geom_point(size = 2) + geom_line(size = 1) + facet_grid(scenario~beta.char) 
gg.bias <- gg.bias + labs(x = "Sample size (per group)", y = "Estimate", color = "") + guides(linetype = "none")
gg.bias

## se
gg.se <- ggplot(dtS.sim, aes(x = n.char, y = sd, linetype = target, color = type.char, group = type.char)) 
gg.se <- gg.se + geom_point(size = 2) + geom_line(size = 1) + facet_grid(scenario~beta.char) 
gg.se <- gg.se + labs(x = "Sample size (per group)", y = "Standard deviation", color = "") + guides(linetype = "none")
gg.se

## cali se
gg.cali <- ggplot(dtS.sim, aes(x = sd, y = average.se, linetype = target, color = type.char, group = type.char)) 
gg.cali <- gg.cali + geom_abline(intercept = 0, slope = 1, color = "brown")
gg.cali <- gg.cali + geom_point(size = 2) + geom_line(size = 1) + facet_grid(scenario~beta.char) 
gg.cali <- gg.cali + labs(x = "Sample size (per group)", y = "Power", color = "") + guides(linetype = "none")
gg.cali

## relative efficiency
table.reff <- dcast(dtS.sim[type!="proportion",.(scenario,beta.char,n.char,type.char,sd)], scenario+beta.char+n.char~type.char, value.var = "sd")
table.reff[,c("pool (se)", "pool (gls)") := .(.SD[["pool (se)"]]/.SD[["pool (average)"]],.SD[["pool (gls)"]]/.SD[["pool (average)"]])]
table.reff[, scenario := paste(scenario,beta.char,sep=": ")]
table.reff[, beta.char := NULL]
table.reff[duplicated(scenario), scenario := ""]
setnames(table.reff, old = c("n.char","pool (average)","pool (se)", "pool (gls)"), new = c("n","sd(pool-av)", "sd(pool-se)/sd(pool-av)","sd(pool-gls)/sd(pool-av)"))

ggtable.reff <- table.reff
ggtable.reff[, scenario := gsub("Null hypothesis | Alternative hypothesis |\\(|\\)","",scenario)]
ggtable.reff[["sd(pool-av)"]] <- round(ggtable.reff[["sd(pool-av)"]],3)
ggtable.reff[["sd(pool-se)/sd(pool-av)"]] <- round(ggtable.reff[["sd(pool-se)/sd(pool-av)"]],3)
ggtable.reff[["sd(pool-gls)/sd(pool-av)"]] <- round(ggtable.reff[["sd(pool-gls)/sd(pool-av)"]],3)


## power
gg.power <- ggplot(dtS.sim, aes(x = n.char, y = power, linetype = target, color = type.char, group = type.char)) 
gg.power <- gg.power + geom_hline(yintercept = 0.05, color = "brown")
gg.power <- gg.power + geom_point(size = 2) + geom_line(size = 1) + facet_grid(scenario~beta.char) 
gg.power <- gg.power + labs(x = "Sample size (per group)", y = "Power", color = "") + guides(linetype = "none")
gg.power


## proportion
dcast(dtS.sim[beta==0 & type == "proportion", .(scenario,n,average)], scenario ~ n)
##      scenario         20         50        100        200        400
## 1: scenario 1 0.08160396 0.07660466 0.06384783 0.06483713 0.06068949
## 2: scenario 2 0.06085783 0.06050266 0.05937698 0.05941852 0.06520626
## 3: scenario 3 0.05254837 0.06156795 0.05841988 0.06206031 0.05113468

dcast(dtS.sim[beta==0 & type == "average", .(scenario,n,power)], scenario ~ n)
##      scenario    20    50   100   200   400
## 1: scenario 1 0.084 0.068 0.040 0.052 0.036
## 2: scenario 2 0.064 0.060 0.044 0.068 0.044
## 3: scenario 3 0.048 0.052 0.068 0.076 0.044
dcast(dtS.sim[beta==0 & type == "fixse", .(scenario,n,power)], scenario ~ n)


## * Export

## bias
png(file.path(path,"figures","simulation-fig-bias.png"))
gg.bias
dev.off()

pdf(file.path(path,"figures","simulation-fig-bias.pdf"))
gg.bias
dev.off()

## standard error
png(file.path(path,"figures","simulation-fig-dispersion.png"))
gg.se
dev.off()

pdf(file.path(path,"figures","simulation-fig-dispersion.pdf"))
gg.se
dev.off()

## calibration se
png(file.path(path,"figures","simulation-fig-calibrationSE.png"))
gg.cali
dev.off()

pdf(file.path(path,"figures","simulation-fig-calibrationSE.pdf"))
gg.cali
dev.off()

# relative efficiency
png(file.path(path,"figures","simulation-table-reff.png"), width = 550, height = 700)
grid.arrange(tableGrob(as.data.frame(ggtable.reff), rows = NULL))
dev.off()

pdf(file.path(path,"figures","simulation-table-reff.pdf"), width = 550, height = 700)
grid.arrange(tableGrob(as.data.frame(ggtable.reff), rows = NULL))
dev.off()

## power
png(file.path(path,"figures","simulation-figure-power.png"))
gg.power
dev.off()

pdf(file.path(path,"figures","simulation-figure-power.pdf"))
gg.power
dev.off()


##----------------------------------------------------------------------
### simulation-BUILD.R ends here
