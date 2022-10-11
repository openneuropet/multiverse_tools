### BUILD_simulation.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 16 2022 (09:44) 
## Version: 
## Last-Updated: okt 10 2022 (18:22) 
##           By: Brice Ozenne
##     Update #: 36
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## cd /projects/biostat01/people/hpl802/pipeline/
## * Path
if(system("whoami",intern=TRUE) == c("unicph\\hpl802")){
    path <- "x:/pipeline/"
}else if(system("whoami",intern=TRUE) == c("hpl802")){
    path <- ""
}
path.results <- file.path(path,"Results")
path.results1 <- file.path(path.results,"simulation-scenario1")
path.results2 <- file.path(path.results,"simulation-scenario2")
path.results3 <- file.path(path.results,"simulation-scenario3")

## * Dependencies
library(ggplot2)
library(gridExtra)
library(mvtnorm)
library(data.table)
library(Matrix)
library(pbapply)
library(xtable)

loadRes <- function(path, tempo.file = FALSE, type = NULL,
                    export.attribute = NULL, trace = TRUE){
    all.files <- list.files(path)
    file.tempo <- grep("(tempo)",all.files,value = TRUE)
    file.final <- setdiff(all.files, file.tempo)

    if(tempo.file){
        file.read <- file.tempo
    }else{
        file.read <- file.final
    }
    if(!is.null(type)){
        file.read <- grep(pattern=type,x=file.read,value=TRUE)
    }

    n.file <- length(file.read)

    myApply <- switch(as.character(as.logical(trace)),
                      "TRUE" = pbapply::pblapply,
                      "FALSE" = lapply)

    ls.out <- do.call(myApply, args = list(X = 1:n.file, FUN = function(iFile){
        iRead <- try(readRDS(file = file.path(path,file.read[iFile])))
        if(inherits(iRead,"try-error")){
            return(NULL)
        }else{
            iOut <- cbind(data.table::as.data.table(iRead),
                          file = file.read[iFile])
        return(iOut)
        }
    }))
    out <- do.call(rbind, ls.out)
    return(out)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## * Load results
dt.sim1 <- loadRes(path.results1, tempo.file = TRUE)
dt.sim2 <- loadRes(path.results2, tempo.file = TRUE)
dt.sim3 <- loadRes(path.results3, tempo.file = TRUE)
## dt.sim1[n == 1000 & type=="average" & beta == 0, .N, by = c("type","n","beta","file")]

dt.sim <- rbind(cbind(scenario = "scenario 1", dt.sim1),
                cbind(scenario = "scenario 2", dt.sim2),
                cbind(scenario = "scenario 3", dt.sim3))

## ** check
## unique(dt.sim$n)
dt.sim1[n==1000 & beta == 0 & type %in% c("average","gls"), .(sd(estimate),mean(se)), by = "type"]
dt.sim2[n==1000 & beta == 0 & type %in% c("average","gls"), .(sd(estimate),mean(se)), by = "type"]
dt.sim3[n==1000 & beta == 0 & type %in% c("average","gls"), .(sd(estimate),mean(se)), by = "type"]

## * Process results
dt.sim[, n.char := factor(n/2, unique(n/2))]
dt.sim[, beta.char := ifelse(beta==0,paste0("Null hypothesis (\u03B2=",beta,")"),paste0("Alternative hypothesis (\u03B2=",beta,")"))]
dt.sim[, beta.char := factor(beta.char,levels=unique(beta.char))]
dt.sim[, type.char := factor(type,levels=c("average","fixse","gls","robust","proportion"),c("pool (average)","pool (se)","pool (gls)","pool (robust gls)","proportion"))]
dt.sim[, target := factor(type,levels=c("average","fixse","gls","robust","proportion"),c("pool","pool","pool","pool","proportion"))]
dt.sim[, GS := beta]
dt.sim[type=="proportion", GS := as.numeric(beta>0)]
dtS.sim <- dt.sim[,.(rep = .N, average = mean(estimate), sd = sd(estimate), average.se = mean(se), bias = mean(estimate)-GS, power = mean(p.value<=0.05)),
                    by = c("scenario","beta","n","type","n.char","beta.char","type.char","target","GS")]



## * display results
## bias
gg.bias <- ggplot(dtS.sim, aes(x = n.char, y = average, linetype = target, color = type.char, group = type.char)) 
gg.bias <- gg.bias + geom_hline(aes(yintercept = beta), color = "brown")
gg.bias <- gg.bias + geom_point(size = 2) + geom_line(size = 1) + facet_grid(beta.char~scenario) 
gg.bias <- gg.bias + labs(x = "Sample size (per group)", y = "Estimate", color = "") + guides(linetype = "none")
gg.bias <- gg.bias + scale_color_manual(values = gg_color_hue(5))
gg.bias <- gg.bias + theme(legend.position="bottom",
                           text = element_text(size=15), 
                           axis.line = element_line(size = 1.25),
                           axis.ticks = element_line(size = 2),
                           axis.ticks.length=unit(.25, "cm"),
                           legend.key.size = unit(3,"line"))
gg.bias

## se
gg.se <- ggplot(dtS.sim, aes(x = n.char, y = sd, linetype = target, color = type.char, group = type.char)) 
gg.se <- gg.se + geom_point(size = 2) + geom_line(size = 1) + facet_grid(beta.char~scenario) 
gg.se <- gg.se + labs(x = "Sample size (per group)", y = "Standard deviation", color = "") + guides(linetype = "none")
gg.se <- gg.se + scale_color_manual(values = gg_color_hue(5))
gg.se <- gg.se + theme(legend.position="bottom",
                           text = element_text(size=15), 
                           axis.line = element_line(size = 1.25),
                           axis.ticks = element_line(size = 2),
                           axis.ticks.length=unit(.25, "cm"),
                           legend.key.size = unit(3,"line"))
gg.se

## cali se
gg.cali <- ggplot(dtS.sim[!is.na(dtS.sim$average.se)], aes(x = sd, y = average.se, linetype = target, color = type.char, group = type.char)) 
gg.cali <- gg.cali + geom_abline(intercept = 0, slope = 1, color = "brown")
gg.cali <- gg.cali + geom_point(size = 2) + geom_line(size = 1) + facet_grid(beta.char~scenario) 
gg.cali <- gg.cali + labs(x = "Empirical standard deviation", y = "Average modeled standard deviation", color = "") + guides(linetype = "none")
gg.cali <- gg.cali + scale_color_manual(values = gg_color_hue(5)[1:3])
gg.cali <- gg.cali + theme(legend.position="bottom",
                           text = element_text(size=15), 
                           axis.line = element_line(size = 1.25),
                           axis.ticks = element_line(size = 2),
                           axis.ticks.length=unit(.25, "cm"),
                           legend.key.size = unit(3,"line"))
gg.cali

## relative efficiency
table.reff <- dcast(dtS.sim[type!="proportion",.(scenario,beta.char,n.char,type.char,sd)], scenario+beta.char+n.char~type.char, value.var = "sd")
table.reff[,c("pool (se)", "pool (gls)", "pool (robust gls)") := .(.SD[["pool (se)"]]/.SD[["pool (average)"]],.SD[["pool (gls)"]]/.SD[["pool (average)"]], .SD[["pool (robust gls)"]]/.SD[["pool (average)"]])]
table.reff[, scenario := paste(scenario,gsub("\u03B2","$\\beta$",beta.char, fixed = TRUE),sep=": ")]
table.reff[, beta.char := NULL]
table.reff[duplicated(scenario), scenario := ""]
setnames(table.reff,
         old = c("n.char","pool (average)","pool (se)", "pool (gls)", "pool (robust gls)"),
         new = c("n","sd(average)", "relative sd(se)","relative sd(gls)","relative sd(robust gls)"))

## table.reffscenario
## table.reff$scenario
ggtable.reff <- copy(table.reff)
ggtable.reff[grepl("Alternative",scenario), scenario := gsub("scenario 1: ","\\hphantom{scenario 1:}",scenario, fixed = TRUE)]
ggtable.reff[grepl("Alternative",scenario), scenario := gsub("scenario 2: ","\\hphantom{scenario 2:}",scenario, fixed = TRUE)]
ggtable.reff[grepl("Alternative",scenario), scenario := gsub("scenario 3: ","\\hphantom{scenario 3:}",scenario, fixed = TRUE)]
ggtable.reff[, scenario := gsub("Null hypothesis|Alternative hypothesis","",scenario)]
ggtable.reff[, scenario := gsub(" \\(|)$","",scenario)]
ggtable.reff[["sd(average)"]] <- round(ggtable.reff[["sd(average)"]],3)
ggtable.reff[["relative sd(se)"]] <- paste0(formatC(100*ggtable.reff[["relative sd(se)"]], digits = 1, format = "f"),"\\%")
ggtable.reff[["relative sd(gls)"]] <- paste0(formatC(100*ggtable.reff[["relative sd(gls)"]], digits = 1, format = "f"),"\\%")
ggtable.reff[["relative sd(robust gls)"]] <- paste(formatC(100*ggtable.reff[["relative sd(robust gls)"]], digits = 1, format = "f"),"\\%")


## power
gg.power <- ggplot(dtS.sim[!is.na(dtS.sim$power)], aes(x = n.char, y = power, linetype = target, color = type.char, group = type.char)) 
gg.power <- gg.power + geom_hline(yintercept = 0.05, color = "brown")
gg.power <- gg.power + geom_point(size = 2) + geom_line(size = 1) + facet_grid(beta.char~scenario) 
gg.power <- gg.power + labs(x = "Sample size (per group)", y = "Rejection rate", color = "") + guides(linetype = "none")
gg.power <- gg.power + scale_color_manual(values = gg_color_hue(5)[1:3])
gg.power <- gg.power + theme(legend.position="bottom",
                           text = element_text(size=15), 
                           axis.line = element_line(size = 1.25),
                           axis.ticks = element_line(size = 2),
                           axis.ticks.length=unit(.25, "cm"),
                           legend.key.size = unit(3,"line"))
gg.power


## proportion
dcast(dtS.sim[beta==0 & type == "proportion", .(scenario,n,average)], scenario ~ n, value.var = "average")
##      scenario         20         50        100        200        400       1000
## 1: scenario 1 0.05816622 0.07008068 0.05517511 0.06282701 0.06543944 0.06069456
## 2: scenario 2 0.06018205 0.05423457 0.06226046 0.06809183 0.06392763 0.07618177
## 3: scenario 3 0.04989776 0.06438369 0.05026791 0.05515174 0.06093257 0.05215220

dcast(dtS.sim[beta==0.5 & type == "proportion", .(scenario,n,average)], scenario ~ n, value.var = "average")
##      scenario         20        50       100       200       400      1000
## 1: scenario 1 0.10688060 0.1688707 0.2772736 0.4704577 0.7267768 0.9834127
## 2: scenario 2 0.08157253 0.1074918 0.1819484 0.2366261 0.4113410 0.6355033
## 3: scenario 3 0.07629947 0.1087106 0.1663567 0.2671230 0.4610557 0.7940794


## type 1 error
dcast(dtS.sim[beta==0 & type == "average", .(scenario,n,power)], scenario ~ n, value.var = "power")
##      scenario         20         50        100        200        400       1000
## 1: scenario 1 0.06111987 0.06151420 0.05993691 0.05165615 0.05638801 0.04455836
## 2: scenario 2 0.07071823 0.05745856 0.05856354 0.05082873 0.04088398 0.04309392
## 3: scenario 3 0.06130108 0.06005004 0.05504587 0.05212677 0.05254379 0.05004170

dcast(dtS.sim[beta==0 & type == "fixse", .(scenario,n,power)], scenario ~ n, value.var = "power")
##      scenario         20         50        100        200        400       1000
## 1: scenario 1 0.06664038 0.06506309 0.06033123 0.05362776 0.05678233 0.04455836
## 2: scenario 2 0.09060773 0.04861878 0.06298343 0.06408840 0.05635359 0.05193370
## 3: scenario 3 0.06588824 0.06005004 0.05546289 0.05504587 0.05462886 0.05129274

## * Export

## bias
png(file.path(path,"figures","simulation-fig-bias.png"))
gg.bias
dev.off()

pdf(file.path(path,"figures","simulation-fig-bias.pdf"), width = 10)
gg.bias
dev.off()

## standard error
png(file.path(path,"figures","simulation-fig-dispersion.png"))
gg.se
dev.off()

pdf(file.path(path,"figures","simulation-fig-dispersion.pdf"), width = 10)
gg.se
dev.off()

## calibration se
png(file.path(path,"figures","simulation-fig-calibrationSE.png"))
gg.cali
dev.off()

pdf(file.path(path,"figures","simulation-fig-calibrationSE.pdf"), width = 10)
gg.cali
dev.off()

## # relative efficiency
## png(file.path(path,"figures","simulation-table-reff.png"), width = 550, height = 700)
## grid.arrange(tableGrob(as.data.frame(ggtable.reff), rows = NULL))
## dev.off()

## pdf(file.path(path,"figures","simulation-table-reff.pdf"), width = 550, height = 700)
## grid.arrange(tableGrob(as.data.frame(ggtable.reff), rows = NULL))
## dev.off()

## power
png(file.path(path,"figures","simulation-fig-power.png"))
gg.power
dev.off()

pdf(file.path(path,"figures","simulation-fig-power.pdf"), width = 10)
gg.power
dev.off()


xtable.reff <- capture.output(print(xtable(ggtable.reff, digits = 3),
                                    include.rownames = FALSE,
                                    sanitize.text.function = identity,
                                    booktabs = TRUE))

index.space <- c(which(grepl("2: ",xtable.reff))-1,which(grepl("3: ",xtable.reff))-1)
index.space2 <- which(grepl("phantom",xtable.reff))-1
xtable.reff[index.space] <- paste(xtable.reff[index.space]," [4mm]")
xtable.reff[index.space2] <- paste(xtable.reff[index.space2]," [2mm]")
cat(sapply(xtable.reff,paste0,"\n"))
## % latex table generated in R 4.2.0 by xtable 1.8-4 package
## % Mon Oct 10 10:12:15 2022
## \begin{table}[ht]
## \centering
## \begin{tabular}{llrrrr}
##   \hline
## scenario & n & sd(average) & relative sd(se) & relative sd(gls) & relative sd(robust gls) \\ 
##   \hline
## scenario 1: $\beta$=0 & 10 & 0.555 & 99.7\% & 7.063 & 4.622 \\ 
##    & 25 & 0.349 & 99.9\% & 2.381 & 2.204 \\ 
##    & 50 & 0.248 & 99.9\% & 1.692 & 1.965 \\ 
##    & 100 & 0.177 & 100.0\% & 1.379 & 1.811 \\ 
##    & 200 & 0.125 & 100.0\% & 1.246 & 1.603 \\ 
##    & 500 & 0.079 & 100.0\% & 1.150 & 1.401 \\ 
##   \hphantom{scenario 1:} $\beta$=0.5 & 10 & 0.550 & 99.6\% & 6.882 & 4.685 \\ 
##    & 25 & 0.347 & 100.0\% & 2.414 & 2.207 \\ 
##    & 50 & 0.249 & 100.0\% & 1.698 & 1.989 \\ 
##    & 100 & 0.175 & 100.0\% & 1.369 & 1.789 \\ 
##    & 200 & 0.123 & 100.0\% & 1.243 & 1.653 \\ 
##    & 500 & 0.079 & 100.0\% & 1.153 & 1.385 \\ 
##   scenario 2: $\beta$=0 & 10 & 0.649 & 80.6\% & 0.902 & 1.338 \\ 
##    & 25 & 0.411 & 79.2\% & 0.795 & 1.228 \\ 
##    & 50 & 0.294 & 78.7\% & 0.765 & 1.199 \\ 
##    & 100 & 0.206 & 79.3\% & 0.766 & 1.176 \\ 
##    & 200 & 0.145 & 79.6\% & 0.766 & 1.170 \\ 
##    & 500 & 0.092 & 79.1\% & 0.755 & 1.168 \\ 
##   \hphantom{scenario 2:} $\beta$=0.5 & 10 & 0.657 & 80.4\% & 0.898 & 1.348 \\ 
##    & 25 & 0.411 & 79.7\% & 0.799 & 1.225 \\ 
##    & 50 & 0.291 & 79.7\% & 0.780 & 1.190 \\ 
##    & 100 & 0.206 & 78.4\% & 0.751 & 1.179 \\ 
##    & 200 & 0.144 & 79.3\% & 0.762 & 1.175 \\ 
##    & 500 & 0.092 & 78.9\% & 0.755 & 1.169 \\ 
##   scenario 3: $\beta$=0 & 10 & 0.697 & 99.9\% & 8.212 & 6.222 \\ 
##    & 25 & 0.442 & 99.9\% & 2.770 & 2.684 \\ 
##    & 50 & 0.314 & 100.0\% & 1.885 & 2.089 \\ 
##    & 100 & 0.224 & 100.5\% & 1.534 & 1.771 \\ 
##    & 200 & 0.156 & 100.7\% & 1.377 & 1.556 \\ 
##    & 500 & 0.099 & 101.0\% & 1.235 & 1.340 \\ 
##   \hphantom{scenario 3:} $\beta$=0.5 & 10 & 0.686 & 100.0\% & 8.422 & 6.479 \\ 
##    & 25 & 0.446 & 101.0\% & 2.647 & 2.523 \\ 
##    & 50 & 0.313 & 100.4\% & 1.927 & 2.164 \\ 
##    & 100 & 0.223 & 100.5\% & 1.540 & 1.813 \\ 
##    & 200 & 0.153 & 100.3\% & 1.365 & 1.574 \\ 
##    & 500 & 0.100 & 100.4\% & 1.240 & 1.334 \\ 
##    \hline
## \end{tabular}
## \end{table}
##----------------------------------------------------------------------
### BUILD_simulation.R ends here
