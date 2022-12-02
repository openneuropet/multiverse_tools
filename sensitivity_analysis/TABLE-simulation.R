### TABLE-simulation.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 12 2022 (11:52) 
## Version: 
## Last-Updated: dec  1 2022 (09:38) 
##           By: Brice Ozenne
##     Update #: 11
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Dependency
library(data.table)
library(xtable)

## * Load data
dtS.sim <- readRDS(file.path("Results","simulation-summary-scenario.rds"))


## * Table efficiency
table.reff <- dcast(dtS.sim[type!="proportion",.(scenario,beta.char,n.char,type.char,sd)], scenario+beta.char+n.char~type.char, value.var = "sd")
table.reff[,c("pool (se)", "pool (gls)", "pool (constrained gls)", "pool (robust gls)") := .(.SD[["pool (se)"]]/.SD[["pool (average)"]],
                                                                                             .SD[["pool (gls)"]]/.SD[["pool (average)"]],
                                                                                             .SD[["pool (constrained gls)"]]/.SD[["pool (average)"]],
                                                                                             .SD[["pool (robust gls)"]]/.SD[["pool (average)"]])]
table.reff[, scenario := paste(scenario,gsub("\u03B2","$\\beta$",beta.char, fixed = TRUE),sep=": ")]
table.reff[, beta.char := NULL]
table.reff[, J := 6]
table.reff[grepl("scenario 1|scenario 3",scenario), J := 20]
table.reff[duplicated(scenario), scenario := ""]
table.reff$n.char <- ifelse(nchar(as.character(table.reff$n.char))==3,as.character(table.reff$n.char),paste0(as.character(table.reff$n.char)," "))
table.reff$n.char <- paste0(table.reff$n.char, " (",round(as.numeric(as.character(table.reff$n))/table.reff$J,2),")")
table.reff$n.char <- factor(table.reff$n.char,unique(table.reff$n))
table.reff[, J := NULL]
setnames(table.reff,
         old = c("n.char","pool (average)","pool (se)", "pool (gls)", "pool (constrained gls)", "pool (robust gls)"),
         new = c("n  (n/J)","sd(average)", "relative sd(se)","relative sd(gls)","relative sd(constrained gls)","relative sd(robust gls)"))

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
ggtable.reff[["relative sd(constrained gls)"]] <- paste(formatC(100*ggtable.reff[["relative sd(constrained gls)"]], digits = 1, format = "f"),"\\%")
ggtable.reff[["relative sd(robust gls)"]] <- paste(formatC(100*ggtable.reff[["relative sd(robust gls)"]], digits = 1, format = "f"),"\\%")

## * Table proportion
table.proportionH0 <- dcast(dtS.sim[beta==0 & type == "proportion", .(scenario,n,average)], scenario ~ n, value.var = "average")
table.proportionH0
##      scenario         20         50        100        200        400       1000
## 1: scenario 1 0.06311754 0.06766368 0.06223540 0.06991387 0.06437179 0.06184707
## 2: scenario 2 0.06329151 0.06174227 0.06234090 0.06331214 0.06430527 0.06332216
## 3: scenario 3 0.05605482 0.06004671 0.05534097 0.06227068 0.05925187 0.05618175
range(as.data.frame(table.proportionH0)[-1])
## [1] 0.05534097 0.06991387

table.proportionH1 <- dcast(dtS.sim[beta==0.5 & type == "proportion", .(scenario,n,average)], scenario ~ n, value.var = "average")
table.proportionH1
##      scenario         20        50       100       200       400      1000
## 1: scenario 1 0.10209858 0.1663022 0.2765966 0.4652930 0.7379401 0.9819058
## 2: scenario 2 0.07982463 0.1165405 0.1713475 0.2580254 0.4017544 0.6342199
## 3: scenario 3 0.07455742 0.1082825 0.1696291 0.2760800 0.4681342 0.7911460
range(as.data.frame(table.proportionH1)[-1])
## [1] 0.07455742 0.98190578

## * Table type 1 error
table.type1Ag <- dcast(dtS.sim[beta==0 & type == "average", .(scenario,n,power)], scenario ~ n, value.var = "power")
table.type1Ag
##      scenario     20     50    100    200    400   1000
## 1: scenario 1 0.0631 0.0599 0.0531 0.0521 0.0507 0.0519
## 2: scenario 2 0.0663 0.0554 0.0542 0.0499 0.0456 0.0472
## 3: scenario 3 0.0641 0.0579 0.0526 0.0516 0.0508 0.0515
range(as.data.frame(table.type1Av)[-1])
## [1] 0.0456 0.0663

table.type1Se <- dcast(dtS.sim[beta==0 & type == "fixse", .(scenario,n,power)], scenario ~ n, value.var = "power")
table.type1Se
##      scenario     20     50    100    200    400   1000
## 1: scenario 1 0.0672 0.0612 0.0545 0.0528 0.0506 0.0522
## 2: scenario 2 0.0802 0.0569 0.0547 0.0533 0.0491 0.0501
## 3: scenario 3 0.0687 0.0618 0.0525 0.0519 0.0503 0.0522
range(as.data.frame(table.type1Se)[-1])
## [1] 0.0491 0.0802

table.type1GLS <- dcast(dtS.sim[beta==0 & type == "gls", .(scenario,n,power)], scenario ~ n, value.var = "power")
table.type1GLS
##      scenario     20     50    100    200    400   1000
## 1: scenario 1 0.7716 0.2481 0.1169 0.0779 0.0657 0.0576
## 2: scenario 2 0.1785 0.0815 0.0674 0.0613 0.0528 0.0515
## 3: scenario 3 0.7709 0.2497 0.1146 0.0788 0.0624 0.0570
range(as.data.frame(table.type1GLS)[-1])
## [1] 0.0515 0.7716

## * export
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
##  % Thu Dec  1 09:33:54 2022
##  \begin{table}[ht]
##  \centering
##  \begin{tabular}{llrllll}
##    \toprule
##  scenario & n  (n/J) & sd(average) & relative sd(se) & relative sd(gls) & relative sd(constrained gls) & relative sd(robust gls) \\ 
##    \midrule
##  scenario 1: $\beta$=0 & 10  (0.5) & 0.554 & 99.8\% & 259.1\% & 104.1 \% & 184.3 \% \\ 
##     & 25  (1.25) & 0.355 & 99.9\% & 113.0\% & 98.5 \% & 113.7 \% \\ 
##     & 50  (2.5) & 0.250 & 99.9\% & 96.4\% & 94.6 \% & 105.6 \% \\ 
##     & 100 (5) & 0.177 & 100.0\% & 90.9\% & 90.8 \% & 99.4 \% \\ 
##     & 200 (10) & 0.125 & 100.0\% & 89.2\% & 89.2 \% & 95.5 \% \\ 
##     & 500 (25) & 0.079 & 100.0\% & 87.9\% & 87.9 \% & 91.3 \% \\   [2mm]
##    \hphantom{scenario 1:}$\beta$=0.5 & 10  (0.5) & 0.557 & 99.8\% & 255.8\% & 104.3 \% & 182.6 \% \\ 
##     & 25  (1.25) & 0.353 & 99.9\% & 113.1\% & 98.8 \% & 113.1 \% \\ 
##     & 50  (2.5) & 0.250 & 99.9\% & 96.6\% & 94.8 \% & 105.5 \% \\ 
##     & 100 (5) & 0.177 & 100.0\% & 91.6\% & 91.5 \% & 99.8 \% \\ 
##     & 200 (10) & 0.125 & 100.0\% & 89.8\% & 89.8 \% & 95.9 \% \\ 
##     & 500 (25) & 0.079 & 100.0\% & 88.0\% & 88.0 \% & 91.3 \% \\   [4mm]
##    scenario 2: $\beta$=0 & 10  (1.67) & 0.659 & 80.2\% & 89.0\% & 87.3 \% & 134.2 \% \\ 
##     & 25  (4.17) & 0.411 & 78.9\% & 78.8\% & 78.6 \% & 122.7 \% \\ 
##     & 50  (8.33) & 0.295 & 79.0\% & 77.1\% & 77.1 \% & 119.5 \% \\ 
##     & 100 (16.67) & 0.206 & 79.6\% & 76.9\% & 76.9 \% & 117.7 \% \\ 
##     & 200 (33.33) & 0.144 & 79.8\% & 76.8\% & 76.8 \% & 117.2 \% \\ 
##     & 500 (83.33) & 0.091 & 79.4\% & 76.1\% & 76.1 \% & 116.7 \% \\   [2mm]
##    \hphantom{scenario 2:}$\beta$=0.5 & 10  (1.67) & 0.652 & 80.4\% & 90.1\% & 88.4 \% & 135.8 \% \\ 
##     & 25  (4.17) & 0.407 & 79.2\% & 79.5\% & 79.3 \% & 122.8 \% \\ 
##     & 50  (8.33) & 0.290 & 79.4\% & 77.8\% & 77.8 \% & 118.6 \% \\ 
##     & 100 (16.67) & 0.207 & 79.1\% & 76.1\% & 76.1 \% & 117.7 \% \\ 
##     & 200 (33.33) & 0.143 & 79.0\% & 75.9\% & 75.9 \% & 117.3 \% \\ 
##     & 500 (83.33) & 0.093 & 78.9\% & 75.3\% & 75.3 \% & 116.8 \% \\   [4mm]
##    scenario 3: $\beta$=0 & 10  (0.5) & 0.698 & 100.0\% & 208.1\% & 100.5 \% & 197.9 \% \\ 
##     & 25  (1.25) & 0.444 & 100.2\% & 91.6\% & 86.7 \% & 113.2 \% \\ 
##     & 50  (2.5) & 0.313 & 100.1\% & 78.1\% & 78.0 \% & 108.6 \% \\ 
##     & 100 (5) & 0.221 & 100.3\% & 74.0\% & 74.0 \% & 106.9 \% \\ 
##     & 200 (10) & 0.157 & 100.2\% & 71.7\% & 71.7 \% & 107.6 \% \\ 
##     & 500 (25) & 0.099 & 100.2\% & 71.2\% & 71.2 \% & 108.8 \% \\   [2mm]
##    \hphantom{scenario 3:}$\beta$=0.5 & 10  (0.5) & 0.694 & 100.2\% & 212.2\% & 101.3 \% & 200.5 \% \\ 
##     & 25  (1.25) & 0.439 & 100.4\% & 92.1\% & 87.1 \% & 114.2 \% \\ 
##     & 50  (2.5) & 0.312 & 100.3\% & 78.7\% & 78.5 \% & 107.9 \% \\ 
##     & 100 (5) & 0.222 & 100.2\% & 73.8\% & 73.8 \% & 107.2 \% \\ 
##     & 200 (10) & 0.156 & 100.3\% & 72.7\% & 72.7 \% & 108.4 \% \\ 
##     & 500 (25) & 0.099 & 100.4\% & 71.6\% & 71.6 \% & 108.2 \% \\ 
##     \bottomrule
##  \end{tabular}
##  \end{table}
##----------------------------------------------------------------------
### TABLE-simulation.R ends here
