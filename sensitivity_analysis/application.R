### analysis-pipeline.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug 24 2022 (17:53) 
## Version: 
## Last-Updated: sep 14 2022 (18:14) 
##           By: Brice Ozenne
##     Update #: 29
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(data.table)
library(ggplot2)
library(readxl)
library(LMMstar)
library(ggpubr)

## * Import data
## 14 regions: amygdala, thalamus, putamen, caudate, anterior cingulate cortex (ACC),
##             hippocampus, orbital frontal cortex, superior frontal cortex, occipital cortex, superior temporal gyrus,
##             insula, inferior temporal gyrus, parietal cortex, entorhinal cortex

ls.dfraw <- c(lapply(1:4, function(iPip){ ## iPip <- 1
    cbind(pipeline = iPip, read_xls(file.path("data","correlatedData.xls"), sheet = iPip, .name_repair = "minimal"))
}),
lapply(5:8, function(iPip){
    cbind(pipeline = iPip, read_xls(file.path("data","uncorrelatedData.xls"), sheet = iPip-4, .name_repair = "minimal"))
})
)

dtraw <- as.data.table(do.call(rbind,ls.dfraw))
names(dtraw)[2] <- "index"
names(dtraw) <- gsub("-","_",names(dtraw))

## * Process data
dtW <- dtraw[pipeline!=5, .(pipeline = pipeline,
                            index = index,
                            amygdala = (Left_Amygdala+Right_Amygdala)/2,
                            thalamus = (Left_Thalamus_Proper+Right_Thalamus_Proper)/2,
                            putamen = (Left_Putamen+Right_Putamen)/2,
                            caudate = (Left_Caudate+Right_Caudate)/2,
                            hippocampus = (Left_Hippocampus+Right_Hippocampus)/2)]

## * Analysis

## ** region-specific
name.region <- c("amygdala","thalamus","putamen","caudate")

ls.mlmmRS <- lapply(name.region, function(iRegion){ ## iRegion <- "amygdala"
    iFormula <- as.formula(paste0(iRegion,"~1"))
    e.mlmm <- mlmm(iFormula, repetition = ~1|index, data = dtW, by = "pipeline", effects = "(Intercept)=2")
})
## any pipeline with an effect
summary(ls.mlmmRS[[1]], method = "single-step", columns = c("estimate","se","df","null","lower","upper","p.value"))

## at least one pipeline with no effect
res.mlmmRS <- model.tables(ls.mlmmRS[[1]], method = "none")
max(res.mlmmRS$p.value)

## proportion of pipelines with an effect
resAdj.mlmmRS <- model.tables(ls.mlmmRS[[1]], method = "single-step")
Wald <- (resAdj.mlmmRS$estimate-2)/resAdj.mlmmRS$se
tc <- attr(resAdj.mlmmRS,"quantile")

mean(resAdj.mlmmRS$p.value <= 0.05)
mean(abs(Wald) >= tc)
## [1] 0.4285714

1-mean(sapply(Wald, function(iStat){pt(tc - iStat, df = resAdj.mlmmRS$df)-pt(-tc - iStat, df = resAdj.mlmmRS$df)}))
## [1] 0.6481363

## ** global
dtL <- melt(dtW, id.vars = c("pipeline","index"))

e.mlmmG <- mlmm(value ~ 0+variable, repetition = ~variable|index, data = dtL, by = "pipeline")
## summary(e.mlmmG)
   ##                                 estimate    se df lower upper p.value    
   ## pipeline=1: variableamygdala       2.118 0.056 59 2.007 2.229  <2e-16 ***
   ## pipeline=1: variablethalamus       2.019 0.032 59 1.955 2.083  <2e-16 ***
   ## pipeline=1: variableputamen        2.284 0.037 59 2.209 2.358  <2e-16 ***
   ## pipeline=1: variablecaudate        1.817 0.034 59 1.748 1.886  <2e-16 ***
   ## pipeline=1: variablehippocampus    0.763 0.018 59 0.727 0.799  <2e-16 ***
## summary(e.mlmmG[[1]])
##    pipeline=1: (Intercept)    2.118 0.056 59 2.007 2.229  <2e-16 ***
## summary(e.mlmmG[[2]])
   ## pipeline=1: (Intercept)    2.019 0.032 59 1.955 2.083  <2e-16 ***


eigen(cov2cor(vcov(e.mlmmG)))$values




## * figures

## ** Region specific
ls.forestRS <- lapply(ls.mlmmRS, function(eMLMM){ ## eMLMM <- ls.mlmm[[1]]
    iGG <- plot(eMLMM, type = "forest", plot = FALSE)$plot + ggtitle(eMLMM$object$outcome)
    iGG + scale_x_discrete(breaks = rownames(eMLMM$univariate),
                           labels = gsub("(Intercept)", "", rownames(eMLMM$univariate), fixed = TRUE))
})

gg.forestRS <- do.call(ggarrange, ls.forestRS)

ls.heatRS <- lapply(ls.mlmmRS, function(eMLMM){ ## eMLMM <- ls.mlmm[[1]]
    iPlot <- plot(eMLMM, type = "heat", plot = FALSE)
    iGG <- iPlot$plot + ggtitle(eMLMM$object$outcome) + xlab("") + ylab("")
    iGG <- iGG + scale_x_discrete(breaks = unique(iPlot$data$row),
                                  labels = gsub("pipeline=", "", unique(iPlot$data$row), fixed = TRUE))
})

gg.heatRS <- do.call(ggarrange, ls.heatRS)

## export
pdf(file.path("figures","heatmapRS.pdf"), width = 12)
gg.heatRS
dev.off()

png(file.path("figures","heatmapRS.png"), width = 1000, heigh = 650)
gg.heatRS
dev.off()

pdf(file.path("figures","forestRS.pdf"), width = 12)
gg.forestRS
dev.off()

png(file.path("figures","forestRS.png"), width = 1000, heigh = 650)
gg.forestRS
dev.off()

## ** Global
pp.heatG <- plot(e.mlmmG, type  = "heat", plot = FALSE)
gg.heatG <- pp.heatG$plot + xlab("") + ylab("")
gg.heatG <- gg.heatG + scale_y_discrete(breaks = unique(pp.heatG$data$row),
                              labels = gsub("variable", "", unique(pp.heatG$data$row), fixed = TRUE))
txt <- gsub(" |\\:","",gsub("pipeline=", "", gsub("variable", "", unique(pp.heatG$data$col), fixed = TRUE), fixed = TRUE))
gg.heatG <- gg.heatG + scale_x_discrete(breaks = unique(pp.heatG$data$col),
                              labels = substr(txt,1,2))
gg.heatG


pdf(file.path("figures","heatmapG.pdf"), width = 12)
gg.heatG
dev.off()
png(file.path("figures","heatmapG.png"), heigh = 650, width = 1000)
gg.heatG
dev.off()


##----------------------------------------------------------------------
### analysis-pipeline.R ends here
