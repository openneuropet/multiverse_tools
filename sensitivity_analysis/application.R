### analysis-pipeline.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug 24 2022 (17:53) 
## Version: 
## Last-Updated: nov  7 2022 (17:31) 
##           By: Brice Ozenne
##     Update #: 39
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

## From martinnoergaard (slack 12-09-2022)
## I have extracted the following 28 regions, that were averaged between left and right.
## This corresponds to the column indexes [6, 13, 1, 8, 3, 10, 2, 9, 39, 73, 5, 12, 27, 61, 24, 58, 41, 75, 43, 77, 48, 82, 22, 56, 42, 76, 19, 53]

## ** pipeline data
dico.region <- matrix(c(6, 13, 1, 8, 3, 10, 2, 9, 39, 73, 5, 12, 27, 61, 24, 58, 41, 75, 43, 77, 48, 82, 22, 56, 42, 76, 19, 53), ncol = 2, byrow = TRUE,
                      dimnames = list(c("amygdala", "thalamus", "putamen", "caudate", "ACC", "hippocampus", "OFC", "SFC", "OC", "STG", "insula", "ITG", "PC", "EC"), c("left","right")))

ls.dfraw <- c(lapply(1:4, function(iPip){ ## iPip <- 1
    cbind(pipeline = iPip, read_xls(file.path("data","correlatedData.xls"), sheet = iPip, .name_repair = "minimal"))
}),
lapply(5:8, function(iPip){
    cbind(pipeline = iPip, read_xls(file.path("data","uncorrelatedData.xls"), sheet = iPip-4, .name_repair = "minimal"))
})
)

dtraw <- as.data.table(do.call(rbind,ls.dfraw))
dtav <- data.table(pipeline = dtraw[[1]],
                  index = dtraw[[2]],
                  (dtraw[,.SD,.SDcols = colnames(dtraw)[2+dico.region[,1]]]+dtraw[,.SD,.SDcols = colnames(dtraw)[2+dico.region[,2]]])/2)
names(dtav)[3:NCOL(dtav)] <- rownames(dico.region)

## ** atlas data
dfAll.atlas <- read.csv(file.path("data","bp_table.csv"))
df.atlas <- dfAll.atlas[-1,1:10]
names(df.atlas) <- c("mean.DASB","sd.DASB","mean.cumi","sd.cumi","mean.az","sd.az","mean.cimbi","sd.cimbi","mean.sb","sd.sb")

type.region <- which(rowSums(df.atlas!="")==0)
reptype.region <- diff(c(type.region,NROW(df.atlas)+1))

dt.atlas <- cbind(region = rownames(df.atlas),
                  type = unlist(lapply(1:length(type.region), function(iR){rep(names(type.region[iR]),reptype.region[iR])})),
                  data.table(do.call(cbind,lapply(df.atlas,as.numeric))))[region!=type]


## * Process data
dtW <- dtav[pipeline!=5]

## * Analysis

## ** region-specific
name.region <- rownames(dico.region)

name.region[name.region %in% tolower(dt.atlas$region) == FALSE]


ls.mlmmRS <- lapply(name.region, function(iRegion){ ## iRegion <- "amygdala"
    
    if(iRegion %in% tolower(dt.atlas$region)){
        null.value <- dt.atlas[tolower(region) == iRegion, mean.DASB]
    }else{
        null.value <- 0
    }

    iFormula <- as.formula(paste0(iRegion,"~1"))
    e.mlmm <- do.call("mlmm", list(iFormula, repetition = ~1|index, data = dtW, by = "pipeline", effects = paste0("(Intercept)=",null.value)))

})
## any pipeline with an effect
summary(ls.mlmmRS[[1]], method = "single-step", columns = c("estimate","se","df","null","lower","upper","p.value"))
##                         estimate    se df lower upper null p.value
## pipeline=1: (Intercept)    2.118 0.056 59 1.995  2.24    2  0.0614
## pipeline=2: (Intercept)    2.102 0.056 59  1.98 2.225    2  0.1127
## pipeline=3: (Intercept)    2.106 0.054 59 1.986 2.227    2  0.0883
## pipeline=4: (Intercept)    2.093 0.054 59 1.974 2.212    2  0.1447
## pipeline=6: (Intercept)    3.157 0.048 59 3.052 3.263    2  <2e-16
## pipeline=7: (Intercept)    2.809 0.039 59 2.724 2.895    2  <2e-16
## pipeline=8: (Intercept)    2.549 0.032 59 2.478 2.621    2  <2e-16

## at least one pipeline with no effect
res.mlmmRS <- model.tables(ls.mlmmRS[[1]], method = "none")
max(res.mlmmRS$p.value)
## [1] 0.09068143

## proportion of pipelines with an effect
proportion(ls.mlmmRS[[1]], method = "single-step", n.sample = 0)
##             estimate        se     lower     upper
## proportion 0.6498292 0.1379316 0.4616139 0.9221975

## ** global
dtL <- melt(dtW, id.vars = c("pipeline","index"))

e.mlmmG <- mlmm(value ~ 0+variable, repetition = ~variable|index, data = dtL, by = "pipeline", df = FALSE)
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
ls.forestRS <- lapply(ls.mlmmRS, function(eMLMM){ ## eMLMM <- ls.mlmmRS[[1]]
    
    iGG <- plot(eMLMM, type = "forest", plot = FALSE)$plot + ggtitle(eMLMM$object$outcome)
    if(eMLMM$object$outcome %in% c("amygdala","thalamus","putamen","caudate")){
      iGG <- iGG + theme(plot.margin=unit(c(0,0,-0.25,0),"cm"))
    }else{
      iGG <- iGG + theme(plot.margin=unit(c(-0.25,0,-0.25,0),"cm"))
    }
    iGG + scale_x_discrete(breaks = rownames(eMLMM$univariate),
                           labels = gsub("pipeline","pip",gsub("(Intercept)", "", rownames(eMLMM$univariate), fixed = TRUE)))
})

gg.forestRS <- do.call(ggarrange, ls.forestRS)

ls.heatRS <- lapply(ls.mlmmRS, function(eMLMM){ ## eMLMM <- ls.mlmmRS[[1]]
    iPlot <- plot(eMLMM, type = "heat", plot = FALSE)
    iGG <- iPlot$plot + ggtitle(eMLMM$object$outcome) + xlab("") + ylab("")
    if(eMLMM$object$outcome %in% c("amygdala","thalamus","putamen","caudate")){
      iGG <- iGG + theme(plot.margin=unit(c(0,0,-0.25,0),"cm"))
    }else{
      iGG <- iGG + theme(plot.margin=unit(c(-0.25,0,-0.25,0),"cm"))
    }
    iGG <- iGG + scale_y_discrete(breaks = unique(iPlot$data$col),
                                  labels = gsub("pipeline", "pip=", unique(iPlot$data$col), fixed = TRUE))
    iGG <- iGG + scale_x_discrete(breaks = unique(iPlot$data$row),
                                  labels = gsub("pipeline=", "", unique(iPlot$data$row), fixed = TRUE))
})

gg.heatRS <- do.call(ggarrange, c(ls.heatRS, list(common.legend = TRUE)))

## export
pdf(file.path("figures","application-figure-heatmapRS.pdf"), width = 12)
gg.heatRS
dev.off()

png(file.path("figures","application-figure-heatmapRS.png"), width = 1000, heigh = 650)
gg.heatRS
dev.off()

pdf(file.path("figures","application-figure-forestRS.pdf"), width = 12)
gg.forestRS
dev.off()

png(file.path("figures","application-figure-forestRS.png"), width = 1000, heigh = 650)
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


pdf(file.path("figures","application-figure-heatmapG.pdf"), width = 12)
gg.heatG
dev.off()
png(file.path("figures","application-figure-heatmapG.png"), heigh = 650, width = 1000)
gg.heatG
dev.off()


##----------------------------------------------------------------------
### analysis-pipeline.R ends here
