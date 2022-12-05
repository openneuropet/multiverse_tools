### analysis-pipeline.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug 24 2022 (17:53) 
## Version: 
## Last-Updated: Dec  5 2022 (11:18) 
##           By: Brice Ozenne
##     Update #: 70
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
library(egg) ## ggarrange
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


## * Dictionnary for regions/pipeline
## 14 regions: amygdala, thalamus, putamen, caudate, anterior cingulate cortex (ACC),
##             hippocampus, orbital frontal cortex, superior frontal cortex, occipital cortex, superior temporal gyrus,
##             insula, inferior temporal gyrus, parietal cortex, entorhinal cortex

## From martinnoergaard (slack 12-09-2022)
## I have extracted the following 28 regions, that were averaged between left and right.
## This corresponds to the column indexes [6, 13, 1, 8, 3, 10, 2, 9, 39, 73, 5, 12, 27, 61, 24, 58, 41, 75, 43, 77, 48, 82, 22, 56, 42, 76, 19, 53]

dicoMatin.region <- matrix(c(6, 13, 1, 8, 3, 10, 2, 9, 39, 73, 5, 12, 27, 61, 24, 58, 41, 75, 43, 77, 48, 82, 22, 56, 42, 76, 19, 53), ncol = 2, byrow = TRUE,
                           dimnames = list(c("amygdala", "thalamus", "putamen", "caudate", "ACC", "hippocampus", "OFC", "SFC", "OC", "STG", "insula", "ITG", "PC", "EC"), c("left","right")))
## Switch around OC and SFC

rename.region <-  as.data.frame(rbind(c(combined = "amygdala",    left = "Left-Amygdala",            right = "Right-Amygdala", full = "Amygdala"),
                                      c(combined = "thalamus",    left = "Left-Thalamus-Proper",     right = "Right-Thalamus-Proper", full = "Thalamus"),
                                      c(combined = "putamen",     left = "Left-Putamen",             right = "Right-Putamen", full = "Putamen"),           
                                      c(combined = "caudate",     left = "Left-Caudate",             right = "Right-Caudate", full = "Caudate"),           
                                      c(combined = "ACC",         left = "rostralanteriorcingulate", right = "rostralanteriorcingulate2", full = "Anterior cingulate cortex"),
                                      c(combined = "hippocampus", left = "Left-Hippocampus",         right = "Right-Hippocampus", full = "Hippocampus"),
                                      c(combined = "OFC",         left = "medialorbitofrontal",      right = "medialorbitofrontal2", full = "Orbital frontal cortex"),     
                                      c(combined = "SFC",         left = "superiorfrontal",          right = "superiorfrontal2", full = "Superior frontal cortex"),        
                                      c(combined = "OC",          left = "lateraloccipital",         right = "lateraloccipital2", full = "Occipital cortex"),         
                                      c(combined = "STG",         left = "superiortemporal",         right = "superiortemporal2", full = "Superior temporal gyrus"),        
                                      c(combined = "insula",      left = "insula",                   right = "insula2", full = "Insula"),                  
                                      c(combined = "ITG",         left = "inferiortemporal",         right = "inferiortemporal2", full = "Inferior temporal gyrus"),        
                                      c(combined = "PC",          left = "superiorparietal",         right = "superiorparietal2", full = "Parietal cortex"),        
                                      c(combined = "EC",          left = "entorhinal",               right = "entorhinal2", full = "Entorhinal cortex")
                                      ))

                    

## * Import data

## ** pipeline data
ls.dfraw <- c(lapply(1:4, function(iPip){ ## iPip <- 1
    cbind(pipeline.index = iPip, read_xls(file.path("data","correlatedData_intervention.xls"), sheet = iPip, .name_repair = "minimal"))
}),
lapply(5:9, function(iPip){
    cbind(pipeline.index = iPip, read_xls(file.path("data","uncorrelatedData_intervention.xls"), sheet = iPip-4, .name_repair = "minimal"))
})
)

## remove pipeline 5 which is the same as pipeline 1
range(ls.dfraw[[1]][,-(1:2)]-ls.dfraw[[5]][,-(1:2)])
ls.dfraw[[5]] <- NULL

## check large values and put NA there
M.checkLarge <- do.call(rbind,lapply(ls.dfraw, function(iData){apply(abs(iData[,colnames(iData) %in% unlist(rename.region[2:3])]),2,max)}))
rownames(M.checkLarge) <- c(paste0("uncorrelated",1:4),paste0("correlated",1:4))
round(M.checkLarge[,colSums(M.checkLarge>1e3)>0])

## collect data from all pipelines
dtraw <- as.data.table(do.call(rbind,ls.dfraw))

## add 2 into column names to avoid duplicated names
names(dtraw)[duplicated(names(dtraw))] <- paste0(names(dtraw)[duplicated(names(dtraw))],"2")

## add 0 in the index variable to facilitate ordering, i.e. D1 -> D01 which in alphabetic order is before D10
for(iN in 1:9){
    dtraw[[2]] <- gsub(paste0("D",iN,"$"),paste0("D0",iN),dtraw[[2]])
}

## average data between hemispheres
all(rename.region$left %in% names(dtraw))
all(rename.region$right %in% names(dtraw))


dtav <- data.table(index.pipeline = dtraw$pipeline.index,
                   pipeline = factor(dtraw$pipeline.index, levels = c(1:4,6:9), paste0("pipeline ", 1:8)),
                   index = dtraw[[2]],
                   0.5*dtraw[,.SD,.SDcols = rename.region$left] + 0.5*dtraw[,.SD,.SDcols = rename.region$right]
                   )
names(dtav)[4:NCOL(dtav)] <- rename.region$combined
setkeyv(dtav,c("index.pipeline","index"))
dtW <- dtav
## table(dtav$pipeline,dtav$index)

dtL <- melt(dtW, id.vars = c("index.pipeline","index","pipeline"))

## ** atlas data
dfAll.atlas <- read.csv(file.path("data","bp_table.csv"))
df.atlas <- dfAll.atlas[-1,1:10]
names(df.atlas) <- c("mean.DASB","sd.DASB","mean.cumi","sd.cumi","mean.az","sd.az","mean.cimbi","sd.cimbi","mean.sb","sd.sb")

type.region <- which(rowSums(df.atlas!="")==0)
reptype.region <- diff(c(type.region,NROW(df.atlas)+1))

dt.atlas <- cbind(region = tolower(rownames(df.atlas)),
                  type = unlist(lapply(1:length(type.region), function(iR){rep(names(type.region[iR]),reptype.region[iR])})),
                  data.table(do.call(cbind,lapply(df.atlas,as.numeric))))[region!=type]

rename.region
## 5          ACC rostralanteriorcingulate rostralanteriorcingulate2
## 7          OFC      medialorbitofrontal      medialorbitofrontal2
## 8          SFC          superiorfrontal          superiorfrontal2
## 9           OC         lateraloccipital         lateraloccipital2
## 10         STG         superiortemporal         superiortemporal2
## 12         ITG         inferiortemporal         inferiortemporal2
## 13          PC         superiorparietal         superiorparietal2
## 14          EC               entorhinal               entorhinal2


dt.atlas$region[dt.atlas$region == "rostral anterior"] <- "ACC"
dt.atlas$region[dt.atlas$region == "medial orbito frontal"] <- "OFC"
dt.atlas$region[dt.atlas$region == "superior frontal"] <- "SFC"
dt.atlas$region[dt.atlas$region == "lateral occipital"] <- "OC"
dt.atlas$region[dt.atlas$region == "superior temporal"] <- "STG"
dt.atlas$region[dt.atlas$region == "inferior temporal"] <- "ITG"
dt.atlas$region[dt.atlas$region == "superior parietal"] <- "PC" ## inferior or superior ? 
dt.atlas$region[dt.atlas$region == "entorhinal"] <- "EC"

dtR.atlas <- dt.atlas[region %in%  colnames(dtav)]


## * Analysis

## ** region-specific
name.region <- rename.region$combined

ls.mlmmRS <- setNames(lapply(name.region, function(iRegion){ ## iRegion <- "amygdala"

    print(iRegion)
    null.value <- dtR.atlas[region == iRegion, mean.DASB]
    iFormula <- as.formula(paste0(iRegion,"~1"))
    e.mlmm <- do.call("mlmm", list(iFormula, repetition = ~1|index, data = dtW, by = "pipeline", effects = paste0("(Intercept)=",null.value), trace = FALSE))    
    return(e.mlmm)

}), name.region)

## *** any pipeline with an effect (slow!)
ls.SmlmmRS <- lapply(ls.mlmmRS, FUN = function(iLMM){confint(iLMM, method = "single-step", columns = c("estimate","se","df","null","lower","upper","p.value"))})

table.estimate <- round(do.call(cbind,lapply(ls.SmlmmRS,"[",1)),3)
table.pvalue <- round(do.call(cbind,lapply(ls.SmlmmRS,"[","p.value")),3)
table.pvalue[table.pvalue==0] <- "<0.001"

colnames(table.estimate) <- names(ls.SmlmmRS)
colnames(table.pvalue) <- names(ls.SmlmmRS)

table.estimate[1:7];table.estimate[8:12]
##                                        amygdala thalamus putamen caudate   ACC hippocampus   OFC
## pipeline 1 (correlated)  : (Intercept)    2.044    2.021   2.232   1.738 0.862       0.754 0.700
## pipeline 2 (correlated)  : (Intercept)    2.035    2.030   2.231   1.747 0.864       0.751 0.696
## pipeline 3 (correlated)  : (Intercept)    2.025    1.961   2.193   1.690 0.839       0.733 0.680
## pipeline 4 (correlated)  : (Intercept)    2.020    1.978   2.198   1.706 0.843       0.733 0.679
## pipeline 5 (uncorrelated): (Intercept)    1.851    2.284   2.150   1.569 0.855       0.783 0.676
## pipeline 6 (uncorrelated): (Intercept)    1.828    2.274   2.134   1.564 0.841       0.821 0.667
## pipeline 7 (uncorrelated): (Intercept)    2.001    2.008   2.208   1.738 0.865       0.769 0.697
## pipeline 8 (uncorrelated): (Intercept)    1.853    2.288   2.149   1.560 0.840       0.778 0.659
##                                          SFC    OC   STG insula   ITG
## pipeline 1 (correlated)  : (Intercept) 0.454 0.524 0.724  1.110 0.446
## pipeline 2 (correlated)  : (Intercept) 0.458 0.523 0.725  1.110 0.440
## pipeline 3 (correlated)  : (Intercept) 0.428 0.516 0.704  1.084 0.441
## pipeline 4 (correlated)  : (Intercept) 0.435 0.516 0.708  1.088 0.437
## pipeline 5 (uncorrelated): (Intercept) 0.446 0.476 0.687  1.125 0.407
## pipeline 6 (uncorrelated): (Intercept) 0.380 0.497 0.687  1.113 0.407
## pipeline 7 (uncorrelated): (Intercept) 0.454 0.545 0.721  1.105 0.444
## pipeline 8 (uncorrelated): (Intercept) 0.406 0.441 0.670  1.123 0.387

table.pvalue[1:7];table.pvalue[8:12]
##                                        amygdala thalamus putamen caudate    ACC hippocampus    OFC
## pipeline 1 (correlated)  : (Intercept)    0.386    0.039   0.031   0.038  0.001       0.001  0.199
## pipeline 2 (correlated)  : (Intercept)    0.528     0.02   0.031   0.019  0.001       0.001  0.104
## pipeline 3 (correlated)  : (Intercept)    0.739    0.765   0.231   0.511 <0.001       0.033  0.006
## pipeline 4 (correlated)  : (Intercept)    0.836    0.412   0.180   0.249 <0.001       0.028  0.003
## pipeline 5 (uncorrelated): (Intercept)    0.001   <0.001   0.931   0.007 <0.001      <0.001  0.003
## pipeline 6 (uncorrelated): (Intercept)   <0.001   <0.001   0.986   0.019 <0.001      <0.001  0.001
## pipeline 7 (uncorrelated): (Intercept)        1    0.196   0.175   0.119  0.006      <0.001  0.162
## pipeline 8 (uncorrelated): (Intercept)    0.001   <0.001   0.943   0.003 <0.001      <0.001 <0.001
##                                           SFC     OC    STG insula    ITG
## pipeline 1 (correlated)  : (Intercept)  0.021  0.894  0.904  0.036  0.002
## pipeline 2 (correlated)  : (Intercept)  0.052  0.814   0.92  0.034 <0.001
## pipeline 3 (correlated)  : (Intercept) <0.001  0.306  0.091 <0.001 <0.001
## pipeline 4 (correlated)  : (Intercept) <0.001  0.313  0.156 <0.001 <0.001
## pipeline 5 (uncorrelated): (Intercept)  0.006 <0.001  0.006  0.284 <0.001
## pipeline 6 (uncorrelated): (Intercept)  0.387  0.006  0.011  0.148 <0.001
## pipeline 7 (uncorrelated): (Intercept)  0.095  0.389  0.777  0.045  0.001
## pipeline 8 (uncorrelated): (Intercept) <0.001 <0.001 <0.001  0.238 <0.001


## *** at least one pipeline with no effect
M.interH0 <- do.call(rbind,lapply(ls.mlmmRS, FUN = function(iLMM){
    max(model.tables(iLMM)$p.value)
}))
##                    [,1]
## amygdala    0.971443188
## thalamus    0.514414387
## putamen     0.874083133
## caudate     0.357576367
## ACC         0.003144652
## hippocampus 0.015694795
## OFC         0.137048669
## SFC         0.144185738
## OC          0.602996379
## STG         0.686771867
## insula      0.172066544
## ITG         0.000577909
## PC          0.791305982
## EC          0.000000000

## *** proportion of pipelines with an effect
ls.prop <- lapply(ls.mlmmRS, proportion, method = "single-step", n.sample = 0)
M.prop <- do.call(rbind,ls.prop)[,1:5]
M.prop[,"estimate"] <- paste0(round(100*M.prop[,"estimate"],2),"%")
##              estimate se df lower upper
## amygdala    0.5121319 NA NA    NA    NA
## thalamus    0.6493860 NA NA    NA    NA
## putamen     0.3790945 NA NA    NA    NA
## caudate     0.5958235 NA NA    NA    NA
## ACC         0.9543778 NA NA    NA    NA
## hippocampus 0.8938824 NA NA    NA    NA
## OFC         0.7262501 NA NA    NA    NA
## SFC         0.7421621 NA NA    NA    NA
## OC          0.5313145 NA NA    NA    NA
## STG         0.5233707 NA NA    NA    NA
## insula      0.6107282 NA NA    NA    NA
## ITG         0.9763883 NA NA    NA    NA
## PC          0.5522278 NA NA    NA    NA
## EC          1.0000000 NA NA    NA    NA


## *** common effect
rbind("average" = model.tables(ls.mlmmRS[[1]], method = "average"),
      "pool.se" = model.tables(ls.mlmmRS[[1]], method = "pool.se"),
      "pool.gls" = model.tables(ls.mlmmRS[[1]], method = "pool.gls")

## ** global
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



## * table
df.xtable <- data.frame("proportion of pipelines with a difference" = M.prop[,"estimate"],
                          "p.value against at least one pipeline with no difference" = M.interH0[,1],
                          "p.value against all pipelines with no difference" = apply(do.call(cbind,lapply(ls.SmlmmRS,"[","p.value")),2,min),
                          check.names = FALSE)
xtable::xtable(df.xtable[c("amygdala","thalamus","OFC"),])



df.forest <- do.call(rbind,lapply(1:length(ls.mlmmRS), function(iM){ ## iM <- 1
    iOut <- confint(ls.mlmmRS[[iM]], method = c("none","average","pool.se","pool.gls","pool.gls1"), columns = c("estimate","se","lower","upper","null","p.value"))
    iOut <- cbind(region =  names(ls.mlmmRS)[iM], type = rownames(iOut), iOut)
    rownames(iOut) <- NULL
    iOut
}))
rownames(df.forest) <- NULL

df.cor <- do.call(rbind,lapply(1:length(ls.mlmmRS), function(iM){
    cbind(region = names(ls.mlmmRS)[iM], reshape2::melt(cov2cor(ls.mlmmRS[[iM]]$vcov)))
}))

write.csv(df.cor, file = "results/data-gg-heatmap.csv")
write.csv(df.forest, file = "results/data-gg-forest.csv")


## * figures
reg
## 14 regions: amygdala, thalamus, putamen, caudate, anterior cingulate cortex (ACC),
##             hippocampus, orbital frontal cortex, superior frontal cortex, occipital cortex, superior temporal gyrus,
##             insula, inferior temporal gyrus, parietal cortex, entorhinal cortex


## ** Data
gg.data <- ggplot(dtL, aes(x =  pipeline, y = value)) + facet_wrap(~variable) + geom_point()

## ** Region specific
shape.forest <- c(rep(19,8),c(8,15,17,11))
color.forest <- c(rep("black",8), gg_color_hue(6)[1:4])
    
ls.forestRS <- lapply(ls.mlmmRS, function(eMLMM){ ## eMLMM <- ls.mlmmRS[[1]]

    iGG <- plot(eMLMM, type = "forest", method = c("none","average","pool.se","pool.gls","pool.gls1"), size.text = 11,
                add.args = list(size.estimate = 2.5, size.ci = 1, width.ci = 0.5, shape = shape.forest, color = color.forest), plot = FALSE)$plot
    iYlabel <- ggplot_build(iGG)$layout$panel_params[[1]]$y$get_labels()
    iGG <- iGG + scale_x_discrete(breaks = iYlabel,
                                  labels = gsub("average", "pool (average)",
                                                gsub("pool.se","pool (se)",
                                                     gsub("pool.gls","pool (gls)",
                                                          gsub("pool.gls1","pool (constrained gls)", iYlabel, fixed = TRUE), fixed = TRUE), fixed = TRUE), fixed = TRUE))
    
    iGG + ggtitle(rename.region$full[rename.region$combined==eMLMM$object$outcome])     ## + scale_y_continuous(breaks = seq(0,4,0.25)) 
})[sort(name.region)]

ls.forestRScrop <- lapply(names(ls.forestRS), function(iName){
    iGG <- ls.forestRS[[iName]]
    
    if(iName %in% sort(name.region)[c(1,5,10)] == FALSE){
        iGG <- iGG  + theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() )
    }
    return(iGG + theme(plot.margin=unit(c(0.2,0,-0.4,0),"cm")))
})
ls.forestRScrop <- c(ls.forestRScrop[1:4], list(ggplot()), ls.forestRScrop[-(1:4)])

gg.forestRS <- do.call(egg::ggarrange,c(ls.forestRScrop, list(nrow = 3)))

ls.heatRS.txt <- lapply(ls.mlmmRS, function(eMLMM){ ## eMLMM <- ls.mlmmRS[[1]]
    iPlot <- plot(eMLMM, type = "heat", plot = FALSE, add.args = list(limits = c(0,1.0001), midpoint = 0.5, value.text = TRUE, value.round = 3, value.size = 4))
    iGG <- iPlot$plot + xlab("") + ylab("")
    iGG <- iGG + theme(plot.margin=unit(c(0,0,0,0),"cm"), axis.text.x = element_text(angle = 90, hjust = 1))
    return(iGG)
})[sort(name.region)]

ls.heatRS <- lapply(ls.mlmmRS, function(eMLMM){ ## eMLMM <- ls.mlmmRS[[1]]
    iPlot <- plot(eMLMM, type = "heat", plot = FALSE, add.args = list(limits = c(0,1.0001), midpoint = 0.5), size.text = 12)
    iGG <- iPlot$plot + ggtitle(rename.region$full[rename.region$combined==eMLMM$object$outcome]) + xlab("") + ylab("")
    iGG <- iGG + theme(plot.margin=unit(c(0,0,0,0),"cm"), axis.text.x = element_text(angle = 90, hjust = 1))
    return(iGG)
})[sort(name.region)]

ls.heatRScrop <- lapply(names(ls.heatRS), function(iName){
    iGG <- ls.heatRS[[iName]] + guides(fill = "none")
    
    if(iName %in% sort(name.region)[c(1,5,10)] == FALSE){
        iGG <- iGG + theme(axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title.y = element_blank() )
    }    
    if(iName %in% sort(name.region)[c(10:14)] == FALSE){
        iGG <- iGG + theme(axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.title.x = element_blank() )
    }    
    return(iGG)
})
ls.heatRScrop <- c(ls.heatRScrop[1:4], list(ggpubr::as_ggplot(ggpubr::get_legend(ls.heatRS[[1]]))), ls.heatRScrop[-(1:4)])

gg.heatRS <- do.call(egg::ggarrange, c(ls.heatRScrop, list(nrow = 3, ncol = 5)))

## export
pdf(file.path("figures","application-figure-data.pdf"), width = 12)
gg.data
dev.off()

png(file.path("figures","application-figure-data.png"), width = 1000, heigh = 650)
gg.data
dev.off()

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

for(iR in names(ls.mlmmRS)){ ## iR <- names(ls.mlmmRS)[1]
    pdf(file.path("figures",paste0("application-figure-",iR,".pdf")), width = 12)
    XXX <- ggarrange(ls.forestRS[[iR]] + theme(text = element_text(size=15)),ls.heatRS.txt[[iR]], ncol = 2, widths = c(0.4,0.6), newpage = FALSE)
    dev.off()
}


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
