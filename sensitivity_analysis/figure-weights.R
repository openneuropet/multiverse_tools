### FIGURE-weights.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct 11 2022 (12:03) 
## Version: 
## Last-Updated: okt 12 2022 (11:51) 
##           By: Brice Ozenne
##     Update #: 9
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(ggplot2)

wAverage <- rep(5,20)
wPool <- c(rep(3.81,15),38.1,1.9,1.2,0.95,0.63)
wGLS <- c(rep(0.57,15),81.31,4.07,2.71,2.03,1.36)

name.pip <- c(paste0("C",1:9,"   (\u03C3\u00B2=2.5)"),paste0("C",10:15," (\u03C3\u00B2=2.5)"),paste0("I",1:5,"    (\u03C3\u00B2=",c(0.25,5,7,10,15),")"))
df <- rbind(data.frame(pipeline = name.pip,
                       estimator = "average",
                       weight = wAverage/100),
            data.frame(pipeline = name.pip,
                       estimator = "pool-se",
                       weight = wPool/100),
            data.frame(pipeline = name.pip,
                       estimator = "GLS",
                       weight = wGLS/100))
df$pipeline <- factor(df$pipeline, rev(name.pip))
df$estimator <- factor(df$estimator, rev(unique(df$estimator)))

gg.bar <- ggplot(df, aes(fill=pipeline, y=weight, x=estimator))
gg.bar <- gg.bar + geom_bar(position="stack", stat="identity")
gg.bar <- gg.bar + scale_fill_manual(values = c(rainbow(5),gray.colors(15)))
gg.bar <- gg.bar + coord_flip() + scale_y_continuous(labels = scales::percent)
gg.bar <- gg.bar + ggtitle("Scenario 3")
gg.bar <- gg.bar + theme(text = element_text(size=15), 
                         legend.text = element_text(size=9), 
                         axis.line.x = element_line(size = 1.25),
                         axis.title = element_blank(),
                         axis.ticks = element_line(size = 2),
                         axis.ticks.length=unit(.25, "cm"),
                         panel.grid.major=element_blank(),
                         panel.grid.minor=element_blank(),
                         plot.background=element_blank(),
                         panel.background=element_blank(),
                         panel.border=element_blank(),
                         legend.key.size = unit(0.75,"line"))
gg.bar

ggsave(gg.bar, filename = "figures/gg-weights-scenario3.pdf", height =  4, width = 10, device = cairo_pdf)
ggsave(gg.bar, filename = "figures/gg-weights-scenario3.png", height =  4, width = 10, device = cairo_pdf)


##----------------------------------------------------------------------
### FIGURE-weights.R ends here
