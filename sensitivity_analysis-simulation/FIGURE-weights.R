### FIGURE-weights.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct 11 2022 (12:03) 
## Version: 
## Last-Updated: okt 28 2024 (17:38) 
##           By: Brice Ozenne
##     Update #: 12
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(ggplot2)
library(Matrix)


Sigma3 <- as.matrix(bdiag((diag(0.05,15,15) + 0.95)*2.5, matrix(0.25), matrix(5), matrix(7.5), matrix(10), matrix(15)))
Sigmatot3 <- Sigma3+1

## average
weight3.average <- 1/NROW(Sigmatot3)
## pool-se
weight3.poolse <- (1/diag(Sigmatot3))/sum((1/diag(Sigmatot3)))
## gls
weight3.gls <- rowSums(solve(Sigmatot3))/sum(solve(Sigmatot3))

name.pip <- c(paste0("C",1:9,"   (\u03C3\u00B2=2.5)"),paste0("C",10:15," (\u03C3\u00B2=2.5)"),paste0("I",1:5,"    (\u03C3\u00B2=",c(0.25,5,7,10,15),")"))
df <- rbind(data.frame(pipeline = name.pip,
                       estimator = "average",
                       weight = weight3.average/100),
            data.frame(pipeline = name.pip,
                       estimator = "pool-se",
                       weight = weight3.poolse/100),
            data.frame(pipeline = name.pip,
                       estimator = "GLS",
                       weight = weight3.gls/100))
df$pipeline <- factor(df$pipeline, rev(name.pip))
df$estimator <- factor(df$estimator, rev(unique(df$estimator)))

gg.bar <- ggplot(df, aes(fill=pipeline, y=weight, x=estimator))
gg.bar <- gg.bar + geom_bar(position="stack", stat="identity")
gg.bar <- gg.bar + scale_fill_manual(values = c(rainbow(5),gray.colors(15)))
gg.bar <- gg.bar + coord_flip() + scale_y_continuous(labels = scales::percent)
gg.bar <- gg.bar + scale_x_discrete(breaks = c("average","pool-se","GLS"),
                                    labels = c(expression(hat(Psi)["average"]),expression(hat(Psi)["pool-se"]),expression(hat(Psi)["GLS"])))
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

ggsave(gg.bar, filename = "figures/figure-simulation-weights-scenario3.pdf", height =  4, width = 10, device = cairo_pdf)
ggsave(gg.bar, filename = "figures/figure-simulation-weights-scenario3.png", height =  4, width = 10, device = cairo_pdf)


##----------------------------------------------------------------------
### FIGURE-weights.R ends here
