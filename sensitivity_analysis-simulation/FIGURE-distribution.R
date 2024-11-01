### FIGURE-distribution.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  1 2024 (10:17) 
## Version: 
## Last-Updated: nov  1 2024 (11:34) 
##           By: Brice Ozenne
##     Update #: 4
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(Matrix)
library(ggplot2)
library(ggpubr)
source("FCT_simulation.R")

Sigma <- as.matrix(bdiag((diag(0.05,15,15) + 0.95)*2.5, matrix(0.25), matrix(5), matrix(7.5), matrix(10), matrix(15)))
Sigmatot <- Sigma+diag(1,NROW(Sigma),NROW(Sigma))
n.obs <- 1e5

set.seed(1)
dt.scenario3 <- simData(n.obs = n.obs, sigma.pipe = Sigma, beta = 0)$wide
dt.scenario4 <- simData(n.obs = n.obs, sigma.pipe = Sigma, beta = 0, df = 3)$wide
dt.scenario5 <- simData(n.obs = n.obs, sigma.pipe = Sigma, beta = 0, half.distribution = TRUE)$wide

dt.dist <- rbind(data.table(scenario = "Scenario 3", type = "Noise", pipeline = "Pipeline 1", value = dt.scenario3$Epip.1),
                 data.table(scenario = "Scenario 3", type = "Observed", pipeline = "Pipeline 1", value = dt.scenario3$Ypip.1),
                 data.table(scenario = "Scenario 3", type = "Noise", pipeline = "Pipeline 16", value = dt.scenario3$Epip.16),
                 data.table(scenario = "Scenario 3", type = "Observed", pipeline = "Pipeline 16", value = dt.scenario3$Ypip.16),
                 data.table(scenario = "Scenario 3", type = "Noise", pipeline = "Pipeline 20", value = dt.scenario3$Epip.20),
                 data.table(scenario = "Scenario 3", type = "Observed", pipeline = "Pipeline 20", value = dt.scenario3$Ypip.20),
                 data.table(scenario = "Scenario 4", type = "Noise", pipeline = "Pipeline 1", value = dt.scenario4$Epip.1),
                 data.table(scenario = "Scenario 4", type = "Observed", pipeline = "Pipeline 1", value = dt.scenario4$Ypip.1),
                 data.table(scenario = "Scenario 4", type = "Noise", pipeline = "Pipeline 16", value = dt.scenario4$Epip.16),
                 data.table(scenario = "Scenario 4", type = "Observed", pipeline = "Pipeline 16", value = dt.scenario4$Ypip.16),
                 data.table(scenario = "Scenario 4", type = "Noise", pipeline = "Pipeline 20", value = dt.scenario4$Epip.20),
                 data.table(scenario = "Scenario 4", type = "Observed", pipeline = "Pipeline 20", value = dt.scenario4$Ypip.20),
                 data.table(scenario = "Scenario 5", type = "Noise", pipeline = "Pipeline 1", value = dt.scenario5$Epip.1),
                 data.table(scenario = "Scenario 5", type = "Observed", pipeline = "Pipeline 1", value = dt.scenario5$Ypip.1),
                 data.table(scenario = "Scenario 5", type = "Noise", pipeline = "Pipeline 16", value = dt.scenario5$Epip.16),
                 data.table(scenario = "Scenario 5", type = "Observed", pipeline = "Pipeline 16", value = dt.scenario5$Ypip.16),
                 data.table(scenario = "Scenario 5", type = "Noise", pipeline = "Pipeline 20", value = dt.scenario5$Epip.20),
                 data.table(scenario = "Scenario 5", type = "Observed", pipeline = "Pipeline 20", value = dt.scenario5$Ypip.20))
dt.dist[, pipeline := as.factor(pipeline)]          

dt.dist[, scenario2 := factor(scenario, levels = paste0("Scenario ",3:5), labels = c("Scenario 3 \n (Gaussian)","Scenario 4 \n (Student)","Scenario 5 \n (half Gaussian)"))]
dt.dist[, pipeline2 := factor(pipeline, levels = paste0("Pipeline ",c(1,16,20)), labels = c("Pipeline 1 \n (correlated, \n moderate noise)","Pipeline 16 \n (low noise)","Pipeline 20 \n (high noise)"))]

gg.dist <- ggplot(dt.dist, aes(x=value, group = scenario2, color = scenario2))
gg.dist <- gg.dist + facet_grid(pipeline2~type, scales = "free")
gg.dist <- gg.dist + theme(#plot.margin = unit(c(0,0.2,0,1), 'lines'),
                           text = element_text(size=15), 
                           legend.text = element_text(size=9), 
                           axis.line.x = element_line(linewidth = 1.25),
                           axis.title = element_blank(),
                           axis.ticks = element_line(size = 2),
                           axis.ticks.length=unit(.25, "cm"),
                           panel.grid.major=element_blank(),
                           panel.grid.minor=element_blank(),
                           plot.background=element_blank(),
                           panel.background=element_blank(),
                           panel.border=element_blank(),
                           legend.key.size = unit(0.75,"line"))
gg.dist <- gg.dist + labs(x = "", y = "Density", color = "", fill = "")

gg.arrdist <- ggarrange(gg.dist %+% dt.dist[pipeline == "Pipeline 16"] + coord_cartesian(xlim = c(-3,3)) + geom_density(alpha = 0.2, linewidth = 1.25, bw = 0.01),
                        gg.dist %+% dt.dist[pipeline == "Pipeline 1"] + coord_cartesian(xlim = c(-5,5)) + geom_density(alpha = 0.2, linewidth = 1.25, bw = 0.01),
                        gg.dist %+% dt.dist[pipeline == "Pipeline 20"] + coord_cartesian(xlim = c(-10,10)) + geom_density(alpha = 0.2, linewidth = 1.25, bw = 0.15),
                        ncol = 1, nrow = 3, legend = "bottom", common.legend = TRUE)
gg.arrdist

ggsave(gg.arrdist, filename = "figures/figure-simulation-distribution.pdf", height =  8, width = 10, device = cairo_pdf)




##----------------------------------------------------------------------
### FIGURE-distribution.R ends here
