### FIGURE-simulation.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 12 2022 (11:48) 
## Version: 
## Last-Updated: nov  1 2024 (11:34) 
##           By: Brice Ozenne
##     Update #: 63
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Dependencies
library(data.table)
library(ggplot2)
library(scales)
library(gridExtra)
library(ggpubr)

## * Data
dtS.sim <- readRDS(file.path("Results","simulation-summary-scenario.rds"))

## ** uncertainty about estimated mean and proportion
dtS.sim[target=="pool",average.upper := average + qnorm(0.975) * sd/sqrt(rep)]
dtS.sim[target=="pool",average.lower := average + qnorm(0.025) * sd/sqrt(rep)]
dtS.sim[target=="proportion",average.upper := average + qnorm(0.975) * sqrt(average*(1-average))/sqrt(rep)]
dtS.sim[target=="proportion",average.lower := average + qnorm(0.025) * sqrt(average*(1-average))/sqrt(rep)]
dtS.sim[,power.upper := power + qnorm(0.975) * sqrt(power*(1-power))/sqrt(rep)]
dtS.sim[,power.lower := power + qnorm(0.025) * sqrt(power*(1-power))/sqrt(rep)]

## UNCERTAINTY ABOUT MEAN
## library(LMMstar)
## X <- rnorm(1e4)
## vcov(lmm(Y ~ 1, data = data.frame(Y=X)))
## var(X)/sqrt(1e4)

## ** uncertainty about estimated sd
dtS.sim[,sd.upper := exp(log(sd) + qnorm(0.975)/sqrt(2*(rep-1)))]
dtS.sim[,sd.lower := exp(log(sd) + qnorm(0.025)/sqrt(2*(rep-1)))]

## library(LMMstar)
## X <- rnorm(1e4)
## vcov(lmm(Y ~ 1, data = data.frame(Y=X)), effects = "all", transform.sigma = "log")
## 1/(2*(NROW(X)-1))

## vcov(lmm(Y ~ 1, data = data.frame(Y=X)), effects = "all", transform.sigma = "square")
## 2*var(X)^2/(NROW(X)-1)

## vcov(lmm(Y ~ 1, data = data.frame(Y=X)), effects = "all", transform.sigma = "none")
## var(X)/(2*(NROW(X)-1))
## sd(X)/sqrt(2*(NROW(X)-1))
## 1/sqrt(2*(NROW(X)-1))

## * Figure 2

dtSpool.sim <- dtS.sim[target=="pool" & type %in% c("average", "fixse", "gls", "gls1")]

## ** bias
ggPool.bias <- ggplot(dtSpool.sim, aes(x = n.char, y = average-GS, linetype = target, group = type.char)) 
ggPool.bias <- ggPool.bias + geom_ribbon(aes(ymin = average.lower-GS, ymax = average.upper-GS), alpha = 0.25)
ggPool.bias <- ggPool.bias + geom_hline(aes(yintercept = 0), color = "brown") + facet_grid(beta.char~scenario) 
ggPool.bias <- ggPool.bias + geom_point(size = 2, aes(shape = type.char, color = type.char)) + geom_line(linewidth = 1, aes(color = type.char))
ggPool.bias <- ggPool.bias + labs(x = "Sample size (per group)", y = "Estimate", shape = "", color = "") + guides(linetype = "none")
ggPool.bias <- ggPool.bias + scale_color_manual(values = unname(palette.colors()[c(7,5,4,3)]))
ggPool.bias <- ggPool.bias + scale_color_manual(values = unname(palette.colors()[c(7,5,4,3)]),
                                              breaks = paste0("pool (",c("average","se","gls","constrained gls"),")"),
                                              labels = c(expression(hat(Psi)["average"]),expression(hat(Psi)["pool-se"]),expression(hat(Psi)["GLS"]),expression(hat(Psi)["constrained GLS"])))
ggPool.bias <- ggPool.bias + scale_shape_manual(values = c(15,17,19,8),
                                            breaks = paste0("pool (",c("average","se","gls","constrained gls"),")"),
                                            labels = c(expression(hat(Psi)["average"]),expression(hat(Psi)["pool-se"]),expression(hat(Psi)["GLS"]),expression(hat(Psi)["constrained GLS"])))
ggPool.bias <- ggPool.bias + theme_minimal() + theme(legend.position="bottom",
                                                     text = element_text(size=15), 
                                                     axis.line = element_line(linewidth = 1.25),
                                                     axis.ticks = element_line(linewidth = 2),
                                                     axis.ticks.length=unit(.25, "cm"),
                                                     legend.key.size = unit(3,"line"))
ggPool.bias + coord_cartesian(ylim = c(-0.02,0.02))
ggPool.bias


## ** se
ggPool.se <- ggplot(dtSpool.sim, aes(x = n.char, y = sd, linetype = target, group = type.char)) 
ggPool.se <- ggPool.se + facet_grid(beta.char~scenario) 
ggPool.se <- ggPool.se + geom_ribbon(aes(ymin = sd.lower, ymax = sd.upper), alpha = 0.25)
ggPool.se <- ggPool.se + geom_point(size = 2.5, aes(shape = type.char, color = type.char)) + geom_line(linewidth = 1.2, aes(color = type.char)) 
ggPool.se <- ggPool.se + labs(x = "Sample size (per group)", y = "Standard deviation", shape = "Global effect estimator", color = "Global effect estimator") + guides(linetype = "none")
ggPool.se <- ggPool.se + scale_color_manual(values = unname(palette.colors()[c(7,5,4,3)]),
                                            breaks = paste0("pool (",c("average","se","gls","constrained gls"),")"),
                                            labels = c(expression(hat(Psi)["average"]),expression(hat(Psi)["pool-se"]),expression(hat(Psi)["GLS"]),expression(hat(Psi)["constrained GLS"])))
ggPool.se <- ggPool.se + scale_shape_manual(values = c(15,17,19,8),
                                            breaks = paste0("pool (",c("average","se","gls","constrained gls"),")"),
                                            labels = c(expression(hat(Psi)["average"]),expression(hat(Psi)["pool-se"]),expression(hat(Psi)["GLS"]),expression(hat(Psi)["constrained GLS"])))
ggPool.se <- ggPool.se + theme_minimal() + theme(legend.position="bottom",
                                                 text = element_text(size=15), 
                                                 axis.line = element_line(linewidth = 1.25),
                                                 axis.ticks = element_line(size = 2),
                                                 axis.ticks.length=unit(.25, "cm"),
                                                 legend.key.size = unit(3,"line"))
ggPool.se


ggPoolLog.se <- ggPool.se + coord_trans(y="log10") + ylab("Standard deviation (log scale)")
ggPoolLog.se

## under the null
ggPoolLog.seH0 <- ggPoolLog.se %+% dtSpool.sim[beta == 0 & scenario %in% paste0("scenario ",1:3),]
ggPoolLog.seH0 <- ggPoolLog.seH0 + scale_y_continuous(breaks = c(0.1,0.5,1,1.5), labels = paste0("     ",c(0.1,0.5,1,1.5)))
ggPoolLog.seH0

## under the alternative
ggPoolLog.seH1 <- ggPoolLog.se %+% dtSpool.sim[beta == 0.5 & scenario %in% paste0("scenario ",1:3),]
ggPoolLog.seH1 <- ggPoolLog.seH1 + scale_y_continuous(breaks = c(0.1,0.5,1,1.5), labels = paste0("     ",c(0.1,0.5,1,1.5)))
ggPoolLog.seH1

## ** rejection rate
ggPool.rejection <- ggplot(dtSpool.sim, aes(x = n.char, y = power, linetype = target, group = type.char)) 
ggPool.rejection <- ggPool.rejection + facet_grid(beta.char~scenario) 
ggPool.rejection <- ggPool.rejection + geom_ribbon(aes(ymin = power.lower, ymax = power.upper), alpha = 0.1)
ggPool.rejection <- ggPool.rejection + geom_hline(yintercept = 0.05, color = "brown")
ggPool.rejection <- ggPool.rejection + geom_point(size = 2.5, aes(shape = type.char, color = type.char)) + geom_line(linewidth = 1.2, aes(color = type.char))
ggPool.rejection <- ggPool.rejection + labs(x = "Sample size (per group)", y = "Rejection rate", color = "Global effect estimator", shape = "Global effect estimator") + guides(linetype = "none")
ggPool.rejection <- ggPool.rejection + scale_color_manual(values = unname(palette.colors()[c(7,5,4,3)]),
                                            breaks = paste0("pool (",c("average","se","gls","constrained gls"),")"),
                                            labels = c(expression(hat(Psi)["average"]),expression(hat(Psi)["pool-se"]),expression(hat(Psi)["GLS"]),expression(hat(Psi)["constrained GLS"])))
ggPool.rejection <- ggPool.rejection + scale_shape_manual(values = c(15,17,19,8),
                                            breaks = paste0("pool (",c("average","se","gls","constrained gls"),")"),
                                            labels = c(expression(hat(Psi)["average"]),expression(hat(Psi)["pool-se"]),expression(hat(Psi)["GLS"]),expression(hat(Psi)["constrained GLS"])))
ggPool.rejection <- ggPool.rejection + scale_y_continuous(labels = scales::percent)
ggPool.rejection <- ggPool.rejection + theme_minimal() + theme(legend.position="bottom",
                                                               text = element_text(size=15), 
                                                               axis.line = element_line(linewidth = 1.25),
                                                               axis.ticks = element_line(size = 2),
                                                               axis.ticks.length=unit(.25, "cm"),
                                                               legend.key.size = unit(3,"line"))
ggPool.rejection

ggPoolLog.rejection <- ggPool.rejection + coord_trans(y="log10") + ylab("Rejection rate (log scale)")
ggPoolLog.rejection

## under the null
ggPoolLog.type1 <- ggPoolLog.rejection %+% dtSpool.sim[beta == 0 & scenario %in% paste0("scenario ",1:3),]
ggPoolLog.type1 <- ggPoolLog.type1 + scale_y_continuous(breaks = c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),
                                                labels = scales::percent)
ggPoolLog.type1

## under the alternative
ggPoolLog.power <- ggPoolLog.rejection %+% dtSpool.sim[beta == 0.5 & scenario %in% paste0("scenario ",1:3),]
ggPoolLog.power <- ggPoolLog.power + scale_y_continuous(breaks = c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),
                                                labels = scales::percent)
ggPoolLog.power

## ** combine
figure2 <- ggarrange(ggPoolLog.seH0 + labs(x = ""),
                     ggPoolLog.type1 + theme(strip.background = element_blank(), strip.text.x = element_blank()),
                     nrow = 2, common.legend = TRUE, legend = "bottom")
figure2

figure2.bis <- ggarrange(ggPoolLog.seH1 + labs(x = ""),
                         ggPoolLog.power + theme(strip.background = element_blank(),strip.text.x = element_blank()),
                         nrow = 2, common.legend = TRUE, legend = "bottom")
figure2.bis

## ** export
cairo_pdf(file.path("figures","figure-simulation-pool-H0.pdf"), width = 10, height = 7)
figure2
dev.off()

png(file.path("figures","figure-simulation-pool-H0.png"), width = 1000, height = 700)
figure2
dev.off()

cairo_pdf(file.path("figures","figure-simulation-pool-H1.pdf"), width = 10, height = 7)
figure2.bis
dev.off()

png(file.path("figures","figure-simulation-pool-H1.png"), width = 1000, height = 700)
figure2.bis
dev.off()

## * Figure 3

dtSprop.sim <- dtS.sim[target=="proportion"]
dtSprop.sim[,beta := as.numeric(beta>0)]

## ** bias
ggProp.bias <- ggplot(dtSprop.sim, aes(x = n.char, y = average, linetype = target, group = type.char)) 
ggProp.bias <- ggProp.bias + geom_ribbon(aes(ymin = average.lower, ymax = average.upper), alpha = 0.25)
ggProp.bias <- ggProp.bias + geom_hline(aes(yintercept = pmax(0.05,beta)), color = "brown") + facet_grid(beta.char~scenario, scales = "free") 
ggProp.bias <- ggProp.bias + geom_point(size = 4, aes(shape = type.char, color = type.char)) + geom_line(linewidth = 1, aes(color = type.char))
ggProp.bias <- ggProp.bias + labs(x = "Sample size (per group)", y = "Average estimate", shape = "Proportion", color = "Proportion") + guides(linetype = "none")
ggProp.bias <- ggProp.bias + scale_color_manual(values = unname(palette.colors()[c(2,6)]), breaks = c("proportion (np)","proportion (p)"),
                                                labels = c(expression(hat(eta)[1]),expression(hat(eta)[Phi])))
ggProp.bias <- ggProp.bias + scale_shape_manual(values = c(18,4), breaks = c("proportion (np)","proportion (p)"),
                                                labels = c(expression(hat(eta)[1]),expression(hat(eta)[Phi])))
ggProp.bias <- ggProp.bias + theme_minimal() + theme(legend.position="bottom",
                                                     text = element_text(size=15), 
                                                     axis.line = element_line(linewidth = 1.25),
                                                     axis.ticks = element_line(size = 2),
                                                     axis.ticks.length=unit(.25, "cm"),
                                                     legend.key.size = unit(3,"line"))
ggProp.bias
## ggProp.bias + coord_cartesian(ylim = c(-0.01,0.01))

## under the null
ggProp.biasH0 <- ggProp.bias %+% dtSprop.sim[beta == 0 & scenario %in% paste0("scenario ",1:3),]
ggProp.biasH0

## under the alternative
ggProp.biasH1 <- ggProp.bias %+% dtSprop.sim[beta > 0 & scenario %in% paste0("scenario ",1:3),]
ggProp.biasH1

## ** se
ggProp.se <- ggplot(dtSprop.sim, aes(x = n.char, y = sd, linetype = target, group = type.char)) 
ggProp.se <- ggProp.se + facet_grid(beta.char~scenario) 
ggProp.se <- ggProp.se + geom_ribbon(aes(ymin = sd.lower, ymax = sd.upper), alpha = 0.25)
ggProp.se <- ggProp.se + geom_point(size = 2.5, aes(shape = type.char, color = type.char)) + geom_line(linewidth = 1.2, aes(color = type.char)) 
ggProp.se <- ggProp.se + labs(x = "Sample size (per group)", y = "Standard deviation", shape = "Proportion", color = "Proportion") + guides(linetype = "none")
ggProp.se <- ggProp.se + scale_color_manual(values = unname(palette.colors()[c(2,6)]), breaks = c("proportion (np)","proportion (p)"),
                                            labels = c(expression(hat(eta)[1]),expression(hat(eta)[Phi])))
ggProp.se <- ggProp.se + scale_shape_manual(values = c(18,4), breaks = c("proportion (np)","proportion (p)"),
                                            labels = c(expression(hat(eta)[1]),expression(hat(eta)[Phi])))
ggProp.se <- ggProp.se + theme_minimal() + theme(legend.position="bottom",
                                                 text = element_text(size=15), 
                                                 axis.line = element_line(size = 1.25),
                                                 axis.ticks = element_line(size = 2),
                                                 axis.ticks.length=unit(.25, "cm"),
                                                 legend.key.size = unit(3,"line"))
ggProp.se


ggPropLog.se <- ggProp.se + coord_trans(y="log10") + ylab("Standard deviation (log scale)")
ggPropLog.se

## under the null
ggPropLog.seH0 <- ggPropLog.se %+% dtSprop.sim[beta == 0 & scenario %in% paste0("scenario ",1:3),]
ggPropLog.seH0 <- ggPropLog.seH0 + scale_y_continuous(breaks = c(0.10,0.08,0.06,0.04), labels = paste0("     ",c(0.10,0.08,0.06,0.04)))
ggPropLog.seH0

## under the alternative
ggPropLog.seH1 <- ggPropLog.se %+% dtSprop.sim[beta > 0 & scenario %in% paste0("scenario ",1:3),]
ggPropLog.seH1 <- ggPropLog.seH1 + scale_y_continuous(breaks = c(0.4,0.3,0.2,0.1,0.05), labels = paste0("     ",c(0.4,0.3,0.2,0.1,0.05)))
ggPropLog.seH1

## ** rejection rate
## ggProp.rejection <- ggplot(dtSprop.sim, aes(x = n.char, y = power, linetype = target, group = type.char)) 
## ggProp.rejection <- ggProp.rejection + facet_grid(beta.char~scenario) 
## ggProp.rejection <- ggProp.rejection + geom_ribbon(aes(ymin = power.lower, ymax = power.upper), alpha = 0.1)
## ggProp.rejection <- ggProp.rejection + geom_hline(yintercept = 0.05, color = "brown")
## ggProp.rejection <- ggProp.rejection + geom_point(size = 2.5, aes(shape = type.char, color = type.char)) + geom_line(linewidth = 1.2, aes(color = type.char))
## ggProp.rejection <- ggProp.rejection + labs(x = "Sample size (per group)", y = "Rejection rate", color = "", shape = "") + guides(linetype = "none")
## ggProp.rejection <- ggProp.rejection + scale_color_manual(values = unname(palette.colors()[c(2,6)]), breaks = c("proportion (np)","proportion (p)"),
##                                                           labels = c(expression(hat(eta)[1]),expression(hat(eta)[Phi])))
## ggProp.rejection <- ggProp.rejection + scale_shape_manual(values = c(18,4), breaks = c("proportion (np)","proportion (p)"),
##                                                           labels = c(expression(hat(eta)[1]),expression(hat(eta)[Phi])))
## ggProp.rejection <- ggProp.rejection + scale_y_continuous(labels = scales::percent)
## ggProp.rejection <- ggProp.rejection + theme_minimal() + theme(legend.position="bottom",
##                                                                text = element_text(size=15), 
##                                                                axis.line = element_line(linewidth = 1.25),
##                                                                axis.ticks = element_line(size = 2),
##                                                                axis.ticks.length=unit(.25, "cm"),
##                                                                legend.key.size = unit(3,"line"))
## ggProp.rejection

## ggPropLog.rejection <- ggProp.rejection + coord_trans(y="log10") + ylab("Rejection rate (log scale)")
## ggPropLog.rejection

## ## under the null
## ggPropLog.type1 <- ggPropLog.rejection %+% dtSprop.sim[beta == 0,]
## ggPropLog.type1 <- ggPropLog.type1 + scale_y_continuous(breaks = c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),
##                                                 labels = scales::percent)
## ggPropLog.type1

## ## under the alternative
## ggPropLog.power <- ggPropLog.rejection %+% dtSprop.sim[beta == 0.5,]
## ggPropLog.power <- ggPropLog.power + scale_y_continuous(breaks = c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),
##                                                 labels = scales::percent)
## ggPropLog.power

## ** combine
figure3 <- ggarrange(ggProp.biasH0 + labs(x = ""),
                     ggPropLog.seH0 + theme(strip.background = element_blank(), strip.text.x = element_blank()),
                     nrow = 2, common.legend = TRUE, legend = "bottom")
figure3

figure3.bis <- ggarrange(ggProp.biasH1 + labs(x = ""),
                         ggPropLog.seH1 + theme(strip.background = element_blank(), strip.text.x = element_blank()),
                         nrow = 2, common.legend = TRUE, legend = "bottom")
figure3.bis


## ** export
cairo_pdf(file.path("figures","figure-simulation-proportion-H0.pdf"), width = 10, height = 7)
figure3
dev.off()

png(file.path("figures","figure-simulation-proportion-H0.png"), width = 1000, height = 700)
figure3
dev.off()

cairo_pdf(file.path("figures","figure-simulation-proportion-H1.pdf"), width = 10, height = 7)
figure3.bis
dev.off()

png(file.path("figures","figure-simulation-proportion-H1.png"), width = 1000, height = 700)
figure3.bis
dev.off()

## * Figure C & D

pdf(file.path("figures","figure-simulation-pool-bias.pdf"), width = 10, height = 7)
ggPool.bias + coord_cartesian(ylim = c(-0.02,0.03)) + ylab("Bias")
dev.off()

xxx <- ggarrange(ggPool.se %+% dtSpool.sim[scenario %in% paste0("scenario ",4:5)], ggPool.rejection %+% dtSpool.sim[scenario %in% paste0("scenario ",4:5)], nrow = 1,
                 legend = "bottom", common.legend = TRUE)
pdf(file.path("figures","figure-simulation-pool-scenario4_5.pdf"), width = 10, height = 7)
xxx
dev.off()


yyy <- ggarrange(ggProp.bias %+% dtSprop.sim[scenario %in% paste0("scenario ",4:5)], ggProp.se %+% dtSprop.sim[scenario %in% paste0("scenario ",4:5)], nrow = 1,
                 legend = "bottom", common.legend = TRUE)
pdf(file.path("figures","figure-simulation-proportion-scenario4_5.pdf"), width = 10, height = 7)
yyy
dev.off()



##----------------------------------------------------------------------
### FIGURE-simulation.R ends here
