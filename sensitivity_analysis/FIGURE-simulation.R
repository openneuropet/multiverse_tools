### FIGURE-simulation.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 12 2022 (11:48) 
## Version: 
## Last-Updated: Dec  5 2022 (13:12) 
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

## * Dependencies
library(data.table)
library(ggplot2)
library(scales)
library(gridExtra)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## * Data
dtS.sim <- readRDS(file.path("results","simulation-summary-scenario.rds"))

## * create table and figure
## ** bias
dtS.sim[target=="pool",average.upper := average + qnorm(0.975) * sd/sqrt(rep)]
dtS.sim[target=="pool",average.lower := average + qnorm(0.025) * sd/sqrt(rep)]
dtS.sim[target=="proportion",average.upper := average + qnorm(0.975) * sqrt(average*(1-average))/sqrt(rep)]
dtS.sim[target=="proportion",average.lower := average + qnorm(0.025) * sqrt(average*(1-average))/sqrt(rep)]

gg.bias <- ggplot(dtS.sim, aes(x = n.char, y = average, linetype = target, group = type.char)) 
gg.bias <- gg.bias + geom_ribbon(aes(ymin = average.lower, ymax = average.upper), alpha = 0.25)
gg.bias <- gg.bias + geom_hline(aes(yintercept = beta), color = "brown") + facet_grid(beta.char~scenario) 
gg.bias <- gg.bias + geom_point(size = 2, aes(shape = type.char, color = type.char)) + geom_line(size = 1, aes(color = type.char))
gg.bias <- gg.bias + labs(x = "Sample size (per group)", y = "Estimate", shape = "", color = "") + guides(linetype = "none")
gg.bias <- gg.bias + scale_color_manual(values = gg_color_hue(6))
gg.bias <- gg.bias + theme_minimal() + theme(legend.position="bottom",
                                             text = element_text(size=15), 
                                             axis.line = element_line(size = 1.25),
                                             axis.ticks = element_line(size = 2),
                                             axis.ticks.length=unit(.25, "cm"),
                                             legend.key.size = unit(3,"line"))
## gg.bias
## gg.bias + coord_cartesian(ylim = c(-0.01,0.01))

## UNCERTAINTY ABOUT MEAN
## library(LMMstar)
## X <- rnorm(1e4)
## vcov(lmm(Y ~ 1, data = data.frame(Y=X)))
## var(X)/sqrt(1e4)


## ** se
dtS.sim[,sd.upper := exp(log(sd) + qnorm(0.975)/sqrt(2*(rep-1)))]
dtS.sim[,sd.lower := exp(log(sd) + qnorm(0.025)/sqrt(2*(rep-1)))]

range(dtS.sim[beta==0, sd]-dtS.sim[beta==0.5, sd])
range((dtS.sim[beta==0, sd]-dtS.sim[beta==0.5, sd])/abs(dtS.sim[beta==0, sd]))

gg.se <- ggplot(dtS.sim, aes(x = n.char, y = sd, linetype = target, group = type.char)) 
gg.se <- gg.se + facet_grid(beta.char~scenario) 
gg.se <- gg.se + geom_ribbon(aes(ymin = sd.lower, ymax = sd.upper), alpha = 0.25)
gg.se <- gg.se + geom_point(size = 2, aes(shape = type.char, color = type.char)) + geom_line(size = 1, aes(color = type.char)) 
gg.se <- gg.se + labs(x = "Sample size (per group)", y = "Standard deviation", shape = "", color = "") + guides(linetype = "none")
gg.se <- gg.se + scale_color_manual(values = gg_color_hue(6))
gg.se <- gg.se + theme_minimal() + theme(legend.position="bottom",
                                         text = element_text(size=15), 
                                         axis.line = element_line(size = 1.25),
                                         axis.ticks = element_line(size = 2),
                                         axis.ticks.length=unit(.25, "cm"),
                                         legend.key.size = unit(3,"line"))
## gg.se

## UNCERTAINTY ABOUT SE
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

gg.seLog <- gg.se + coord_trans(y="log10", ylim = c(0.05,1.5)) + ylab("Standard deviation (log scale)")
gg.seLog <- gg.seLog + scale_y_continuous(breaks = c(1.5,1,0.5,0.1))

gg.se4 <- gg.se %+% dtS.sim[type.char %in% c("proportion","pool (robust gls)") == FALSE,]
gg.seLog4 <- gg.seLog %+% dtS.sim[type.char %in% c("proportion","pool (robust gls)") == FALSE,]
gg.seLog4.H0 <- gg.seLog %+% dtS.sim[type.char %in% c("proportion","pool (robust gls)") == FALSE & beta == 0,]

## ** cali se
gg.cali <- ggplot(dtS.sim[!is.na(dtS.sim$average.se)], aes(x = sd, y = average.se, linetype = target, color = type.char, group = type.char)) 
gg.cali <- gg.cali + geom_abline(intercept = 0, slope = 1, color = "brown")
gg.cali <- gg.cali + geom_point(size = 2) + geom_line(size = 1) + facet_grid(beta.char~scenario) 
gg.cali <- gg.cali + labs(x = "Empirical standard deviation", y = "Average modeled standard deviation", color = "") + guides(linetype = "none")
gg.cali <- gg.cali + scale_color_manual(values = gg_color_hue(6)[1:4])
gg.cali <- gg.cali + theme_minimal() + theme(legend.position="bottom",
                                             text = element_text(size=15), 
                                             axis.line = element_line(size = 1.25),
                                             axis.ticks = element_line(size = 2),
                                             axis.ticks.length=unit(.25, "cm"),
                                             legend.key.size = unit(3,"line"))
## gg.cali

## ** power
dtS.sim[,power.upper := power + qnorm(0.975) * sqrt(power*(1-power))/sqrt(rep)]
dtS.sim[,power.lower := power + qnorm(0.025) * sqrt(power*(1-power))/sqrt(rep)]

gg.power <- ggplot(dtS.sim[!is.na(dtS.sim$power)], aes(x = n.char, y = power, linetype = target, group = type.char)) 
gg.power <- gg.power + facet_grid(beta.char~scenario) 
gg.power <- gg.power + geom_ribbon(aes(ymin = power.lower, ymax = power.upper), alpha = 0.2)
gg.power <- gg.power + geom_hline(yintercept = 0.05, color = "brown")
gg.power <- gg.power + geom_point(size = 2, aes(color = type.char)) + geom_line(size = 1, aes(color = type.char))
gg.power <- gg.power + labs(x = "Sample size (per group)", y = "Rejection rate", color = "") + guides(linetype = "none")
gg.power <- gg.power + scale_color_manual(values = gg_color_hue(6)[1:4])
gg.power <- gg.power + scale_y_continuous(labels = scales::percent)
gg.power <- gg.power + theme_minimal() + theme(legend.position="bottom",
                                               text = element_text(size=15), 
                                               axis.line = element_line(size = 1.25),
                                               axis.ticks = element_line(size = 2),
                                               axis.ticks.length=unit(.25, "cm"),
                                               legend.key.size = unit(3,"line"))
## gg.power

gg.power4.H0 <- gg.power %+% dtS.sim[type.char %in% c("proportion","pool (robust gls)") == FALSE & beta == 0,]
gg.power4.H0 <- gg.power4.H0 + coord_trans(y="log10", ylim = c(0.025,0.8))  + scale_y_continuous(breaks = c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),labels = scales::percent)
gg.power4.H0 <- gg.power4.H0 + ylab("Rejection rate (log scale)")

## * Export

## bias
png(file.path("figures","simulation-fig-bias.png"))
gg.bias
dev.off()

pdf(file.path("figures","simulation-fig-bias.pdf"), width = 10)
gg.bias
dev.off()

## standard error
png(file.path("figures","simulation-fig-dispersion.png"))
gg.se
dev.off()

pdf(file.path("figures","simulation-fig-dispersion.pdf"), width = 10)
gg.se
dev.off()

png(file.path("figures","simulation-fig-dispersion4.png"))
gg.se4
dev.off()

pdf(file.path("figures","simulation-fig-dispersion4.pdf"), width = 10)
gg.se4
dev.off()

## standard error (log scale)
png(file.path("figures","simulation-fig-dispersion-log.png"))
gg.seLog
dev.off()

pdf(file.path("figures","simulation-fig-dispersion-log.pdf"), width = 10)
gg.seLog
dev.off()

png(file.path("figures","simulation-fig-dispersion4-log.png"))
gg.seLog4
dev.off()

pdf(file.path("figures","simulation-fig-dispersion4-log.pdf"), width = 10)
gg.seLog4
dev.off()


## calibration se
png(file.path("figures","simulation-fig-calibrationSE.png"))
gg.cali
dev.off()

pdf(file.path("figures","simulation-fig-calibrationSE.pdf"), width = 10)
gg.cali
dev.off()

## power
png(file.path("figures","simulation-fig-power.png"))
gg.power
dev.off()

pdf(file.path("figures","simulation-fig-power.pdf"), width = 10)
gg.power
dev.off()

## standard error and power
png(file.path("figures","simulation-fig-sdpower4.png"))
egg::ggarrange(gg.seLog4.H0 + guides(color = "none", shape = "none")+labs(x = ""),gg.power4.H0, nrow = 2)
dev.off()

cairo_pdf(file.path("figures","simulation-fig-sdpower4.pdf"), width = 10, height = 7)
egg::ggarrange(gg.seLog4.H0 + guides(color = "none", shape = "none")+labs(x = "")+theme(plot.margin = margin(0,0,-0.5,0, "cm")),
               gg.power4.H0 + theme(plot.margin = margin(0,0,-0.5,0.25, "cm")) , nrow = 2,newpage = FALSE)
dev.off()



##----------------------------------------------------------------------
### FIGURE-simulation.R ends here
