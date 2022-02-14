## Load packages
library(lavaan) ## Run models
library(tidyverse) ## Data Manipulation
library(rethinking) ## Parallel Simulations
library(ggplot2)
library(ggrepel)


## Load Scripts
source("scripts/clpm2_c.R")
source("scripts/clpm3_c.R")
source("scripts/clpm5_c.R")
source("scripts/gen_starts.R")
source("scripts/clpm10_c.R")
source("scripts/ri_clpm10_c.R")
source("scripts/ri_clpm3_c.R")
source("scripts/starts_c.R")
source("scripts/run_sim.R")


results2 <- readRDS("saved/2WaveSimulation.rds")
## Changes names for plot
names(results2) <- c("N", "r","Reliability", "Autoregressive", "power", "estimatex","estimatey")
## Create r labels for plot



## p <- ggplot(data = results2, aes(x = N, y = power, group = r)) +
##     geom_line() +
## #    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
##     facet_grid(cols = vars(Reliability),
##                rows = vars(Autoregressive),
##                labeller=label_both) +
##     geom_text_repel(data=labels, aes(label=r), size=3) + theme_bw() +
##     ylab("Probability of One or More Significant Cross-Lagged Effects")
## p

## ggsave("images/2WaveSimulation.png", width=6.5, height=8.5, units="in")

results2 %>%
    filter(N>25, Reliability==.7, Autoregressive==1) %>%
    ggplot(aes(x = N, y = power, group = r)) +
    scale_x_continuous(breaks=c(50,100,250,500,1000)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) +
    theme_bw() +
    theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(linetype = "Correlation Between Stable Traits") +
    ylab("Probability of One or More Spurious Cross-Lagged Effects")

ggsave("images/grantFigure1.png", width=5, height=5, units="in")


results5 <- readRDS("saved/5WaveSimulation.rds")
## Changes names for plot
names(results5) <- c("N", "r","Reliability", "Autoregressive", "power", "estimatex","estimatey")
## Create r labels for plot



## p <- ggplot(data = results2, aes(x = N, y = power, group = r)) +
##     geom_line() +
## #    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
##     facet_grid(cols = vars(Reliability),
##                rows = vars(Autoregressive),
##                labeller=label_both) +
##     geom_text_repel(data=labels, aes(label=r), size=3) + theme_bw() +
##     ylab("Probability of One or More Significant Cross-Lagged Effects")
## p

## ggsave("images/2WaveSimulation.png", width=6.5, height=8.5, units="in")

results5 %>%
    filter(N>25, Reliability==.7, Autoregressive==1) %>%
    ggplot(aes(x = N, y = power, group = r)) +
    scale_x_continuous(breaks=c(50,100,250,500,1000)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) +
    theme_bw() +
    theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(linetype = "Correlation Between Stable Traits") +
    ylab("Probability of One or More Spurious Cross-Lagged Effects")

ggsave("images/grantFigure2.png", width=5, height=5, units="in")




################################################################################
## 10-Wave Simulation
################################################################################

## Actual Simulation Values
nValues <- c(25,50,100,250,500,1000)
rValues <- c(.1,.3,.5,.7)
reliabilities <- c(.7)
arValues <- c(1)



## Run Sim
##
loopRow <- 1
results <- data.frame(N = numeric(),
                      r = numeric(),
                      reliability = numeric(),
                      AR_Var = numeric(),
                      power = numeric(),
                      estimatex = numeric(),
                      estimatey = numeric())
for (i in 1:length(nValues)) {
    for (j in 1:length(rValues)) {
        for (k in 1:length(reliabilities)) {
            for (l in 1:length(arValues)) {
                nValue <- nValues[i]
                rValue <- rValues[j]
                rel <- reliabilities[k]
                arValue <- arValues[l]
                sims <- data.frame(t(mcreplicate(n=1000, run_sim_clpm(waves = 10,
                                                                      studyN=nValue,      # N to generate
                                                                      ri_x=1,     # Random intercept variance for X
                                                                      ri_y=1,     # Random intercept variance for Y
                                                                      cor_i=rValue,   # Correlation between intercepts
                                                                      x=arValue,        # AR variance for X
                                                                      y=arValue,        # AR variance for Y
                                                                      stab_x=.5,  # Stability of X
                                                                      stab_y=.5,  # Stability of Y
                                                                      yx=0,      # Cross lag (Y on X)
                                                                      xy=0,      # Cross lag (X on Y)
                                                                      cor_xy=.5,  # Correlation between X and Y
                                                                      reliability_x=rel,       # Measurement error for X
                                                                      reliability_y=rel       # Measurement error for Y
                                                                      ), mc.cores=14)))
                results[loopRow,1] <- nValue
                results[loopRow,2] <- rValue
                results[loopRow,3] <- rel
                results[loopRow,4] <- arValue
                results[loopRow,5] <- sum(sims$X2<.05 | sims$X4<.05)/1000
                results[loopRow,6] <- mean(sims$X1)
                results[loopRow,7] <- mean(sims$X3)
                loopRow <- loopRow+1
            }
        }
    }
}

saveRDS(results, "saved/10WaveSimulation.rds")

results10 <- readRDS("saved/10WaveSimulation.rds")
## Changes names for plot
names(results10) <- c("N", "r","Reliability", "Autoregressive", "power", "estimatex","estimatey")
## Create r labels for plot



## p <- ggplot(data = results2, aes(x = N, y = power, group = r)) +
##     geom_line() +
## #    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
##     facet_grid(cols = vars(Reliability),
##                rows = vars(Autoregressive),
##                labeller=label_both) +
##     geom_text_repel(data=labels, aes(label=r), size=3) + theme_bw() +
##     ylab("Probability of One or More Significant Cross-Lagged Effects")
## p

## ggsave("images/2WaveSimulation.png", width=6.5, height=8.5, units="in")

results10 %>%
    filter(N>25, Reliability==.7, Autoregressive==1) %>%
    ggplot(aes(x = N, y = power, group = r)) +
    scale_x_continuous(breaks=c(50,100,250,500,1000)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) +
    theme_bw() +
    theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(linetype = "Correlation Between Stable Traits") +
    ylab("Probability of One or More Spurious Cross-Lagged Effects")

ggsave("images/grantFigure10.png", width=5, height=5, units="in")


