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


################################################################################
## 2-Wave CLMP Simulation, no Stable Trait to get baseline error rates
################################################################################

## Actual Simulation Values
nValues <- c(1000)
rValues <- c(0)
reliabilities <- c(1)
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
                sims <- data.frame(t(mcreplicate(n=10000, run_sim_clpm(waves = 2,
                                                                      studyN=nValue,      # N to generate
                                                                      ri_x=0,     # Random intercept variance for X
                                                                      ri_y=0,     # Random intercept variance for Y
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
                results[loopRow,5] <- sum(sims$X2<.05 | sims$X4<.05)/10000
                results[loopRow,6] <- mean(sims$X1)
                results[loopRow,7] <- mean(sims$X3)
                loopRow <- loopRow+1
            }
        }
    }
}
results


################################################################################
## 2-Wave Simulation
################################################################################

## Actual Simulation Values
nValues <- c(25,50,100,250,500,1000)
rValues <- c(.1,.3,.5,.7)
reliabilities <- c(.5, .7, .9)
arValues <- c(0, .5, 1, 2)



## Run Sim
##
set.seed(121)
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
                sims <- data.frame(t(mcreplicate(n=1000, run_sim_clpm(waves = 2,
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

saveRDS(results, "saved/2WaveSimulation.rds")

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

p <- ggplot(data = results2, aes(x = N, y = power, group = r)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) +
    ##    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
    scale_x_continuous(breaks = c(25,50,100,250,500,1000)) +
    facet_grid(cols = vars(Reliability),
               rows = vars(Autoregressive),
               labeller=label_both) +
    ##    geom_text_repel(data=labels, aes(label=r), size=3) +
    theme_bw() +
    theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(linetype = "Correlation Between Stable Traits") +
    ylab("Probability of One or More Significant Cross-Lagged Effects")
p

p+ geom_hline(aes(yintercept=.0975), color="gray")

ggsave("images/2WaveSimulation.png", width=6.5, height=8.5, units="in")


## scratch
temp <- results2 %>%
    group_by(r, Reliability, Autoregressive) %>%
    summarize(estimate=mean(estimatex))

temp %>%
    pivot_wider(names_from = Reliability, values_from = estimate)


################################################################################
## 5-Wave Simulation
################################################################################

## Actual Simulation Values
nValues <- c(25,50,100,250,500,1000)
rValues <- c(.1,.3,.5,.7)
reliabilities <- c(.5, .7, .9)
arValues <- c(0, .5, 1, 2)



## Run Sim
##
set.seed(122)
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
                sims <- data.frame(t(mcreplicate(n=1000, run_sim_clpm(waves = 5,
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

saveRDS(results, "saved/5WaveSimulation.rds")

results5 <- readRDS("saved/5WaveSimulation.rds")
## Changes names for plot
names(results5) <- c("N", "r","Reliability", "Autoregressive", "power", "estimatex","estimatey")



p <- ggplot(data = results5, aes(x = N, y = power, group = r)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) +
    ##    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
    scale_x_continuous(breaks = c(25,50,100,250,500,1000)) +
    facet_grid(cols = vars(Reliability),
               rows = vars(Autoregressive),
               labeller=label_both) +
    ##    geom_text_repel(data=labels, aes(label=r), size=3) +
    theme_bw() +
    theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(linetype = "Correlation Between Stable Traits") +
    ylab("Probability of One or More Significant Cross-Lagged Effects")
p

ggsave("images/5WaveSimulation.png", width=6.5, height=8.5, units="in")


################################################################################
## 2-Wave Simulation with True CL Paths
## This is original script to get estimates (but not power)
################################################################################

## Actual Simulation Values
nValues <- c(10000)
rValues <- c(.1,.3,.5,.7)
ar_corValues <- c(.1, .3, .5, .7)
reliabilities <- c(.7, 1)
arValues <- c(.5, 1, 2)

## Run Sim
##
loopRow <- 1
results <- data.frame(ar_corValue = numeric(),
                      r = numeric(),
                      reliability = numeric(),
                      AR_Var = numeric(),
                      estimatex = numeric(),
                      estimatey = numeric())
for (i in 1:length(ar_corValues)) {
    for (j in 1:length(rValues)) {
        for (k in 1:length(reliabilities)) {
            for (l in 1:length(arValues)) {
                ar_corValue <- ar_corValues[i]
#                nValue <- nValues[i]
                rValue <- rValues[j]
                rel <- reliabilities[k]
                arValue <- arValues[l]
                sims <- run_sim_clpm(waves = 2,
                                     studyN=10000,      # N to generate
                                     ri_x=1,     # Random intercept variance for X
                                     ri_y=1,     # Random intercept variance for Y
                                     cor_i=rValue,   # Correlation between intercepts
                                     x=arValue,        # AR variance for X
                                     y=arValue,        # AR variance for Y
                                     stab_x=.5,  # Stability of X
                                     stab_y=.5,  # Stability of Y
                                     yx=.5,      # Cross lag (Y on X)
                                     xy=.2,      # Cross lag (X on Y)
                                     cor_xy=ar_corValue,  # Correlation between X and Y
                                     reliability_x=rel,       # Measurement error for X
                                     reliability_y=rel       # Measurement error for Y
                                     )
                results[loopRow,1] <- ar_corValue
                results[loopRow,2] <- rValue
                results[loopRow,3] <- rel
                results[loopRow,4] <- arValue
                results[loopRow,5] <- sims[[1]]
                results[loopRow,6] <- sims[[3]]
                loopRow <- loopRow+1
            }
        }
    }
}

## Changes names for plot
names(results) <- c("rAR", "rST","Reliability", "Autoregressive", "estimate_x","estimate_y")

temp <- results %>%
    filter(Reliability==1) %>%
    select(-estimate_x, -Reliability) %>%
    pivot_wider(names_from = rST, values_from = estimate_y)

names(temp) <- c("AR r", "AR Ratio", "0.1", "0.5", "0.7", "0.9")

papaja::apa_table(temp,
                  midrules=c(3,6,9),
                  align=rep("r", 10),
                  col_spanners=list("Stable Trait Correlation"=c(3,6)),
                  caption="Average Estimated Cross-Lagged Paths In Each Simulation Condition")


temp2 <- results %>%
    filter(Reliability==.7) %>%
    select(-estimate_x, -Reliability) %>%
    pivot_wider(names_from = rST, values_from = estimate_y)

names(temp2) <- c("AR r", "AR Ratio", "0.1", "0.5", "0.7", "0.9")


################################################################################
## 2-wave CLPM simulation with real effects
## This script gets power 
################################################################################


## Actual Simulation Values
nValues <- c(50,100,250,500,1000)
rValues <- c(.1,.3,.5,.7)
reliabilities <- c(.5, .7, .9)
arValues <- c(.5, 1, 2)
clValues <- c(.1, .3, .5)

## Run Sim
##
set.seed(121)
loopRow <- 1
results <- data.frame(N = numeric(),
                      r = numeric(),
                      reliability = numeric(),
                      AR_Var = numeric(),
                      cl_Var = numeric(),
                      powerx = numeric(),
                      powery = numeric(),
                      estimatex = numeric(),
                      estimatey = numeric())
for (i in 1:length(nValues)) {
    for (j in 1:length(rValues)) {
        for (k in 1:length(reliabilities)) {
            for (l in 1:length(arValues)) {
                for (m in 1:length(clValues)) {
                    nValue <- nValues[i]
                    rValue <- rValues[j]
                    rel <- reliabilities[k]
                    arValue <- arValues[l]
                    clValue <- clValues[m]
                    sims <- data.frame(t(mcreplicate(n=1000, run_sim_clpm(waves = 2,
                                                                          studyN=nValue,      # N to generate
                                                                          ri_x=1,     # Random intercept variance for X
                                                                          ri_y=1,     # Random intercept variance for Y
                                                                          cor_i=rValue,   # Correlation between intercepts
                                                                          x=arValue,        # AR variance for X
                                                                          y=arValue,        # AR variance for Y
                                                                          stab_x=.5,  # Stability of X
                                                                          stab_y=.5,  # Stability of Y
                                                                          yx=clValue,      # Cross lag (Y on X)
                                                                          xy=0,      # Cross lag (X on Y)
                                                                          cor_xy=.5,  # Correlation between X and Y
                                                                          reliability_x=rel,       # Measurement error for X
                                                                          reliability_y=rel       # Measurement error for Y
                                                                          ), mc.cores=14)))
                    results[loopRow,1] <- nValue
                    results[loopRow,2] <- rValue
                    results[loopRow,3] <- rel
                    results[loopRow,4] <- arValue
                    results[loopRow,5] <- clValue
                    results[loopRow,6] <- sum(sims$X2<.05)/1000
                    results[loopRow,7] <- sum(sims$X4<.05)/1000
                    results[loopRow,8] <- mean(sims$X1)
                    results[loopRow,9] <- mean(sims$X3)
                    loopRow <- loopRow+1
                }
            }
        }
    }
}

saveRDS(results, "saved/2WaveSimulationPower.rds")

resultsCl <- readRDS("saved/2WaveSimulationPower.rds")
## Changes names for plot
names(resultsCl) <- c("N", "r","Reliability", "Autoregressive", "clValue", "powerx","powery", "estimatex","estimatey")
## Create r labels for plot


## Need different plots because of different CL Values
resultsCl %>%
    filter(clValue==.3) %>%
    ggplot(aes(x = N, y = powery, group = r)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) +
    ##    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
    scale_x_continuous(breaks = c(50,100,250,500,1000)) +
    facet_grid(cols = vars(Reliability),
               rows = vars(Autoregressive),
               labeller=label_both) +
    ##    geom_text_repel(data=labels, aes(label=r), size=3) +
    theme_bw() +
    theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(linetype = "Correlation Between Stable Traits") +
    ylab("Probability of One or More Significant Cross-Lagged Effects")


resultsCl %>%
    filter(clValue==.1) %>%
    ggplot(aes(x = N, y = powery, group = r)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) +
    ##    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
    scale_x_continuous(breaks = c(50,100,250,500,1000)) +
    facet_grid(cols = vars(Reliability),
               rows = vars(Autoregressive),
               labeller=label_both) +
    ##    geom_text_repel(data=labels, aes(label=r), size=3) +
    theme_bw() +
    theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(linetype = "Correlation Between Stable Traits") +
    ylab("Probability of One or More Significant Cross-Lagged Effects")


################################################################################
## 3-wave CLPM simulation with real effects
## This script gets power 
################################################################################


## Actual Simulation Values
nValues <- c(50,100,250,500,1000)
rValues <- c(.1,.3,.5,.7)
reliabilities <- c(.5, .7, .9)
arValues <- c(.5, 1, 2)
clValues <- c(.1, .3, .5)

## Run Sim
##
#set.seed(121)
loopRow <- 1
results <- data.frame(N = numeric(),
                      r = numeric(),
                      reliability = numeric(),
                      AR_Var = numeric(),
                      cl_Var = numeric(),
                      powerx = numeric(),
                      powery = numeric(),
                      estimatex = numeric(),
                      estimatey = numeric())
for (i in 1:length(nValues)) {
    for (j in 1:length(rValues)) {
        for (k in 1:length(reliabilities)) {
            for (l in 1:length(arValues)) {
                for (m in 1:length(clValues)) {
                    nValue <- nValues[i]
                    rValue <- rValues[j]
                    rel <- reliabilities[k]
                    arValue <- arValues[l]
                    clValue <- clValues[m]
                    sims <- data.frame(t(mcreplicate(n=1000, run_sim_clpm(waves = 3,
                                                                          studyN=nValue,      # N to generate
                                                                          ri_x=1,     # Random intercept variance for X
                                                                          ri_y=1,     # Random intercept variance for Y
                                                                          cor_i=rValue,   # Correlation between intercepts
                                                                          x=arValue,        # AR variance for X
                                                                          y=arValue,        # AR variance for Y
                                                                          stab_x=.5,  # Stability of X
                                                                          stab_y=.5,  # Stability of Y
                                                                          yx=clValue,      # Cross lag (Y on X)
                                                                          xy=0,      # Cross lag (X on Y)
                                                                          cor_xy=.5,  # Correlation between X and Y
                                                                          reliability_x=rel,       # Measurement error for X
                                                                          reliability_y=rel       # Measurement error for Y
                                                                          ), mc.cores=14)))
                    results[loopRow,1] <- nValue
                    results[loopRow,2] <- rValue
                    results[loopRow,3] <- rel
                    results[loopRow,4] <- arValue
                    results[loopRow,5] <- clValue
                    results[loopRow,6] <- sum(sims$X2<.05)/1000
                    results[loopRow,7] <- sum(sims$X4<.05)/1000
                    results[loopRow,8] <- mean(sims$X1)
                    results[loopRow,9] <- mean(sims$X3)
                    loopRow <- loopRow+1
                }
            }
        }
    }
}

saveRDS(results, "saved/3WaveCLPMSimulationPower.rds")

resultsCl <- readRDS("saved/3WaveCLPMSimulationPower.rds")
## Changes names for plot
names(resultsCl) <- c("N", "r","Reliability", "Autoregressive", "clValue", "powerx","powery", "estimatex","estimatey")
## Create r labels for plot


## Need different plots because of different CL Values
resultsCl %>%
    filter(clValue==.3) %>%
    ggplot(aes(x = N, y = powery, group = r)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) +
    ##    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
    scale_x_continuous(breaks = c(50,100,250,500,1000)) +
    facet_grid(cols = vars(Reliability),
               rows = vars(Autoregressive),
               labeller=label_both) +
    ##    geom_text_repel(data=labels, aes(label=r), size=3) +
    theme_bw() +
    theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(linetype = "Correlation Between Stable Traits") +
    ylab("Probability of One or More Significant Cross-Lagged Effects")


resultsCl %>%
    filter(clValue==.1) %>%
    ggplot(aes(x = N, y = powery, group = r)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) +
    ##    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
    scale_x_continuous(breaks = c(50,100,250,500,1000)) +
    facet_grid(cols = vars(Reliability),
               rows = vars(Autoregressive),
               labeller=label_both) +
    ##    geom_text_repel(data=labels, aes(label=r), size=3) +
    theme_bw() +
    theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(linetype = "Correlation Between Stable Traits") +
    ylab("Probability of One or More Significant Cross-Lagged Effects")

ggsave("images/3WaveSimulation.clpm.power.cl1.png", width=6.5, height=8.5, units="in")

################################################################################
## 3-Wave RI-CLPM Simulation (with measurement error)
################################################################################

## Actual Simulation Values
nValues <- c(50,100,250,500,1000)
rValues <- c(.1,.3,.5,.7)
reliabilities <- c(.5, .7, .9)
arValues <- c(0, .5, 1, 2)

## Run Sim
##
set.seed(120)
loopRow <- 1
results <- data.frame(N = numeric(),
                      r = numeric(),
                      reliability = numeric(),
                      AR_Var = numeric(),
                      power = numeric(),
                      estimatex = numeric(),
                      estimatey = numeric(),
                      problems = numeric())
for (i in 1:length(nValues)) {
    for (j in 1:length(rValues)) {
        for (k in 1:length(reliabilities)) {
            for (l in 1:length(arValues)) {
                nValue <- nValues[i]
                rValue <- rValues[j]
                rel <- reliabilities[k]
                arValue <- arValues[l]
                sims <- data.frame(t(mcreplicate(n=1000, run_sim_riclpm(waves = 3,
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
                results[loopRow,5] <- sum(sims$pvalue<.05 | sims$pvalue.1<.05, na.rm = TRUE)/(1000-sum(is.na(sims$pvalue)))
                results[loopRow,6] <- mean(sims$est)
                results[loopRow,7] <- mean(sims$est.1)
                results[loopRow,8] <- sum(is.na(sims$pvalue))
                loopRow <- loopRow+1
            }
        }
    }
}


saveRDS(results, "saved/riSimulation3.rds")                


results <- readRDS("saved/riSimulation3.rds")
## Changes names for plot
names(results) <- c("N", "r","Reliability", "Autoregressive", "power", "estimatex","estimatey", "problems")

ggplot(data = results, aes(x = N, y = power, group = r)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) +
    scale_x_continuous(breaks = c(50,100,250,500,1000)) +
    facet_grid(cols = vars(Reliability),
               rows = vars(Autoregressive),
               labeller=label_both) +
    ##    geom_text_repel(data=labels, aes(label=r), size=3) +
    theme_bw() +
    theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(linetype = "Correlation Between Stable Traits") +
    ylab("Probability of One or More Significant Cross-Lagged Effects")

ggsave("images/3WaveSimulationRICLPM.png", width=6.5, height=8.5, units="in")

results %>%
    filter(Reliability==.5, N==1000, Autoregressive==0)


################################################################################
## 3-wave RI-CLPM with true cl paths
## This script varies correlation between stable traits
################################################################################

## Actual Simulation Values
nValues <- c(50,100,250,500,1000)
rValues <- c(.1,.3,.5,.7)
reliabilities <- c(.5, .7, .9)
arValues <- c(.5, 1, 2)
clValues <- c(.1, .2, .3)

## Run Sim
##
#set.seed(121)
loopRow <- 1
results <- data.frame(N = numeric(),
                      r = numeric(),
                      reliability = numeric(),
                      AR_Var = numeric(),
                      cl_Var = numeric(),
                      powerx = numeric(),
                      powery = numeric(),
                      estimatex = numeric(),
                      estimatey = numeric(),
                      problems = numeric())
for (i in 1:length(nValues)) {
    for (j in 1:length(rValues)) {
        for (k in 1:length(reliabilities)) {
            for (l in 1:length(arValues)) {
                for (m in 1:length(clValues)) {
                    nValue <- nValues[i]
                    rValue <- rValues[j]
                    rel <- reliabilities[k]
                    arValue <- arValues[l]
                    clValue <- clValues[m]
                    sims <- data.frame(t(mcreplicate(n=1000, run_sim_riclpm(waves = 3,
                                                                          studyN=nValue,      # N to generate
                                                                          ri_x=1,     # Random intercept variance for X
                                                                          ri_y=1,     # Random intercept variance for Y
                                                                          cor_i=rValue,   # Correlation between intercepts
                                                                          x=arValue,        # AR variance for X
                                                                          y=arValue,        # AR variance for Y
                                                                          stab_x=.5,  # Stability of X
                                                                          stab_y=.5,  # Stability of Y
                                                                          yx=clValue,      # Cross lag (Y on X)
                                                                          xy=0,      # Cross lag (X on Y)
                                                                          cor_xy=.5,  # Correlation between X and Y
                                                                          reliability_x=rel,       # Measurement error for X
                                                                          reliability_y=rel       # Measurement error for Y
                                                                          ), mc.cores=14)))
                    results[loopRow,1] <- nValue
                    results[loopRow,2] <- rValue
                    results[loopRow,3] <- rel
                    results[loopRow,4] <- arValue
                    results[loopRow,5] <- clValue
                    results[loopRow,6] <- sum(sims$pvalue<.05, na.rm = TRUE)/(1000-sum(is.na(sims$pvalue)))
                    results[loopRow,7] <- sum(sims$pvalue.1<.05, na.rm = TRUE)/(1000-sum(is.na(sims$pvalue.1)))
                    results[loopRow,8] <- mean(sims$est, na.rm=TRUE)
                    results[loopRow,9] <- mean(sims$est.1, na.rm=TRUE)
                    results[loopRow,10] <- sum(is.na(sims$pvalue))
                    loopRow <- loopRow+1
                }
            }
        }
    }
}


saveRDS(results, "saved/3WaveSimulationPowerRI.rds")

resultsClRi <- readRDS("saved/3WaveSimulationPowerRI.rds")
## Changes names for plot
names(resultsClRi) <- c("N", "r","Reliability", "Autoregressive", "clValue", "powerx","powery", "estimatex","estimatey", "problems")
## Create r labels for plot

resultsClRi %>%
    filter(clValue==.1) %>%
    ggplot(aes(x = N, y = powerx, group = r)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) +
    ##    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
    scale_x_continuous(breaks = c(50,100,250,500,1000)) +
    facet_grid(cols = vars(Reliability),
               rows = vars(Autoregressive),
               labeller=label_both) +
    ##    geom_text_repel(data=labels, aes(label=r), size=3) +
    theme_bw() +
    theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(linetype = "Correlation Between Stable Traits") +
    ylab("Probability of One or More Significant Cross-Lagged Effects")

resultsClRi %>%
    filter(clValue==.1, Autoregressive==1, Reliability==.7, r==.3)

resultsCl %>%
    filter(clValue==.1, Autoregressive==1, Reliability==.7, r==.3)


################################################################################
## 3-Wave RI-CLPM Simulation (with measurement error)
## Test effect of varying r between AR on false positives
################################################################################

## Actual Simulation Values
nValues <- c(50,100,250,500,1000)
rValues <- c(0,.1,.3,.5,.7)
reliabilities <- c(.5,.7,.9)
arValues <- c(0,.5,1,2)

## Run Sim
##
#set.seed(120)
loopRow <- 1
results <- data.frame(N = numeric(),
                      r = numeric(),
                      reliability = numeric(),
                      AR_Var = numeric(),
                      power = numeric(),
                      estimatex = numeric(),
                      estimatey = numeric(),
                      problems = numeric())
for (i in 1:length(nValues)) {
    for (j in 1:length(rValues)) {
        for (k in 1:length(reliabilities)) {
            for (l in 1:length(arValues)) {
                nValue <- nValues[i]
                rValue <- rValues[j]
                rel <- reliabilities[k]
                arValue <- arValues[l]
                sims <- data.frame(t(mcreplicate(n=1000, run_sim_riclpm(waves = 3,
                                                                        studyN=nValue,      # N to generate
                                                                        ri_x=1,     # Random intercept variance for X
                                                                        ri_y=1,     # Random intercept variance for Y
                                                                        cor_i=.3,   # Correlation between intercepts
                                                                        x=arValue,        # AR variance for X
                                                                        y=arValue,        # AR variance for Y
                                                                        stab_x=.5,  # Stability of X
                                                                        stab_y=.5,  # Stability of Y
                                                                        yx=0,      # Cross lag (Y on X)
                                                                        xy=0,      # Cross lag (X on Y)
                                                                        cor_xy=rValue,  # Correlation between X and Y
                                                                        reliability_x=rel,       # Measurement error for X
                                                                        reliability_y=rel       # Measurement error for Y
                                                                        ), mc.cores=14)))
                results[loopRow,1] <- nValue
                results[loopRow,2] <- rValue
                results[loopRow,3] <- rel
                results[loopRow,4] <- arValue
                results[loopRow,5] <- sum(sims$pvalue<.05 | sims$pvalue.1<.05, na.rm = TRUE)/(1000-sum(is.na(sims$pvalue)))
                results[loopRow,6] <- mean(sims$est)
                results[loopRow,7] <- mean(sims$est.1)
                results[loopRow,8] <- sum(is.na(sims$pvalue))
                loopRow <- loopRow+1
            }
        }
    }
}


saveRDS(results, "saved/riSimulation3_rAr.rds")                


results <- readRDS("saved/riSimulation3_rAr.rds")
## Changes names for plot
names(results) <- c("N", "r","Reliability", "Autoregressive", "power", "estimatex","estimatey", "problems")

ggplot(data = results, aes(x = N, y = power, group = r)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) +
    scale_x_continuous(breaks = c(50,100,250,500,1000)) + expand_limits(y=c(0,1)) +
    facet_grid(cols = vars(Reliability),
               rows = vars(Autoregressive),
               labeller=label_both) +
    ## ##    geom_text_repel(data=labels, aes(label=r), size=3) +
    theme_bw() +
    theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(linetype = "Correlation Between AR Traits") +
    ylab("Probability of One or More Significant Cross-Lagged Effects")

ggsave("images/3WaveSimulationRI.arCor.png", width=6.5, height=8.5, units="in")

results %>%
    filter(Reliability==.5, N==1000, Autoregressive==0)


################################################################################
## 3-wave RI-CLPM with true cl paths
################################################################################

## Actual Simulation Values
nValues <- c(50,100,250,500,1000)
rValues <- c(.1,.3,.5,.7)
reliabilities <- c(.5, .7, .9)
arValues <- c(.5, 1, 2)
clValues <- c(.1, .2, .3)

## Run Sim
##
#set.seed(121)
loopRow <- 1
results <- data.frame(N = numeric(),
                      r = numeric(),
                      reliability = numeric(),
                      AR_Var = numeric(),
                      cl_Var = numeric(),
                      powerx = numeric(),
                      powery = numeric(),
                      estimatex = numeric(),
                      estimatey = numeric(),
                      problems = numeric())
for (i in 1:length(nValues)) {
    for (j in 1:length(rValues)) {
        for (k in 1:length(reliabilities)) {
            for (l in 1:length(arValues)) {
                for (m in 1:length(clValues)) {
                    nValue <- nValues[i]
                    rValue <- rValues[j]
                    rel <- reliabilities[k]
                    arValue <- arValues[l]
                    clValue <- clValues[m]
                    sims <- data.frame(t(mcreplicate(n=1000, run_sim_riclpm(waves = 3,
                                                                          studyN=nValue,      # N to generate
                                                                          ri_x=1,     # Random intercept variance for X
                                                                          ri_y=1,     # Random intercept variance for Y
                                                                          cor_i=rValue,   # Correlation between intercepts
                                                                          x=arValue,        # AR variance for X
                                                                          y=arValue,        # AR variance for Y
                                                                          stab_x=.5,  # Stability of X
                                                                          stab_y=.5,  # Stability of Y
                                                                          yx=clValue,      # Cross lag (Y on X)
                                                                          xy=0,      # Cross lag (X on Y)
                                                                          cor_xy=.5,  # Correlation between X and Y
                                                                          reliability_x=rel,       # Measurement error for X
                                                                          reliability_y=rel       # Measurement error for Y
                                                                          ), mc.cores=14)))
                    results[loopRow,1] <- nValue
                    results[loopRow,2] <- rValue
                    results[loopRow,3] <- rel
                    results[loopRow,4] <- arValue
                    results[loopRow,5] <- clValue
                    results[loopRow,6] <- sum(sims$pvalue<.05, na.rm = TRUE)/(1000-sum(is.na(sims$pvalue)))
                    results[loopRow,7] <- sum(sims$pvalue.1<.05, na.rm = TRUE)/(1000-sum(is.na(sims$pvalue.1)))
                    results[loopRow,8] <- mean(sims$est, na.rm=TRUE)
                    results[loopRow,9] <- mean(sims$est.1, na.rm=TRUE)
                    results[loopRow,10] <- sum(is.na(sims$pvalue))
                    loopRow <- loopRow+1
                }
            }
        }
    }
}


saveRDS(results, "saved/3WaveSimulationPowerRI.rds")

resultsClRi <- readRDS("saved/3WaveSimulationPowerRI.rds")
## Changes names for plot
names(resultsClRi) <- c("N", "r","Reliability", "Autoregressive", "clValue", "powerx","powery", "estimatex","estimatey", "problems")
## Create r labels for plot

resultsClRi %>%
    filter(clValue==.1) %>%
    ggplot(aes(x = N, y = powerx, group = r)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) + expand_limits(y=c(0,1))+
    ##    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
    scale_x_continuous(breaks = c(50,100,250,500,1000)) +
    facet_grid(cols = vars(Reliability),
               rows = vars(Autoregressive),
               labeller=label_both) +
    ##    geom_text_repel(data=labels, aes(label=r), size=3) +
    theme_bw() +
    theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(linetype = "Correlation Between Stable Traits") +
    ylab("Probability of One or More Significant Cross-Lagged Effects")

ggsave("images/3WaveSimulation.ri.power.cl1.png", width=6.5, height=8.5, units="in")

resultsClRi %>%
    filter(clValue==.1, Autoregressive==1, Reliability==.7, r==.1)

resultsCl %>%
    filter(clValue==.1, Autoregressive==1, Reliability==.7, r==.1)

################################################################################
## 3-wave RI-CLPM with true cl paths
## Power
## Vary AR Cor
################################################################################

## Actual Simulation Values
nValues <- c(50,100,250,500,1000)
rValues <- c(.1,.3,.5,.7)
reliabilities <- c(.5, .7, .9)
arValues <- c(.5, 1, 2)
clValues <- c(.1, .2, .3)

## Run Sim
##
#set.seed(121)
loopRow <- 1
results <- data.frame(N = numeric(),
                      r = numeric(),
                      reliability = numeric(),
                      AR_Var = numeric(),
                      cl_Var = numeric(),
                      powerx = numeric(),
                      powery = numeric(),
                      estimatex = numeric(),
                      estimatey = numeric(),
                      problems = numeric())
for (i in 1:length(nValues)) {
    for (j in 1:length(rValues)) {
        for (k in 1:length(reliabilities)) {
            for (l in 1:length(arValues)) {
                for (m in 1:length(clValues)) {
                    nValue <- nValues[i]
                    rValue <- rValues[j]
                    rel <- reliabilities[k]
                    arValue <- arValues[l]
                    clValue <- clValues[m]
                    sims <- data.frame(t(mcreplicate(n=1000, run_sim_riclpm(waves = 3,
                                                                          studyN=nValue,      # N to generate
                                                                          ri_x=1,     # Random intercept variance for X
                                                                          ri_y=1,     # Random intercept variance for Y
                                                                          cor_i=.3,   # Correlation between intercepts
                                                                          x=arValue,        # AR variance for X
                                                                          y=arValue,        # AR variance for Y
                                                                          stab_x=.5,  # Stability of X
                                                                          stab_y=.5,  # Stability of Y
                                                                          yx=clValue,      # Cross lag (Y on X)
                                                                          xy=0,      # Cross lag (X on Y)
                                                                          cor_xy=rValue,  # Correlation between X and Y
                                                                          reliability_x=rel,       # Measurement error for X
                                                                          reliability_y=rel       # Measurement error for Y
                                                                          ), mc.cores=14)))
                    results[loopRow,1] <- nValue
                    results[loopRow,2] <- rValue
                    results[loopRow,3] <- rel
                    results[loopRow,4] <- arValue
                    results[loopRow,5] <- clValue
                    results[loopRow,6] <- sum(sims$pvalue<.05, na.rm = TRUE)/(1000-sum(is.na(sims$pvalue)))
                    results[loopRow,7] <- sum(sims$pvalue.1<.05, na.rm = TRUE)/(1000-sum(is.na(sims$pvalue.1)))
                    results[loopRow,8] <- mean(sims$est, na.rm=TRUE)
                    results[loopRow,9] <- mean(sims$est.1, na.rm=TRUE)
                    results[loopRow,10] <- sum(is.na(sims$pvalue))
                    loopRow <- loopRow+1
                }
            }
        }
    }
}


saveRDS(results, "saved/3WaveSimulationPowerRI.varyAR.rds")

resultsRiVaryAR <- readRDS("saved/3WaveSimulationPowerRI.varyAR.rds")
## Changes names for plot
names(resultsRiVaryAR) <- c("N", "r","Reliability", "Autoregressive", "clValue", "powerx","powery", "estimatex","estimatey", "problems")


resultsRiVaryAR %>%
    filter(clValue==.3) %>%
    ggplot(aes(x = N, y = powerx, group = r)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) + expand_limits(y=c(0,1))+
    ##    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
    scale_x_continuous(breaks = c(50,100,250,500,1000)) +
    facet_grid(cols = vars(Reliability),
               rows = vars(Autoregressive),
               labeller=label_both) +
    ##    geom_text_repel(data=labels, aes(label=r), size=3) +
    theme_bw() +
    theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(linetype = "Correlation Between Stable Traits") +
    ylab("Probability of One or More Significant Cross-Lagged Effects")

ggsave("images/3WaveSimulation.ri.power.cl1.png", width=6.5, height=8.5, units="in")

resultsRiVaryAR %>%
    filter(clValue==.2, Autoregressive==1, Reliability==.7, r==.5)

resultsCl %>%
    filter(clValue==.1, Autoregressive==1, Reliability==.7, r==.1)
