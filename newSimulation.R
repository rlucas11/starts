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
source("scripts/starts_c.R")
source("scripts/run_sim.R")


## Function to simulate clpm data
run_sim_clpm <- function(waves=10,
                         studyN,      
                         ri_x,     
                         ri_y,     
                         cor_i,   
                         x,        
                         y,        
                         stab_x, 
                         stab_y, 
                         yx,      
                         xy,      
                         cor_xy, 
                         reliability_x,       
                         reliability_y) {
    sx <- sum(ri_x, x)/reliability_x - sum(ri_x, x)
    sy <- sum(ri_y, y)/reliability_y - sum(ri_y, y)
    data <- gen_starts(n=studyN,    
                       nwaves=waves, 
                       ri_x=ri_x,    
                       ri_y=ri_y,    
                       cor_i=cor_i,  
                       x=x,       
                       y=y,       
                       stab_x=stab_x, 
                       stab_y=stab_y, 
                       yx=yx,     
                       xy=xy,     
                       cor_xy=cor_xy, 
                       xr=sx,      
                       yr=sy       
                       )
    ##
    if (waves==10) {
        clpm <- clpm10_c
    } else if (waves==5) {
        clpm <- clpm5_c
    } else if (waves==3) {
        clpm <- clpm3_c
    } else if (waves==2) {
        clpm <- clpm2_c
    } else {
        stop("No model defined for that many waves")
    }
    fit_clpm <- lavaan(clpm, data = data)
    output <- c(parameterEstimates(fit_clpm)[3,5],
                parameterEstimates(fit_clpm)[3,8],
                parameterEstimates(fit_clpm)[4,5],
                parameterEstimates(fit_clpm)[4,8])
    output
}


## Testing
waves=2
studyN=25      # N to generate
nwaves=10   # Number of waves
ri_x=1     # Random intercept variance for X
ri_y=1     # Random intercept variance for Y
cor_i=.5   # Correlation between intercepts
x=0        # AR variance for X
y=0        # AR variance for Y
stab_x=.5  # Stability of X
stab_y=.5  # Stability of Y
yx=0      # Cross lag (Y on X)
xy=0      # Cross lag (X on Y)
cor_xy=.5  # Correlation between X and Y
cor_xyr=.25 # Correlation between X and Y residuals
xr=.25       # Measurement error for X
yr=.25       # Measurement error for Y
reliability_x <- .80
reliability_y <- .80

run_sim_clpm(waves=5,
             studyN=10000,      # N to generate
             ri_x=1,     # Random intercept variance for X
             ri_y=1,     # Random intercept variance for Y
             cor_i=0,   # Correlation between intercepts
             x=1,        # AR variance for X
             y=1,        # AR variance for Y
             stab_x=.5,  # Stability of X
             stab_y=.5,  # Stability of Y
             yx=0,      # Cross lag (Y on X)
             xy=0,      # Cross lag (X on Y)
             cor_xy=.5,  # Correlation between X and Y
             reliability_x=1,       # Measurement error for X
             reliability_y=1       # Measurement error for Y
             )


dataS <- gen_starts(n=studyN,
                    nwaves = nwaves,
                    ri_x = ri_x,
                    ri_y = ri_y,
                    cor_i = cor_i,
                    x = x,
                    y = y,
                    stab_x = stab_x,
                    stab_y = stab_y,
                    yx = yx,
                    xy = xy,
                    cor_xy = cor_xy,
                    xr = xr,
                    yr = yr
                    )
temp <- lavaan(clpm10_c, dataS)
parameterEstimates(temp)

sims <- data.frame(t(mcreplicate(n=1000, run_sim_clpm(waves = 2,
                                                    studyN=25,      # N to generate
                                                    ri_x=1,     # Random intercept variance for X
                                                    ri_y=1,     # Random intercept variance for Y
                                                    cor_i=.5,   # Correlation between intercepts
                                                    x=0,        # AR variance for X
                                                    y=0,        # AR variance for Y
                                                    stab_x=.5,  # Stability of X
                                                    stab_y=.5,  # Stability of Y
                                                    yx=0,      # Cross lag (Y on X)
                                                    xy=0,      # Cross lag (X on Y)
                                                    cor_xy=.5,  # Correlation between X and Y
                                                    reliability_x=.80,       # Measurement error for X
                                                    reliability_y=.80       # Measurement error for Y
                                                    ),
                                 mc.cores=16)))

sum(sims$X2<.05 | sims$X4<.05)/1000



run_sim_clpm(50, 1, 1, .5, 1, 1, .5, .5, 0, 0, .5, .8, .8)





################################################################################
## Simulation 10 Waves
################################################################################

## Test Simulation Values

nValues <- c(100,250)
rValues <- c(.3, .5)
reliabilities <- c(.5, .7)
## sValues1 <- 1/reliabilities-1
## sValues2 <- 1.5/reliabilities-1.5
## sValues3 <- 2/reliabilities-2
## sValues4 <- 3/reliabilities-3
arValues <- c(0,1)

## Actual Simulation Values
nValues <- c(25,50,100,250,500,1000)
rValues <- c(.1,.3,.5,.7)
reliabilities <- c(.5, .7, .9)
sValues1 <- 1/reliabilities-1
sValues2 <- 1.5/reliabilities-1.5
sValues3 <- 2/reliabilities-2
sValues4 <- 3/reliabilities-3
arValues <- c(0, .5, 1, 2)



## Run Sim
##
set.seed(119)
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
                sims <- data.frame(t(mcreplicate(n=1000, run_sim_clpm(waves=10,
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

saveRDS(results, "saved/simulation1.rds")                


## Plot

results <- readRDS("saved/simulation1.rds")
## Changes names for plot
names(results) <- c("N", "r","Reliability", "Autoregressive", "power", "estimatex","estimatey")
## Create r labels for plot
labels <- results %>%
  filter(N==25)



p <- ggplot(data = results, aes(x = N, y = power, group = r)) +
    geom_line() +
    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
    facet_grid(cols = vars(Reliability),
               rows = vars(Autoregressive),
               labeller=label_both) +
    geom_text_repel(data=labels, aes(label=r), size=3) + theme_bw() +
    ylab("Probability of One or More Significant Cross-Lagged Effects")
p

p <- ggplot(data = results, aes(x = N, y = power, group = r)) +
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


ggsave("images/10WaveSimulation.png", width=6.5, height=8.5, units="in")


### Scratch
    mcreplicate(n=1000, run_sim_clpm(studyN=nValue,      # N to generate
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
                                     ), mc.cores=14)            

p <- ggplot(data = results, aes(x = N, y = power, group = r)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) +
#    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
    facet_grid(cols = vars(Reliability),
               rows = vars(Autoregressive),
               labeller=label_both) +
    geom_text_repel(data=labels, aes(label=r), size=3) + theme_bw() +
    ylab("Probability of One or More Significant Cross-Lagged Effects")
p


        
    


################################################################################
## 3-waves
################################################################################

## Function to simulate clpm data
## This is not needed anymore. Function is general.
run_sim_clpm3 <- function(studyN,      
                       ri_x,     
                       ri_y,     
                       cor_i,   
                       x,        
                       y,        
                       stab_x, 
                       stab_y, 
                       yx,      
                       xy,      
                       cor_xy, 
                       reliability_x,       
                       reliability_y) {
    sx <- sum(ri_x, x)/reliability_x - sum(ri_x, x)
    sy <- sum(ri_y, y)/reliability_y - sum(ri_y, y)
    data <- gen_starts(n=studyN,    
                       nwaves=3, 
                       ri_x=ri_x,    
                       ri_y=ri_y,    
                       cor_i=cor_i,  
                       x=x,       
                       y=y,       
                       stab_x=stab_x, 
                       stab_y=stab_y, 
                       yx=yx,     
                       xy=xy,     
                       cor_xy=cor_xy, 
                       xr=sx,      
                       yr=sy       
                       )
    ##
    fit_clpm10 <- lavaan(clpm3_c, data = data)
    output <- c(parameterEstimates(fit_clpm10)[3,5],
                parameterEstimates(fit_clpm10)[3,8],
                parameterEstimates(fit_clpm10)[4,5],
                parameterEstimates(fit_clpm10)[4,8])
    output
}




## Actual Simulation Values
nValues <- c(25,50,100,250,500,1000)
rValues <- c(.1,.3,.5,.7)
reliabilities <- c(.5, .7, .9)
sValues1 <- 1/reliabilities-1
sValues2 <- 1.5/reliabilities-1.5
sValues3 <- 2/reliabilities-2
sValues4 <- 3/reliabilities-3
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
                      estimatey = numeric())
for (i in 1:length(nValues)) {
    for (j in 1:length(rValues)) {
        for (k in 1:length(reliabilities)) {
            for (l in 1:length(arValues)) {
                nValue <- nValues[i]
                rValue <- rValues[j]
                rel <- reliabilities[k]
                arValue <- arValues[l]
                sims <- data.frame(t(mcreplicate(n=1000, run_sim_clpm(waves = 3,
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

saveRDS(results, "saved/3WaveSimulation.rds")


results3 <- readRDS("saved/3WaveSimulation.rds")
## Changes names for plot
names(results3) <- c("N", "r","Reliability", "Autoregressive", "power", "estimatex","estimatey")
## Create r labels for plot
labels <- results3 %>%
  filter(N==25)



p <- ggplot(data = results3, aes(x = N, y = power, group = r)) +
    geom_line() +
#    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
    facet_grid(cols = vars(Reliability),
               rows = vars(Autoregressive),
               labeller=label_both) +
    geom_text_repel(data=labels, aes(label=r), size=3) + theme_bw() +
    ylab("Probability of One or More Significant Cross-Lagged Effects")
p

p <- ggplot(data = results3, aes(x = N, y = power, group = r)) +
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


ggsave("images/3WaveSimulation.png", width=6.5, height=8.5, units="in")

p <- ggplot(data = results3, aes(x = N, y = power, group = r)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) +
#    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
    facet_grid(cols = vars(Reliability),
               rows = vars(Autoregressive),
               labeller=label_both) +
    geom_text_repel(data=labels, aes(label=r), size=3) + theme_bw() +
    ylab("Probability of One or More Significant Cross-Lagged Effects")
p




###

mcreplicate(n=1000, run_sim_clpm(studyN=nValue,      # N to generate
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
                                 ), mc.cores=14)

temp <- data.frame(t(mcreplicate(n=1000, run_sim_clpm3(studyN=100,      # N to generate
                                 ri_x=1,     # Random intercept variance for X
                                 ri_y=1,     # Random intercept variance for Y
                                 cor_i=.5,   # Correlation between intercepts
                                 x=1,        # AR variance for X
                                 y=0,        # AR variance for Y
                                 stab_x=.5,  # Stability of X
                                 stab_y=.5,  # Stability of Y
                                 yx=0,      # Cross lag (Y on X)
                                 xy=0,      # Cross lag (X on Y)
                                 cor_xy=.5,  # Correlation between X and Y
                                 reliability_x=.7,       # Measurement error for X
                                 reliability_y=.7       # Measurement error for Y
                                 ), mc.cores=14)))
sum(temp$X2<.05 | temp$X4<.05)/1000


    
################################################################################
## 2-Wave Simulation
################################################################################

## Actual Simulation Values
nValues <- c(25,50,100,250,500,1000)
rValues <- c(.1,.3,.5,.7)
reliabilities <- c(.5, .7, .9)
sValues1 <- 1/reliabilities-1
sValues2 <- 1.5/reliabilities-1.5
sValues3 <- 2/reliabilities-2
sValues4 <- 3/reliabilities-3
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
labels <- results2 %>%
  filter(N==25)



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
sValues1 <- 1/reliabilities-1
sValues2 <- 1.5/reliabilities-1.5
sValues3 <- 2/reliabilities-2
sValues4 <- 3/reliabilities-3
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
## Create r labels for plot
labels <- results5 %>%
  filter(N==25)



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
## 2-Wave Simulation, no ST to get baseline error rates
################################################################################

## Actual Simulation Values
nValues <- c(1000)
rValues <- c(0)
reliabilities <- c(1)
sValues1 <- 1/reliabilities-1
sValues2 <- 1.5/reliabilities-1.5
sValues3 <- 2/reliabilities-2
sValues4 <- 3/reliabilities-3
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
                results[loopRow,5] <- sum(sims$X2<.05 | sims$X4<.05)/1000
                results[loopRow,6] <- mean(sims$X1)
                results[loopRow,7] <- mean(sims$X3)
                loopRow <- loopRow+1
            }
        }
    }
}
results

################################################################################
## 2-Wave Simulation, with Error to show effects of unreliability
################################################################################

################################################################################
## 2-Wave Simulation, no ST to get baseline error rates
################################################################################

## Actual Simulation Values
nValues <- c(500)
rValues <- c(0)
reliabilities <- c(.8)
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
                results[loopRow,5] <- sum(sims$X2<.05 | sims$X4<.05)/1000
                results[loopRow,6] <- mean(sims$X1)
                results[loopRow,7] <- mean(sims$X3)
                loopRow <- loopRow+1
            }
        }
    }
}
results


################################################################################
## 2-Wave Simulation with True CL Paths
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

## scratch

