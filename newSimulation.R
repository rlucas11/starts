## Load packages
library(lavaan) ## Run models
library(tidyverse) ## Data Manipulation
library(rethinking) ## Parallel Simulations
library(ggplot2)
library(ggrepel)


## Load Scripts
source("scripts/clpm3_c.R")
source("scripts/gen_starts.R")
source("scripts/clpm10_c.R")
source("scripts/ri_clpm10_c.R")
source("scripts/starts_c.R")
source("run_sim.R")

## Function to simulate clpm data
run_sim_clpm <- function(studyN,      
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
                       nwaves=10, 
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
    fit_clpm10 <- lavaan(clpm10_c, data = data)
    output <- c(parameterEstimates(fit_clpm10)[3,5],
                parameterEstimates(fit_clpm10)[3,8],
                parameterEstimates(fit_clpm10)[4,5],
                parameterEstimates(fit_clpm10)[4,8])
    output
}


## Testing
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

sims <- data.frame(t(replicate(n=1000, run_sim_clpm(studyN=25,      # N to generate
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
                               simplify = TRUE)))


run_sim_clpm(50, 1, 1, .5, 1, 1, .5, .5, 0, 0, .5, .8, .8)





################################################################################
## Simulation
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
                sims <- data.frame(t(mcreplicate(n=1000, run_sim_clpm(studyN=nValue,      # N to generate
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

ggsave("images/simulation1.png", width=6.5, height=8.5, units="in")


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





################################################################################
## 3-waves
################################################################################

## Function to simulate clpm data
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
                sims <- data.frame(t(mcreplicate(n=1000, run_sim_clpm3(studyN=nValue,      # N to generate
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

ggsave("images/3WaveSimulation.png", width=6.5, height=8.5, units="in")



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
