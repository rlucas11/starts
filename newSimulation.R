## Load packages
library(lavaan) ## Run models
library(tidyverse) ## Data Manipulation
library(rethinking) ## Parallel Simulations

## Load Scripts
source("scripts/gen_starts.R")
source("scripts/clpm10_c.R")
source("scripts/ri_clpm10_c.R")
source("scripts/starts_c.R")
source("run_sim.R")

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



## Initialize Matrix

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
loopRow <- 1
results <- data.frame(value1 = numeric(),
                      value2 = numeric(),
                      value3 = numeric(),
                      value4 = numeric(),
                      value5 = numeric(),
                      power = numeric(),
                      estimatex = numeric(),
                      estimatey = numeric())
for (i in 1:length(nValues)) {
    for (j in 1:length(rValues)) {
        for (k in 1:length(reliabilities)) {
            for (l in 1:length(arValues) {
                nValue <- values1[i]
                rValue <- values2[j]
                rel <- reliabilities[k]
                arValue <- arValues[l]
                sims <- data.frame(t(mcreplicate(n=1000, run_sim_clpm(studyN=nValue,      # N to generate
                                                                      nwaves=10,   # Number of waves
                                                                      ri_x=1,     # Random intercept variance for X
                                                                      ri_y=1,     # Random intercept variance for Y
                                                                      cor_i=rValue,   # Correlation between intercepts
                                                                      x=1,        # AR variance for X
                                                                      y=1,        # AR variance for Y
                                                                      stab_x=.5,  # Stability of X
                                                                      stab_y=.5,  # Stability of Y
                                                                      yx=0,      # Cross lag (Y on X)
                                                                      xy=0,      # Cross lag (X on Y)
                                                                      cor_xy=.5,  # Correlation between X and Y
                                                                      xr=.5,       # Measurement error for X
                                                                      yr=.5       # Measurement error for Y
                                                                      ), mc.cores=14)))
                loopRow <- loopRow+1
                results[loopRow,1] <- value1
                results[loopRow,2] <- value2
                results[loopRow,3] <- sum(sims$X2<.05 | sims$X4<.05)/1000
                results[loopRow,4] <- mean(sims$X1)
                results[loopRow,5] <- mean(sims$X3)
            }
            }



library(ggplot2)

ggplot(aes(x=value1, y=power, group=value2), data=results) + geom_line()

+ scale_x_log10()
