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

temp <- lavaan(clpm10_c, data=gen_starts(n=25))



## Function to simulate clpm data
run_sim_clpm <- function(studyN,      
                       nwaves,   
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
                       xr,       
                       yr) {       
    data <- gen_starts(n=studyN,    
                       nwaves=nwaves, 
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
                       xr=xr,      
                       yr=yr       
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
                                                   nwaves=10,   # Number of waves
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
                                                   xr=.25,       # Measurement error for X
                                                   yr=.25       # Measurement error for Y
                                                   ),
                               simplify = TRUE)))



## Initialize Matrix
results <- data.frame(value = numeric(),
                      power = numeric(),
                      estimatex = numeric(),
                      estimatey = numeric())

testValues <- c(25,50,75,100,250,500,1000)

for (i in 1:length(testValues )) {
    value <- testValues[i]
    sims <- data.frame(t(mcreplicate(n=1000, run_sim_clpm(studyN=value,      # N to generate
                                                          nwaves=10,   # Number of waves
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
                                                          xr=.25,       # Measurement error for X
                                                          yr=.25       # Measurement error for Y
                                                          ), mc.cores=16)))
    results[i,1] <- value
    results[i,2] <- sum(sims$X2<.05 | sims$X4<.05)/1000
    results[i,3] <- mean(sims$X1)
    results[i,4] <- mean(sims$X3)
}


library(ggplot2)

ggplot(aes(x=value, y=power), data=results) + geom_line()
