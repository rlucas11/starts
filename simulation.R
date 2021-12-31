################################################################################
## Scratch pad for clpm/starts project
################################################################################



library(lavaan)
library(tidyverse)

source("gen_starts.R")
source("clpm10.R")
#source("clpm5.R")
source("clpm3.R")
source("clpm2.R")
source("ri-clpm10.R")
source("starts.R")

source("clpm10_c.R")
source("clpm3_c.R")
source("ri_clpm10_c.R")
source("starts_c.R")

source("run_sim.R")


data <- gen_starts(nwaves=10, n=500, yx=0, xy=0, xr=.2, yr=.2, ri_x=1, ri_y=1, stab_x = .5, stab_y = .50, cor_i = .5)
fit_clpm2 <- lavaan(clpm2, data = data)
summary(fit_clpm2)

fit_clpm3 <- lavaan(clpm3, data = data)
fit_clpm5 <- lavaan(clpm5, data = data)
fit_clpm10 <- lavaan(clpm10, data = data)
fit_riclpm <- lavaa(

summary(fit_clpm3)
summary(fit_clpm5)
summary(fit_clpm10)


data <- gen_starts(nwaves=10, n=500, yx=0, xy=0, xr=.2, yr=.2, ri_x=.5, ri_y=.5, stab_x = 1, stab_y = .50, cor_i = .5)

data <- gen_starts(nwaves = 10, n = 500)
fit_clpm2 <- lavaan(clpm2, data = data)

summary(fit_clpm2)


fit_ri_clpm <- lavaan(m_riclpm, data = data)
summary(fit_ri_clpm)



data <- gen_starts(n=10000,      # N to generate
                    nwaves=10,   # Number of waves
                    ri_x=1,     # Random intercept variance for X
                    ri_y=1,     # Random intercept variance for Y
                    cor_i=.5,   # Correlation between intercepts
                    x=1,        # AR variance for X
                    y=1,        # AR variance for Y
                    stab_x=.5,  # Stability of X
                    stab_y=.5,  # Stability of Y
                    yx=0,      # Cross lag (Y on X)
                    xy=0,      # Cross lag (X on Y)
                    cor_xy=.5,  # Correlation between X and Y
                    cor_xyr=.1, # Correlation between X and Y residuals
                    xr=1,       # Measurement error for X
                    yr=1        # Measurement error for Y
                   )

fit_clpm2 <- lavaan(clpm2, data = data)
summary(fit_clpm2)
standardizedSolution(fit_clpm2)

fit_clpm3 <- lavaan(clpm3, data = data)
summary(fit_clpm3)
standardizedSolution(fit_clpm3)


fit_clpm10 <- lavaan(clpm10, data = data)
summary(fit_clpm10)

fit_clpm10c <- lavaan(clpm10_c, data = data)
summary(fit_clpm10c)

fit_ri <- lavaan(ri_clpm10, data = data)
summary(fit_ri)

fit_ri_c <- lavaan(ri_clpm10_c, data = data)
summary(fit_ri_c)

fit_starts <- lavaan(starts, data = data)
summary(fit_starts)

fit_starts_c <- lavaan(starts_c, data = data)
summary(fit_starts_c)

################################################################################
## Run Simulation
################################################################################

run_sim <- function() {
    data <- gen_starts(n=n,      # N to generate
                       nwaves=nwaves,   # Number of waves
                       ri_x=ri_x,     # Random intercept variance for X
                       ri_y=ri_y,     # Random intercept variance for Y
                       cor_i=cor_i,   # Correlation between intercepts
                       x=x,        # AR variance for X
                       y=y,        # AR variance for Y
                       stab_x=stab_x,  # Stability of X
                       stab_y=stab_y,  # Stability of Y
                       yx=yx,      # Cross lag (Y on X)
                       xy=xy,      # Cross lag (X on Y)
                       cor_xy=cor_xy,  # Correlation between X and Y
                       cor_xyr=cor_xyr, # Correlation between X and Y residuals
                       xr=xr,       # Measurement error for X
                       yr=yr        # Measurement error for Y
                       )
    ##
    fit_clpm2 <- lavaan(clpm2, data = data)
    fit_clpm5 <- lavaan(clpm5, data = data)
    fit_clpm10 <- lavaan(clpm10, data = data)
    fit_riclpm <- lavaan(ri_clpm10, data = data)
    output <- list(clpm2=parameterEstimates(fit_clpm2),
                   clpm5=parameterEstimates(fit_clpm5),
                   clpm10=parameterEstimates(fit_clpm10),
                   riclpm=parameterEstimates(fit_riclpm))
    output
}


run_sim <- function() {
    data <- gen_starts(n=n,      # N to generate
                       nwaves=nwaves,   # Number of waves
                       ri_x=ri_x,     # Random intercept variance for X
                       ri_y=ri_y,     # Random intercept variance for Y
                       cor_i=cor_i,   # Correlation between intercepts
                       x=x,        # AR variance for X
                       y=y,        # AR variance for Y
                       stab_x=stab_x,  # Stability of X
                       stab_y=stab_y,  # Stability of Y
                       yx=yx,      # Cross lag (Y on X)
                       xy=xy,      # Cross lag (X on Y)
                       cor_xy=cor_xy,  # Correlation between X and Y
                       cor_xyr=cor_xyr, # Correlation between X and Y residuals
                       xr=xr,       # Measurement error for X
                       yr=yr        # Measurement error for Y
                       )
    ##
    #fit_clpm2 <- lavaan(clpm2, data = data)
    #fit_clpm5 <- lavaan(clpm5, data = data)
    fit_clpm10 <- lavaan(clpm10, data = data)
    fit_riclpm <- lavaan(ri_clpm10, data = data)
    output <- list(#clpm2=parameterEstimates(fit_clpm2),
                   #clpm5=parameterEstimates(fit_clpm5),
                   clpm10=parameterEstimates(fit_clpm10),
                   riclpm=parameterEstimates(fit_riclpm))
    output
}


n=500      # N to generate
nwaves=10   # Number of waves
ri_x=1     # Random intercept variance for X
ri_y=1     # Random intercept variance for Y
cor_i=.5   # Correlation between intercepts
x=1        # AR variance for X
y=1        # AR variance for Y
stab_x=.5  # Stability of X
stab_y=.5  # Stability of Y
yx=.4      # Cross lag (Y on X)
xy=.1      # Cross lag (X on Y)
cor_xy=.5  # Correlation between X and Y
cor_xyr=.25 # Correlation between X and Y residuals
xr=.0       # Measurement error for X
yr=.0       # Measurement error for Y



sims <- replicate(n=100, run_sim(), simplify = FALSE)

sims2 <- replicate(n=100, run_sim2(), simplify = FALSE)
sims5 <- replicate(n=100, run_sim5(), simplify = FALSE)
sims10 <- replicate(n=100, run_sim10(), simplify = FALSE)

simsr <- replicate(n=100, run_simri(), simplify = FALSE)

#sims <- mcreplicate(n=100, run_sim(), mc.cores = 16)



sims %>%
    map("clpm2") %>%
    map(filter, lhs=="y2", rhs=="x1") %>%
    map("est") %>%
    unlist %>%
    mean
sims %>%
    map("clpm2") %>%
    map(filter, lhs=="x2", rhs=="y1") %>%
    map("est") %>%
    unlist %>%
    mean
sims %>%
    map("clpm10") %>%
    map(filter, lhs=="y10", rhs=="x9") %>%
    map("est") %>%
    unlist %>%
    mean
sims %>%
    map("clpm10") %>%
    map(filter, lhs=="x10", rhs=="y9") %>%
    map("est") %>%
    unlist %>%
    mean


sims %>%
    map("clpm2") %>%
    map(filter, lhs=="y2", rhs=="x1") %>%
    map("est") %>%
    unlist %>%
    hist()



sims %>%
    map("clpm10") %>%
    map(filter, lhs=="y2", rhs=="x1") %>%
    map("est") %>%
    unlist %>%
    hist()


temp <- sims %>%
    map("clpm10") %>%
    map("est") %>%
    abind(along = 2)

temp <- sims %>%
map("clpm10") %>%
map(filter, lhs=="y2", rhs=="x1") %>%


    
################################################################################
## scratch
################################################################################


sims %>%


lapply(lapply(sims, '[[', 1), '[[', "est")

sims %>%
    lapply('[[', "clpm2") %>%
    lapply('[[', "est") %>%
    lapply('[[', 2)

## Generate CLPM Data
temp <- gen_starts(n = 10000, nwaves = 10, ri_x = 0, ri_y = 0, xr = 0, yr = 0, yx = .25, xy = .25)
temp_ar <- temp
## Generate latent X and Y intercepts
temp_st <- rmnorm(10000, varcov = matrix(c(.50, .25, .25, .50), nrow=2))
temp_both <- cbind(temp_ar, temp_st)

temp_x <- data.frame(sapply(seq(1,10), function(x) rowSums(temp_both[,c(x,21)])))
temp_y <- data.frame(sapply(seq(11,20), function(x) rowSums(temp_both[,c(x,22)])))
names(temp_y) <- paste0("y", 1:10)

data <- cbind(temp_x, temp_y)
names(data) <- tolower(names(data))

fit_ri <- lavaan(ri_clpm10_c, data = data)
summary(fit_ri)

fit_clpm <- lavaan(clpm10_c, data = temp)
summary(fit_clpm)

fit_starts <- lavaan(starts_c, data = data)
summary(fit_starts)


temp <- gen_starts(n = 10000, nwaves = 10, ri_x = 0, ri_y = 0, xr = 0, yr = 0, yx = .25, xy = .25, stab_x = .8, stab_y = .8, cor_xy = .5)
