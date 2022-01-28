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


fit_clpm <- lavaan(clpm10_c, data = clpm_data)
summary(fit_clpm)

lavaanPlot(model = fit_clpm, graph_options=list(rankdir = "LR", splines=FALSE), coefs=TRUE)


## Models with no st
fit_ri <- lavaan(ri_clpm10_c, data = clpm_data)
summary(fit_ri)

lavaanPlot(model = fit_ri, graph_options=list(rankdir = "LR", splines=FALSE), coefs=TRUE)




n=10000      # N to generate
nwaves=10   # Number of waves
ri_x=0     # Random intercept variance for X
ri_y=0     # Random intercept variance for Y
cor_i=.5   # Correlation between intercepts
x=1        # AR variance for X
y=1        # AR variance for Y
stab_x=.5  # Stability of X
stab_y=.5  # Stability of Y
yx=.1      # Cross lag (Y on X)
xy=.1      # Cross lag (X on Y)
cor_xy=.5  # Correlation between X and Y
cor_xyr=.25 # Correlation between X and Y residuals
xr=0       # Measurement error for X
yr=0       # Measurement error for Y


wxr <- 1 - ((x*stab_x^2 + y*xy^2 + stab_x*xy*cor_xy)) #x residual

((x*stab_x^2 + y*xy^2 + stab_x*xy*cor_xy))

wyr <- 1 - ((y*stab_y^2 + x*yx^2)/(1-cor_xy^2)) #y residual


n=10000      # N to generate
nwaves=10   # Number of waves
ri_x=1     # Random intercept variance for X
ri_y=1     # Random intercept variance for Y
cor_i=.5   # Correlation between intercepts
x=1        # AR variance for X
y=1        # AR variance for Y
stab_x=.5  # Stability of X
stab_y=.5  # Stability of Y
yx=.3      # Cross lag (Y on X)
xy=0      # Cross lag (X on Y)
cor_xy=.5  # Correlation between X and Y
#cor_xyr=.2 # Correlation between X and Y residuals
xr=.5       # Measurement error for X
yr=.5       # Measurement error for Y
##
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
                   xr=xr,       # Measurement error for X
                   yr=yr        # Measurement error for Y
                   )

fit_test <- lavaan(clpm2_c, data = data)

fit_ri <- lavaan(ri_clpm10_c, data = data)
summary(fit_ri)

fit_clpm <- lavaan(clpm10_c, data = data)
summary(fit_clpm)


fit_starts <- lavaan(starts_c, data = data)
summary(fit_starts)



Cxx <- matrix(c(1,.5,.5,1),nrow=2)
b <- matrix(c(.5,.2), nrow=2)
t(b) %*% Cxx %*% b


semPaths(semPlotModel_lavaanModel(ri_clpm10_c))


set.seed(1224)
clpm_data <- gen_starts(n = 10000, 
                        nwaves = 10,
                        x = 1, # Variance in Autoregressive Component for X
                        y = 1, # Variance in Autoregressive Component for Y
                        ri_x = 0, # Random Intercept Variance for X
                        ri_y = 0, # Random Intercept Variance for Y
                        xr = 0, # Measurement Error for X
                        yr = 0, # Measurement Error for Y
                        stab_x = .50, # Stability of Autoregressive Component for X
                        stab_y = .50, # Stability of Autoregressive Component for Y
                        yx = -.50, # Crosslagged path, Y regressed on X
                        xy = .00, # Crosslagged path, X regressed on Y
                        cor_xy = -.5, # Correlation between Autoregressive Component for X and Y
                        )
## Fit CLPM Model
fit_clpm <- lavaan(clpm10_c, data = clpm_data)
## Find indidivual means
clpm_means <- data.frame(x=rowMeans(clpm_data[,1:10]), y=rowMeans(clpm_data[,11:20]))
## Create mean-deviated scores
temp <- cbind(clpm_data, clpm_means)
temp_x <- sweep(temp[,1:10], MARGIN = 1, STATS = rowMeans(temp[,1:10]))
temp_y <- sweep(temp[,11:20], MARGIN = 1, STATS = rowMeans(temp[,11:20]))
md_clpm_data <- cbind(temp_x, temp_y)
## Run CLPM on mean-deviated scores
fit_clpm_md <- lavaan(clpm10_c, data = md_clpm_data)
## summary(fit_clpm_md)
## Also run RI-CLPM on Original Data for Comparison
fit_ri_clpm <- lavaan(ri_clpm10_c, data = clpm_data)

library(ggplot2)
library(ggrepel)

results <- readRDS("saved/simulation1.rds")

labels <- results %>%
  filter(N==25)

p <- ggplot(data = results, aes(x = N, y = power, group = r)) +
    geom_line() +
    scale_x_log10(breaks=c(25, 50,100,250,500,1000)) +
    facet_grid(cols = vars(reliability), rows = vars(AR_Var), labeller=label_both) +
  geom_text_repel(data=labels, aes(label=r), size=3)
p



varCovAR <- matrix(c(1,    .5,   .5^2,  .5^3,  .5^4,
                     .5,   1,    .5,  .5^2,  .5^3,
                     .5^2, .5,   1,  .5,  .5^2,
                     .5^3, .5^2, .5,  1,  .5,
                     .5^4, .5^3, .5^2,  .5,  1),
                   nrow=5, ncol=5, byrow = TRUE)

arVar <- rmnorm(n = 10000, mean=rep(1, 5), varcov=varCovAR)
stVar <- rnorm(n = 10000, mean=5, sd=1)
tempData <- data.frame(cbind(arVar, stVar))

data <- data.frame(sapply(seq(1,5), function(x) rowSums(tempData[,c(x,6)])))
names(data) <- tolower(names(data))

meanModel <- '
x =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5
#
ax1 =~ 1*x1
ax2 =~ 1*x2
ax3 =~ 1*x3
ax4 =~ 1*x4
ax5 =~ 1*x5
#
x1 ~ 0
x2 ~ 0
x3 ~ 0
x4 ~ 0
x5 ~ 0
x ~ 1
#
x ~~ x
ax1 ~~ ax1
ax2 ~~ ax2
ax3 ~~ ax3
ax4 ~~ ax4
ax5 ~~ ax5
x1 ~~ 0*x1
x2 ~~ 0*x2
x3 ~~ 0*x3
x4 ~~ 0*x4
x5 ~~ 0*x5
'

arModel <- '
x =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5
#
ax1 =~ 1*x1
ax2 =~ 1*x2
ax3 =~ 1*x3
ax4 =~ 1*x4
ax5 =~ 1*x5
#
ax2 ~ ax1
ax3 ~ ax2
ax4 ~ ax3
ax5 ~ ax4
#
x1 ~ 0
x2 ~ 0
x3 ~ 0
x4 ~ 0
x5 ~ 0
x ~ 1
ax1 ~ 1
#
x ~~ x
ax1 ~~ ax1
ax2 ~~ ax2
ax3 ~~ ax3
ax4 ~~ ax4
ax5 ~~ ax5
x1 ~~ 0*x1
x2 ~~ 0*x2
x3 ~~ 0*x3
x4 ~~ 0*x4
x5 ~~ 0*x5
'

arModel2 <- '
x =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5
#
ax1 =~ 1*x1
ax2 =~ 1*x2
ax3 =~ 1*x3
ax4 =~ 1*x4
ax5 =~ 1*x5
#
ax2 ~ ax1
ax3 ~ ax2
ax4 ~ ax3
ax5 ~ ax4
#
x ~~ x
ax1 ~~ ax1
ax2 ~~ ax2
ax3 ~~ ax3
ax4 ~~ ax4
ax5 ~~ ax5
x1 ~~ 0*x1
x2 ~~ 0*x2
x3 ~~ 0*x3
x4 ~~ 0*x4
x5 ~~ 0*x5
'


temp <- lavaan(meanModel, data=data, meanstructure=TRUE)
temp2 <- lavaan(arModel, data = data, meanstructure=TRUE)
temp3 <- lavaan(arModel2, data = data, meanstructure=FALSE)

tempAr <- lavaan(arModel, data=arVar, meanstructure=TRUE)


n=1000      # N to generate
nwaves=3   # Number of waves
ri_x=0   # Random intercept variance for X
ri_y=0     # Random intercept variance for Y
cor_i=0   # Correlation between intercepts
x=1        # AR variance for X
y=1        # AR variance for Y
stab_x=.5  # Stability of X
stab_y=.5  # Stability of Y
yx=.00      # Cross lag (Y on X)
xy=.00      # Cross lag (X on Y)
cor_xy=0  # Correlation between X and Y
#cor_xyr=.2 # Correlation between X and Y residuals
xr=0       # Measurement error for X
yr=0       # Measurement error for Y
##

set.seed(434)
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
                   xr=xr,       # Measurement error for X
                   yr=yr        # Measurement error for Y
                   )
cor(data)

meanModel <- '
x =~ 1*x1 + 1*x2 + 1*x3
#
ax1 =~ 1*x1
ax2 =~ 1*x2
ax3 =~ 1*x3
#
x ~~ x
ax1 ~~ ax1
ax2 ~~ ax2
ax3 ~~ ax3
x1 ~~ 0*x1
x2 ~~ 0*x2
x3 ~~ 0*x3
'


riModel <- '
x =~ 1*x1 + 1*x2 + 1*x3
#
ax1 =~ 1*x1
ax2 =~ 1*x2
ax3 =~ 1*x3
#
ax2 ~ ax1
ax3 ~ ax2
#
x ~~ x
ax1 ~~ ax1
ax2 ~~ ax2
ax3 ~~ ax3
x1 ~~ 0*x1
x2 ~~ 0*x2
x3 ~~ 0*x3
'


meanFit <- lavaan(meanModel, data = data)
riFit <- lavaan(riModel, data = data)
summary(meanFit)
summary(riFit)


################################################################################
## True Mean-Deviated Model???
################################################################################

n=1000      # N to generate
nwaves=10   # Number of waves
ri_x=1   # Random intercept variance for X
ri_y=1     # Random intercept variance for Y
cor_i=.5   # Correlation between intercepts
x=1        # AR variance for X
y=1        # AR variance for Y
stab_x=.5  # Stability of X
stab_y=.5  # Stability of Y
yx=0      # Cross lag (Y on X)
xy=0      # Cross lag (X on Y)
cor_xy=.5  # Correlation between X and Y
#cor_xyr=.2 # Correlation between X and Y residuals
xr=0       # Measurement error for X
yr=0       # Measurement error for Y
##

set.seed(434)
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
                   xr=xr,       # Measurement error for X
                   yr=yr        # Measurement error for Y
                   )
cor(data)


new_model <- '
## Random Intercepts
##
ri_x =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5 + 1*x6 + 1*x7 + 1*x8 + 1*x9 + 1*x10
ri_y =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 + 1*y5 + 1*y6 + 1*y7 + 1*y8 + 1*y9 + 1*y10
##
## Create within-person centered variables
  wx1 =~ 1*x1
  wx2 =~ 1*x2
  wx3 =~ 1*x3 
  wx4 =~ 1*x4
  wx5 =~ 1*x5
  wx6 =~ 1*x6
  wx7 =~ 1*x7
  wx8 =~ 1*x8 
  wx9 =~ 1*x9
  wx10 =~ 1*x10
  wy1 =~ 1*y1
  wy2 =~ 1*y2
  wy3 =~ 1*y3
  wy4 =~ 1*y4
  wy5 =~ 1*y5
  wy6 =~ 1*y6
  wy7 =~ 1*y7
  wy8 =~ 1*y8
  wy9 =~ 1*y9
  wy10 =~ 1*y10
##
## Regressions
##
##
## Stabilities
##
## None!
##
## Cross-lags
wy2 ~ c*wx1
wy3 ~ c*wx2
wy4 ~ c*wx3
wy5 ~ c*wx4
wy6 ~ c*wx5
wy7 ~ c*wx6
wy8 ~ c*wx7
wy9 ~ c*wx8
wy10 ~ c*wx9
wx2 ~ d*wy1
wx3 ~ d*wy2
wx4 ~ d*wy3
wx5 ~ d*wy4
wx6 ~ d*wy5
wx7 ~ d*wy6
wx8 ~ d*wy7
wx9 ~ d*wy8
wx10 ~ d*wy9
##
## Variances
ri_x ~~ ri_x
ri_y ~~ ri_y
wx1 ~~ wx1
wy1 ~~ wy1
wx2 ~~ wx2
wy2 ~~ wy2
wx3 ~~ wx3
wy3 ~~ wy3
wx4 ~~ wx4
wy4 ~~ wy4
wx5 ~~ wx5
wy5 ~~ wy5
wx6 ~~ wx6
wy6 ~~ wy6
wx7 ~~ wx7
wy7 ~~ wy7
wx8 ~~ wx8
wy8 ~~ wy8
wx9 ~~ wx9
wy9 ~~ wy9
wx10 ~~ wx10
wy10 ~~ wy10
##
## Covariances
ri_x ~~ ri_y
wx1 ~~ wy1
wx2 ~~ wy2
wx3 ~~ wy3
wx4 ~~ wy4
wx5 ~~ wy5
wx6 ~~ wy6
wx7 ~~ wy7
wx8 ~~ wy8
wx9 ~~ wy9
wx10 ~~ wy10
'

new_model2 <- '
## Random Intercepts
##
ri_x =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5 + 1*x6 + 1*x7 + 1*x8 + 1*x9 + 1*x10
ri_y =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 + 1*y5 + 1*y6 + 1*y7 + 1*y8 + 1*y9 + 1*y10
##
## Create within-person centered variables
  wx1 =~ 1*x1
  wx2 =~ 1*x2
  wx3 =~ 1*x3 
  wx4 =~ 1*x4
  wx5 =~ 1*x5
  wx6 =~ 1*x6
  wx7 =~ 1*x7
  wx8 =~ 1*x8 
  wx9 =~ 1*x9
  wx10 =~ 1*x10
  wy1 =~ 1*y1
  wy2 =~ 1*y2
  wy3 =~ 1*y3
  wy4 =~ 1*y4
  wy5 =~ 1*y5
  wy6 =~ 1*y6
  wy7 =~ 1*y7
  wy8 =~ 1*y8
  wy9 =~ 1*y9
  wy10 =~ 1*y10
##
## Regressions
##
##
## Stabilities
##
## None!
##
## Cross-lags
wy2 ~ c*wx1
wy3 ~ c*wx2
wy4 ~ c*wx3
wy5 ~ c*wx4
wy6 ~ c*wx5
wy7 ~ c*wx6
wy8 ~ c*wx7
wy9 ~ c*wx8
wy10 ~ c*wx9
wx2 ~ d*wy1
wx3 ~ d*wy2
wx4 ~ d*wy3
wx5 ~ d*wy4
wx6 ~ d*wy5
wx7 ~ d*wy6
wx8 ~ d*wy7
wx9 ~ d*wy8
wx10 ~ d*wy9
##
## Variances
ri_x ~~ ri_x
ri_y ~~ ri_y
wx1 ~~ wx1
wy1 ~~ wy1
wx2 ~~ wx2
wy2 ~~ wy2
wx3 ~~ wx3
wy3 ~~ wy3
wx4 ~~ wx4
wy4 ~~ wy4
wx5 ~~ wx5
wy5 ~~ wy5
wx6 ~~ wx6
wy6 ~~ wy6
wx7 ~~ wx7
wy7 ~~ wy7
wx8 ~~ wx8
wy8 ~~ wy8
wx9 ~~ wx9
wy9 ~~ wy9
wx10 ~~ wx10
wy10 ~~ wy10
##
## Covariances
ri_x ~~ ri_y
## wx1 ~~ wy1
## wx2 ~~ wy2
## wx3 ~~ wy3
## wx4 ~~ wy4
## wx5 ~~ wy5
## wx6 ~~ wy6
## wx7 ~~ wy7
## wx8 ~~ wy8
## wx9 ~~ wy9
## wx10 ~~ wy10
'


stable_trait <- '
## Random Intercepts
##
ri_x =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5 + 1*x6 + 1*x7 + 1*x8 + 1*x9 + 1*x10
##
##
## Create within-person centered variables
  wx1 =~ 1*x1
  wx2 =~ 1*x2
  wx3 =~ 1*x3 
  wx4 =~ 1*x4
  wx5 =~ 1*x5
  wx6 =~ 1*x6
  wx7 =~ 1*x7
  wx8 =~ 1*x8 
  wx9 =~ 1*x9
  wx10 =~ 1*x10
##
## Regressions
##
##
## Stabilities
##
## None!
##
## Cross-lags
##
## Variances
ri_x ~~ ri_x
wx1 ~~ wx1
wx2 ~~ wx2
wx3 ~~ wx3
wx4 ~~ wx4
wx5 ~~ wx5
wx6 ~~ wx6
wx7 ~~ wx7
wx8 ~~ wx8
wx9 ~~ wx9
wx10 ~~ wx10
##
## Covariances
'


n=10000      # N to generate
nwaves=10   # Number of waves
ri_x=1   # Random intercept variance for X
ri_y=1     # Random intercept variance for Y
cor_i=.5   # Correlation between intercepts
x=1        # AR variance for X
y=1        # AR variance for Y
stab_x=.5  # Stability of X
stab_y=.5  # Stability of Y
yx=.5      # Cross lag (Y on X)
xy=0      # Cross lag (X on Y)
cor_xy=.5  # Correlation between X and Y
#cor_xyr=.2 # Correlation between X and Y residuals
xr=1       # Measurement error for X
yr=1       # Measurement error for Y
##

set.seed(434)

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
                   xr=xr,       # Measurement error for X
                   yr=yr        # Measurement error for Y
                   )
cor(data)


test1 <- lavaan(ri_clpm10_c, data=data)
summary(test1)

test2 <- lavaan(new_model, data=data)
summary(test2)

test3 <- lavaan(stable_trait, data=data)
summary(test3)

test4 <- lavaan(new_model2, data=data)
summary(test4)

test5 <- lavaan(clpm10_c, data=data)
summary(test5)


temp <- lavaan(ri_clpm10_c, data=data)
parameterEstimates(temp)
temp2 <- lavaan(clpm10_c, data=data)
parameterEstimates(temp2)


parameterEstimates(temp) %>% filter(lhs=="wy2", op=="~", rhs=="wx1") %>% select(est)


names(resultsCl)

resultsCl %>%
    filter(Reliability==.7, Autoregressive==1,r==.1, clValue==.3)
