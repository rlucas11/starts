library(tidyverse)
library(lavaan)
source("scripts/gen_starts.R")

## Generate AR Data
data <- gen_starts(
    n = 1000000,
    nwaves = 3,
    ri_x = 0,
    ri_y = 0,
    cor_i = 0,
    x = 1,
    y = 1,
    stab_x = .5,
    stab_y = .5,
    yx = 0,
    xy = 0,
    cor_xy = 0,
    xr = 0,
    yr = 0
)

## Select X variables
data.x <- data |>
    select(starts_with("x")) |>
    rowwise() |>
    mutate(x.mean = mean(c(x1, x2, x3)))
cor(data.x)

lvmodel <- "
ri =~ 1*x1 + 1*x2 + 1*x3
"

lvFit <- sem(lvModel, data = data.x)
summary(lvFit)


riModel <- "
ri =~ 1*x1 + 1*x2 + 1*x3
x2 ~ x1
x3 ~ x2
"
riFit <- sem(riModel, data = data.x)
summary(riFit)




lvModelCor <- "
ri = ~ 1 * x1 + 1 * x2 + 1 * x3
ri ~~ x.mean
"




lvModelCorFit <- cfa(lvModelCor, data = data.x)

summary(lvModelCorFit)
standardizedSolution((lvModelCorFit))

## Generate Random Variable
ri <- rnorm(10000)
