library(lavaan)
library(tidyverse)
library(knitr)
library(rethinking)
library(papaja)

## Load scripts and models
source("scripts/gen_starts.R") ## Generate data
source("scripts/clpm2_c.R") ## Lavaan model for 2-wave clpm with constraints
source("scripts/clpm3_c.R") ## Lavaan model for 3-wave clpm with constraints
source("scripts/clpm5_c.R") ## Lavaan model for 5-wave clpm with constraints
source("scripts/clpm10_c.R") ## Lavaan model for 10-wave clpm with constraints
source("scripts/ri_clpm3_c.R") ## Lavaan model for 3-wave ri-clpm with constraints
source("scripts/ri_clpm10_c.R") ## Lavaan model for 10-wave ri-clpm with constraints
source("scripts/starts_c.R") ## Lavaan model for 10-wave starts with constraints
source("scripts/run_sim.R") ## Script to run simulations

data <- gen_starts(n=10000, cor_i = 0, yx = 0, xy = 0, cor_xy=.7)

test2 <- lavaan(clpm2_c, data = data)
summary(test2)

test5 <- lavaan(clpm5_c, data = data)
summary(test5)

test10 <- lavaan(clpm10_c, data = data)
summary(test10)
