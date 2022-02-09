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


