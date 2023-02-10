## Load packages
library(lavaan) ## Run models
library(tidyverse) ## Data Manipulation
library(rethinking) ## Parallel Simulations
library(ggplot2)
library(ggrepel)
library(directlabels)

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


results5 <- readRDS("saved/5WaveSimulation.rds")
## Changes names for plot
names(results5) <- c("N", "r","Reliability", "Autoregressive", "power", "estimatex","estimatey")
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

results5 %>%
    filter(N>25, Reliability==.7, Autoregressive==1) %>%
    ggplot(aes(x = N, y = power, group = r)) +
    scale_x_continuous(breaks=c(50,100,250,500,1000)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) +
    theme_bw() +
    theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(linetype = "Correlation Between Stable Traits") +
    ylab("Probability of One or More Spurious Cross-Lagged Effects")

ggsave("images/grantFigure2.png", width=5, height=5, units="in")




################################################################################
## 10-Wave Simulation
################################################################################

## Actual Simulation Values
nValues <- c(25,50,100,250,500,1000)
rValues <- c(.1,.3,.5,.7)
reliabilities <- c(.7)
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
                sims <- data.frame(t(mcreplicate(n=1000, run_sim_clpm(waves = 10,
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

saveRDS(results, "saved/10WaveSimulation.rds")

results10 <- readRDS("saved/10WaveSimulation.rds")
## Changes names for plot
names(results10) <- c("N", "r","Reliability", "Autoregressive", "power", "estimatex","estimatey")
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

results10 %>%
    filter(N>25, Reliability==.7, Autoregressive==1) %>%
    ggplot(aes(x = N, y = power, group = r)) +
    scale_x_continuous(breaks=c(50,100,250,500,1000)) +
    geom_line(aes(linetype=as.factor(r)),color="black", size=.5) +
    theme_bw() +
    theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(linetype = "Correlation Between Stable Traits") +
    ylab("Probability of One or More Spurious Cross-Lagged Effects")

ggsave("images/grantFigure10.png", width=5, height=5, units="in")


################################################################################
## Hilda Data
################################################################################

nYears <- 10

################################################################################
## Plot AR-1 Data
################################################################################
## Function to calculate decline in correlations as lag increases
arTrend <- function(x, n) {
    corVec <- vector(length=n)
    for (i in 1:n) {
        corVec[i] <- x^i
    }
    return(corVec)
}

data <- data.frame(
    Years = rep(1:nYears, 3),
    Correlation = c(
        arTrend(.75, nYears),
        arTrend(.50, nYears),
        arTrend(.25, nYears)
    ),
    Stability = rep(c(.75, .50, .25), each = nYears)
)

arPlot <- ggplot(aes(x=Years, y=Correlation, group=Stability),
                 data=data) +
    geom_line(lwd=.2) +
        theme_grey(base_size = 8) +
        ylim(-.1, 1) +
        scale_x_continuous(labels=c(1:10), breaks=c(1:10)) +
    theme(
        panel.background = element_rect(fill='transparent', color=NA),
        plot.background = element_rect(fill='transparent', color=NA),
		axis.line = element_line(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    ggtitle('Panel A')

################################################################################
## Plot AR-1 Data with Cross-Lagged Paths
################################################################################
source("~/Projects/starts/scripts/usefulFunctions.R")
source("~/Projects/starts/scripts/gen_starts.R")

cl <- .20 ## 90th Percentile, from Orth
dataCL.1 <- data.frame(
    Years = rep(1:nYears, 3),
    Correlation = c(
        summarizeR(cor(gen_starts(
            n = 10000,
            nwaves = nYears+1,
            ri_x = 0,
            ri_y = 0,
            cor_i = 0,
            x = 1,
            y = 1,
            stab_x = .73,
            stab_y = .73,
            yx = cl,
            xy = cl,
            cor_xy = .5,
            xr = 0,
            yr = 0
        )[, 1:11])),
        summarizeR(cor(gen_starts(
            n = 10000,
            nwaves = nYears+1,
            ri_x = 0,
            ri_y = 0,
            cor_i = 0,
            x = 1,
            y = 1,
            stab_x = .5,
            stab_y = .5,
            yx = cl,
            xy = cl,
            cor_xy = .5,
            xr = 0,
            yr = 0
        )[, 1:11])),
        summarizeR(cor(gen_starts(
            n = 10000,
            nwaves = nYears+1,
            ri_x = 0,
            ri_y = 0,
            cor_i = 0,
            x = 1,
            y = 1,
            stab_x = .25,
            stab_y = .25,
            yx = cl,
            xy = cl,
            cor_xy = .5,
            xr = 0,
            yr = 0
        )[, 1:11]))
    ),
    Stability = rep(c(.75, .50, .25), each = nYears)
)

clPlot <- ggplot(aes(x=Years, y=Correlation, group=Stability),
                 data=dataCL.1) +
    geom_line(lwd=.2) +
        theme_grey(base_size = 8) +
        ylim(-.1, 1) +
        scale_x_continuous(labels=c(1:10), breaks=c(1:10)) +
    theme(
        panel.background = element_rect(fill='transparent', color=NA),
        plot.background = element_rect(fill='transparent', color=NA),
		axis.line = element_line(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        ## axis.text.y=element_blank(),
        ## axis.title.y=element_blank()
    ) 

stab73 <- dataCL.1 %>%
    filter(Stability == .75) %>%
    select(
        Years,
        Correlation
    ) %>%
    mutate(Source = 1,
           Variable = "") 


################################################################################
## Plot Actual Data from HILDA
################################################################################

hilda <- read_csv("data/hildaCors.csv")
names(hilda) <- c("Years", "Correlation", "Variable")
hilda$Source <- 0

hilda2 <- rbind(hilda, stab73)

hildaPlot <- ggplot(aes(x=Years, y=Correlation, fill=Variable),
                 data=hilda[hilda$Years <= 10,]) +
    geom_line(lwd=.2) +
	theme_grey(base_size=20) +
    theme(
        panel.background = element_rect(fill='transparent', color=NA),
        plot.background = element_rect(fill='transparent', color=NA),
		axis.line = element_line(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'),
        ## axis.text.y=element_blank(),
        ## axis.title.y=element_blank(),
    ) +
    ylim(-.1, 1) +
#    xlim(1,12) +
    scale_x_continuous(labels=c(1:10), breaks=c(1:10), limits=c(1,11)) +
    scale_color_manual(values = rep("black", 10))

my.dl <- list(fill="white", "draw.rects")
direct.label(
    hildaPlot,
    list(
        cex = .6,
        "far.from.others.borders",
        "calc.boxes",
        "enlarge.box",
        "my.dl"
    )
)

ggsave("images/hildaGrant.png", scale=8)


## With Generated Data
hildaPlot2 <- ggplot(aes(
    x=Years,
    y=Correlation,
    fill=Variable,
    linetype = factor(Source)
),
data=hilda2[hilda2$Years <= 10,]) +
    geom_line(lwd=.2) +
    theme_grey(base_size=20) +
    theme(
        panel.background = element_rect(fill='transparent', color=NA),
        plot.background = element_rect(fill='transparent', color=NA),
        axis.line = element_line(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'),
        ## axis.text.y=element_blank(),
        ## axis.title.y=element_blank(),
    ) +
    ylim(-.1, 1) +
        #    xlim(1,12) +
        scale_x_continuous(labels = c(1:10), breaks = c(1:10), limits = c(1, 11)) +
        scale_color_manual(values = rep("black", 10))
hildaPlot2    

my.dl <- list(fill="white", "draw.rects")
direct.label(
    hildaPlot2,
    list(
        cex = .6,
        "far.from.others.borders",
        "calc.boxes",
        "enlarge.box",
        "my.dl"
    )
)



################################################################################
## Combine Plots
################################################################################

grid.arrange(arPlot,
             clPlot,
##             direct.label(hildaPlot, list(cex=.3, "angled.boxes")),
             direct.label(hildaPlot, list(cex=.3,
                             "far.from.others.borders",
                             "calc.boxes",
                             "enlarge.box",
                             "my.dl")),
             nrow=1)


################################################################################
## Revision Figures
################################################################################

hildaR <- hilda %>%
    filter(Variable == "LifeSat" |
           Variable == "Activity" |
           Variable == "Pain" |
           Variable == "Health" |
           Variable == "Weight")


## Generate corresponding CLPM data
cl <- .20 ## 90th Percentile, from Orth
dataCL.1 <- data.frame(
    Years = rep(1:nYears, 5),
    Correlation = c(
        summarizeR(cor(gen_starts(
            n = 10000,
            nwaves = nYears+1,
            ri_x = 0,
            ri_y = 0,
            cor_i = 0,
            x = 1,
            y = 1,
            stab_x = .60,
            stab_y = .60,
            yx = cl,
            xy = cl,
            cor_xy = .5,
            xr = 0,
            yr = 0
        )[, 1:11])),
        summarizeR(cor(gen_starts(
            n = 10000,
            nwaves = nYears+1,
            ri_x = 0,
            ri_y = 0,
            cor_i = 0,
            x = 1,
            y = 1,
            stab_x = .57,
            stab_y = .57,
            yx = cl,
            xy = cl,
            cor_xy = .5,
            xr = 0,
            yr = 0
        )[, 1:11])),
        summarizeR(cor(gen_starts(
            n = 10000,
            nwaves = nYears+1,
            ri_x = 0,
            ri_y = 0,
            cor_i = 0,
            x = 1,
            y = 1,
            stab_x = .64,
            stab_y = .64,
            yx = cl,
            xy = cl,
            cor_xy = .5,
            xr = 0,
            yr = 0
        )[, 1:11])),
        summarizeR(cor(gen_starts(
            n = 10000,
            nwaves = nYears+1,
            ri_x = 0,
            ri_y = 0,
            cor_i = 0,
            x = 1,
            y = 1,
            stab_x = .95,
            stab_y = .95,
            yx = .03, ## Had to adjust
            xy = .03,
            cor_xy = .5,
            xr = 0,
            yr = 0
            )[, 1:11])),
        summarizeR(cor(gen_starts(
            n = 10000,
            nwaves = nYears+1,
            ri_x = 0,
            ri_y = 0,
            cor_i = 0,
            x = 1,
            y = 1,
            stab_x = .72,
            stab_y = .72,
            yx = cl,
            xy = cl,
            cor_xy = .5,
            xr = 0,
            yr = 0
            )[, 1:11]))
    ),
    Stability = rep(c(.60, .57, .64, .95, .72), each = nYears)
)

stabCl <- dataCL.1 %>%
    select(
        Years,
        Correlation
    ) %>%
    mutate(Source = 1)
stabCl$Variable <- rep(c("LifeSat", "Activity", "Pain", "Health", "Weight"), each=nYears)

hildaR <- rbind(hildaR, stabCl)


hildaPlot <- ggplot(aes(x=Years, y=Correlation, fill=Variable, linetype=factor(Source)),
                 data=hildaRf[hildaRf$Years <= 15,]) +
    geom_line(lwd=.2) +
	theme_grey(base_size=20) +
    theme(
        panel.background = element_rect(fill='transparent', color=NA),
        plot.background = element_rect(fill='transparent', color=NA),
		axis.line = element_line(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.position = "none"
        ## axis.text.y=element_blank(),
        ## axis.title.y=element_blank(),
    ) +
    ylim(-.1, 1) +
#    xlim(1,12) +
    scale_x_continuous(labels=c(1:10), breaks=c(1:10), limits=c(1,11)) +
    scale_color_manual(values = rep("black", 10))

my.dl <- list(fill="white", "draw.rects")
direct.label(
    hildaPlot,
    list(
        cex = .6,
        "far.from.others.borders",
        "calc.boxes",
        "enlarge.box",
        "my.dl"
    )
)
