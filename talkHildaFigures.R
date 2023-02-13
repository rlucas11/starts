library(tidyverse)
library(ggplot2)
source("scripts/usefulFunctions.R")

commute <- read_csv("data/commute.csv")
commuteWideLog <- commute %>%
    mutate(logCommute = log(commute + 1)) %>%
    select(xwaveid, year, logCommute) %>%
    pivot_wider(
        names_from = year,
        names_prefix = "commute",
        values_from = logCommute
    )
commuteR <- data.frame(
    Lag = 1:18,
    r = summarizeR(cor(commuteWideLog[, -1], use = "pair")),
    Variable = "Commute"
)

income <- read_csv("data/income.csv")
incHWide <- income %>%
    mutate(logInc = log(income + 1)) %>%
    select(xwaveid, hhid, year, logInc) %>%
    group_by(year, hhid) %>%
    filter(xwaveid == min(xwaveid)) %>%
    ungroup() %>%
    select(xwaveid, year, logInc) %>%
    arrange(year, xwaveid) %>%
    pivot_wider(names_from = year, names_prefix = "income", values_from = logInc)
incomeR <- data.frame(
    Lag = 1:19,
    r = summarizeR(cor(incHWide[, -1], use = "pair")),
    Variable = "Income"
)


## Household Wages
wages <- read_csv("data/wages.csv")
wHWide <- wages %>%
    mutate(logW = log(wages + 1)) %>%
    select(xwaveid, hhid, year, logW) %>%
    group_by(year, hhid) %>%
    filter(xwaveid == min(xwaveid)) %>%
    ungroup() %>%
    select(xwaveid, year, logW) %>%
    arrange(year, xwaveid) %>%
    pivot_wider(names_from = year, names_prefix = "wages", values_from = logW)
wagesR <- data.frame(
    Lag = 1:19,
    r = summarizeR(cor(wHWide[, -1], use = "pair")),
    Variable = "Wages"
)


#### Life Satisfation
## Get Data
ls <- read_csv("data/ls.csv")
lsWide <- ls %>%
    pivot_wider(names_from = year, names_prefix = "ls", values_from = ls)
lsR <- data.frame(
    Lag = 1:19,
    r = summarizeR(cor(lsWide[, -1], use = "pair")),
    Variable = "LS"
)

#### sf-36 pain
sf36Pain <- read_csv("data/sf36Pain.csv")
sf36PainWide <- sf36Pain %>%
    pivot_wider(names_from = year, names_prefix = "sf36Pain", values_from = sf36Pain)
painR <- data.frame(
    Lag = 1:19,
    r = summarizeR(cor(sf36PainWide[, -1], use = "pair")),
    Variable = "Pain"
)

#### physical activity
physical <- read_csv("data/physical.csv")
physicalWide <- physical %>%
    pivot_wider(names_from = year, names_prefix = "physical", values_from = physical)
physicalR <- data.frame(
    Lag = 1:19,
    r = summarizeR(cor(physicalWide[, -1], use = "pair")),
    Variable = "Activity"
)

#### Feeling rushed or pressed for time
rushed <- read_csv("data/rushed.csv")
rushedWide <- rushed %>%
    pivot_wider(names_from = year, names_prefix = "rushed", values_from = rushed)
rushedR <- data.frame(
    Lag = 1:19,
    r = summarizeR(cor(rushedWide[, -1], use = "pair")),
    Variable = "Rushed"
)

#### Weight
weight <- read_csv("data/weight.csv")
weightWide <- weight %>%
    pivot_wider(names_from = year, names_prefix = "weight", values_from = weight)
weightR <- data.frame(
    Lag = 1:14,
    r = summarizeR(cor(weightWide[, -1], use = "pair")),
    Variable = "Weight"
)

#### Social Support
ss <- read_csv("data/ss.csv")
ssWide <- ss %>%
    pivot_wider(names_from = year, names_prefix = "ss", values_from = ss)
supportR <- data.frame(
    Lag = 1:19,
    r = summarizeR(cor(ssWide[, -1], use = "pair")),
    Variable = "Support"
)

#### General Health
genHealth <- read_csv("data/genHealth.csv")
genHealthWide <- genHealth %>%
    pivot_wider(
        names_from = year,
        names_prefix = "genHealth",
        values_from = genHealth
    )
healthR <- data.frame(
    Lag = 1:19,
    r = summarizeR(cor(genHealthWide[, -1], use = "pair")),
    Variable = "Health"
)


data <- rbind(
    commuteR,
    incomeR,
    wagesR,
    lsR,
    physicalR,
    painR,
    rushedR,
    weightR,
    supportR,
    healthR
)

write_csv(data, "saved/hildaCors.csv")
