library(tidyverse)
library(haven)
library(psych)


################################################################################
## Set basic info
################################################################################

## Because the HILDA can't be shared, and because it is used for many projects,
## I keep it in a separate directory. Set path to the data here:
rawPath <- file.path("/home","rich","data","HILDA","HILDA20")

## Basic details of HILDA version
maxWave <- 20
waveLetters <- letters[1:maxWave]
## List of files to read
in_files <- list.files(path=rawPath, pattern="Rperson", full.names=TRUE)
in_files_h <- list.files(path=rawPath, pattern="Household", full.names = TRUE)


## List of waves and corresponding years
yearsWaves <- data.frame(waveNo = 1:maxWave,
                         wave = letters[1:maxWave],
                         year = 2001:(2001+(maxWave-1)),
                         file = in_files,
                         h_file = in_files_h)


################################################################################
## Get Master List
################################################################################
if (!file.exists("data/master.csv")) {
  master <- read_dta(file.path(rawPath,
                               paste0("Master_", letters[maxWave],
                                      maxWave, "0c.dta")))
    master <- master %>%
      select(xwaveid, sex, yrenter, ends_with("status"),
             ends_with("hhrhid")) %>%
      filter(if_any(contains("status"), ~. == 1 | . == 2))
               write_csv(master, "data/master.csv")
}
master <- read_csv("data/master.csv")



################################################################################
## Function to extract variables (modified for HILDA)
################################################################################

extractHilda <- function(varNameStem, newName, waves) {
        long  <- data.frame(xwaveid = character(),
                            year = numeric(),
                            value = numeric())
        names(long) <- c("xwaveid", "year", newName)
        for (i in waves) {
            currentVarName <- paste0(i, varNameStem)
            currentYear <- yearsWaves[yearsWaves$wave==i, "year"]
            current <- read_dta(yearsWaves[yearsWaves$wave==i, "file"])
            new <- data.frame(xwaveid = current$xwaveid,
                              year = currentYear,
                              value = current[,currentVarName])
            names(new) <- c("xwaveid", "year", newName)
            long <- bind_rows(long, new)
        }
        return(long)
}


## Extract data from household files (must link to individual)
extractHildaH <- function(varNameStem, newName, waves) {
  master <- read_csv("data/master.csv", show_col_types = FALSE)
  long <- data.frame(
    xwaveid = character(),
    hhid = numeric(),
    year = numeric(),
    value = numeric()
  )
  names(long) <- c("xwaveid", "hhid", "year", newName)
  for (i in waves) {
    currentVarName <- paste0(i, varNameStem)
    currentYear <- yearsWaves[yearsWaves$wave == i, "year"]
    currentHid <- paste0(i, "hhrhid")
    currentStatus <- paste0(i, "fstatus")
    current <- read_dta(yearsWaves[yearsWaves$wave == i, "h_file"]) %>%
      select({{ currentHid }}, {{ currentVarName }}) %>%
      mutate({{ currentHid }} := as.numeric(.data[[currentHid]]))
    link <- master %>%
      filter(.data[[currentStatus]] == 1 | .data[[currentStatus]] == 2) %>%
      mutate(
        year = {{ currentYear }},
        {{ currentHid }} := as.numeric(.data[[currentHid]])
      ) %>%
      left_join(current, by = {{ currentHid }}) %>%
      select(xwaveid, {{ currentHid }}, year, {{ currentVarName }})
    names(link) <- c("xwaveid", "hhid", "year", newName)
    long <- bind_rows(long, link)
  }
  return(long)
}

## (requires psych library)
extractHildaMultiple <- function(varNameStems, newNames, waves) {
    counter <- 1
    for (i in waves) {
        currentVarNames <- paste0(i, varNameStems)
        currentYear <- yearsWaves[yearsWaves$wave==i, "year"]
        current <- read_dta(yearsWaves[yearsWaves$wave==i, "file"])
        if (counter == 1) {
            long <- current[, c("xwaveid", currentVarNames)]
            names(long) <- c("xwaveid", newNames)
            long$year <- currentYear
        } else {
            new <- current[, c("xwaveid", currentVarNames)]
            names(new) <- c("xwaveid", newNames)
            new$year <- currentYear
            long <- bind_rows (long, new)
        }
        counter <- counter+1
    }
    return(long)
}



################################################################################
## Life Satisfaction (updated 2022)
################################################################################
if(!file.exists("data/ls.csv")) {
    ls <- extractHilda("losat", "ls", letters[1:20])
    ls <- mutate(ls, ls = replace(ls, which(ls<0), NA))
    ## Change format of year to deal with readr
    ls$year <- as.integer(ls$year) ## Needed so write_csv writes this correctly
    write_csv(ls, "data/ls.csv")
}
ls <- read_csv("data/ls.csv")
    
lsWide <- ls %>%
    pivot_wider(names_from = year, names_prefix = "ls", values_from = ls)


################################################################################
## Weight (updated 2022)
################################################################################

if(!file.exists("data/weight.csv")) {
    weight <- extractHilda("bmwtkg", "weight", letters[6:20])
    weight <- mutate(weight, weight = replace(weight, which(weight<0), NA))
    weight$year <- as.integer(weight$year)
    write_csv(weight, "data/weight.csv")
}
weight <- read_csv("data/weight.csv")

weightWide <- weight %>%
    pivot_wider(names_from = year, names_prefix = "weight", values_from = weight)


################################################################################
## Social Support (updated 2022)
################################################################################

if(!file.exists("data/ss.csv")) {
    ss_list=c("lssuppv","lssupnh","lssuplf","lssupac","lssuplt","lssupcd","lssupvl","lssuppi","lssuptp","lssupsh")
    ss <- extractHildaMultiple(ss_list, ss_list, waves=letters[1:20])
    ss <- ss %>%
        mutate_at(ss_list, ~ replace(., .<0, NA))
    ssKeys <- list(ss = c("-lssuppv", "-lssupnh", "lssuplf", "-lssupac", "-lssuplt", "lssupcd",
                          "-lssupvl", "lssuppi", "lssuptp", "lssupsh"))
    scaleScores <- scoreItems(ssKeys, ss, impute="none")
    ss$ss <- scaleScores$scores
    ss$year <- as.integer(ss$year)
    ss$ss <- as.numeric(ss$ss)
    ss <- select(ss, xwaveid, year, ss)
    write_csv(ss, "data/ss.csv")
}
ss <- read_csv("data/ss.csv")

ssWide <- ss %>%
    pivot_wider(names_from = year, names_prefix = "ss", values_from = ss)



################################################################################
## Health (updated 2022)
################################################################################

if (!file.exists("data/genHealth.csv")) {
    genHealth <- extractHilda("gh1", "genHealth", letters[1:20])
    genHealth <- mutate(genHealth,
                        genHealth = replace(genHealth, which(genHealth<0), NA))
    genHealth$year <- as.integer(genHealth$year)
    write_csv(genHealth, "data/genHealth.csv")
}
genHealth <- read_csv("data/genHealth.csv")

genHealthWide <- genHealth %>%
  pivot_wider(names_from = year,
              names_prefix = "genHealth",
              values_from = genHealth)



################################################################################
## Income (gross regular income) (updated 2022)
################################################################################

if (!file.exists("data/income.csv")) {
    income <- extractHildaH("hifefp", "income", letters[1:20])
    income <- mutate(income, income = replace(income, which(income<0), NA))
    income$year <- as.integer(income$year)
    write_csv(income, "data/income.csv")
}
income <- read_csv("data/income.csv")


################################################################################
## Wages (updated 2022)
################################################################################

if (!file.exists("data/wages.csv")) {
    wages <- extractHildaH("hiwsfei", "wages", letters[1:20])
    wages <- mutate(wages, wages = replace(wages, which(wages<0), NA))
    wages$year <- as.integer(wages$year)
    write_csv(wages, "data/wages.csv")
}
wages <- read_csv("data/wages.csv")


################################################################################
## SF-36 Bodily Pain (2022)
################################################################################

if (!file.exists("data/sf36Pain.csv")) {
    sf36Pain <- extractHilda("ghbp", "sf36Pain", letters[1:20])
    sf36Pain <- mutate(sf36Pain, sf36Pain = replace(sf36Pain, which(sf36Pain<0), NA))
    sf36Pain$year <- as.integer(sf36Pain$year)
    write_csv(sf36Pain, "data/sf36Pain.csv")
}
sf36Pain <- read_csv("data/sf36Pain.csv")


################################################################################
## Physical Activity (2022)
################################################################################

if (!file.exists("data/physical.csv")) {
    physical <- extractHilda("lspact", "physical", letters[1:20])
    physical <- mutate(physical,
                       physical = replace(physical, which(physical<0), NA))
    physical$year <- as.integer(physical$year)
    write_csv(physical, "data/physical.csv")
}
physical <- read_csv("data/physical.csv")

################################################################################
## Commuting (2022)
################################################################################

if (!file.exists("data/commute.csv")) {
    commute <- extractHilda("lshrcom", "commute", letters[2:20])
    commute <- mutate(commute, commute = replace(commute, which(commute<0), NA))
    commute$year <- as.integer(commute$year)
    write_csv(commute, "data/commute.csv")
}
commute <- read_csv("data/commute.csv")



################################################################################
## Feeling Rushed for Time (2022)
################################################################################

if (!file.exists("data/rushed.csv")) {
    rushed <- extractHilda("lsrush", "rushed", letters[1:20])
    rushed <- mutate(rushed, rushed = replace(rushed, which(rushed < 0), NA))
    rushed$year <- as.integer(rushed$year)
    write_csv(rushed, "data/rushed.csv")
}
rushed <- read_csv("data/rushed.csv")


################################################################################
## STARTS Models for Extracted Variables
################################################################################

library(STARTS)
source("scripts/usefulFunctions.R")

#### commute
commute <- read_csv("data/commute.csv")
commuteWideLog <- commute %>%
  mutate(logCommute = log(commute+1)) %>%
  select(xwaveid, year, logCommute) %>%
  pivot_wider(names_from = year,
              names_prefix = "commute",
              values_from = logCommute)
commuteStarts <- starts_uni_estimate(commuteWideLog[, c(2:20)])
summary(commuteStarts)

#### income (household level)
incHWide <- income %>%
    mutate(logInc = log(income+1)) %>%
    select(xwaveid, hhid, year, logInc) %>%
    group_by(year, hhid) %>%
    filter(xwaveid==min(xwaveid))%>%
    ungroup() %>%
    select(xwaveid, year, logInc) %>%
    arrange(year, xwaveid) %>%
    pivot_wider(names_from = year, names_prefix = "income", values_from = logInc)

## plotCors(summarizeR(cor(incHWide[,-1], use="pair")))

incStartsH <- starts_uni_estimate(incHWide[,c(2:21)])
summary(incStartsH)

## Household Wages
wages <- read_csv("data/wages.csv")
wHWide <- wages %>%
    mutate(logW = log(wages+1)) %>%
    select(xwaveid, hhid, year, logW) %>%
    group_by(year, hhid) %>%
    filter(xwaveid==min(xwaveid))%>%
    ungroup() %>%
    select(xwaveid, year, logW) %>%
    arrange(year, xwaveid) %>%
    pivot_wider(names_from = year, names_prefix = "wages", values_from = logW)

## plotCors(summarizeR(cor(wHWide[,-1], use="pair")))

wageStartsH <- starts_uni_estimate(wHWide[,-1])
summary(wageStartsH)

## wageStartsH2 <- starts_uni_estimate(wHWide[,-1], estimator = "MCMC")

#### Life Satisfation
## Get Data
ls <- read_csv("data/ls.csv")
lsWide <- ls %>%
    pivot_wider(names_from = year, names_prefix = "ls", values_from = ls)

lsStarts <- starts_uni_estimate(lsWide[, c(2:21)])
summary(starts1)

#### sf-36 pain
sf36Pain <- read_csv("data/sf36Pain.csv")
sf36PainWide <- sf36Pain %>%
    pivot_wider(names_from = year, names_prefix = "sf36Pain", values_from = sf36Pain)
sf36PainStarts <- starts_uni_estimate(sf36PainWide[, c(2:21)])
summary(sf36PainStarts)


#### physical activity
physical <- read_csv("data/physical.csv")
physicalWide <- physical %>%
    pivot_wider(names_from = year, names_prefix = "physical", values_from = physical)
physicalStarts <- starts_uni_estimate(physicalWide[, c(2:21)])
summary(physicalStarts)

#### Feeling rushed or pressed for time
rushed <- read_csv("data/rushed.csv")
rushedWide <- rushed %>%
    pivot_wider(names_from = year, names_prefix = "rushed", values_from = rushed)
rushedStarts <- starts_uni_estimate(rushedWide[, c(2:21)])
summary(rushedStarts)

#### Weight
weight <- read_csv("data/weight.csv")
weightWide <- weight %>%
    pivot_wider(names_from = year, names_prefix = "weight", values_from = weight)
weightStarts <- starts_uni_estimate(weightWide[, c(2:16)])
summary(weightStarts)


#### Social Support
ss <- read_csv("data/ss.csv")
ssWide <- ss %>%
    pivot_wider(names_from = year, names_prefix = "ss", values_from = ss)
ssStarts <- starts_uni_estimate(ssWide[, c(2:21)])


#### General Health
genHealth <- read_csv("data/genHealth.csv")
genHealthWide <- genHealth %>%
  pivot_wider(names_from = year,
              names_prefix = "genHealth",
              values_from = genHealth)
genHealthStarts <- starts_uni_estimate(genHealthWide[, c(2:21)])


decompTable <- matrix(nrow = 10, ncol = 4)
decompTable[1, ] <- c("Life Satisfaction", lsStarts$var_prop$est)
decompTable[2, ] <- c("Social Support", ssStarts$var_prop$est)
decompTable[3, ] <- c("General Health", genHealthStarts$var_prop$est)
decompTable[4, ] <- c("SF-36 Pain", sf36PainStarts$var_prop$est)
decompTable[5, ] <- c("Weight", weightStarts$var_prop$est)
decompTable[6, ] <- c("Physical Activity", physicalStarts$var_prop$est)
decompTable[7, ] <- c("Pressed for Time", rushedStarts$var_prop$est)
decompTable[8, ] <- c("Household Income", incStartsH$var_prop$est)
decompTable[9, ] <- c("Household Wages", wageStartsH$var_prop$est)
decompTable[10, ] <- c("Minutes Commuting", commuteStarts$var_prop$est)


decompTable <- as.data.frame(decompTable)
names(decompTable) <- c("Variable", "Stable Trait", "Autoregressive Trait", "State")

write_csv(decompTable, "saved/hildaDecompTable.csv")
