################################################################################
## Function to create and run mplus STARTS models
## Can include phantom variables for missing waves
## 
## "waves" = total number of possible waves (e.g., 20 for HILDA, 37 for SOEP)
## "actual waves" = waves with data (e.g., c(1,2,3,5,7,9))
################################################################################

library(MplusAutomation)

## Build trait part of input file
buildTrait <- function(waves) {
    wavesList <- 1:waves
    title <- "!!! Stable Trait;\n"
    ## First do "BY" statements
    traitModel <- paste0(title,"trait by lx", wavesList[1],"@1;\n")
    for (w in wavesList[-1]) {
        traitModel <- paste0(traitModel,
                         paste0("trait by lx", w, "@1;\n")
                         )
    }
    ## Then set variance
    traitModel <- paste0(traitModel, "\n!!! Trait Variance\ntrait (t);\n")
    return(traitModel)
}

## Build autoregressive part of input file
buildAR <- function(waves) {
    wavesList <- 1:waves
    title <- "!!! Autoregressive Part;\n"
    subtitle1 <- "\n!!! Indicator Statements;\n"
    ## First do "BY" statements linking AR to latent occasion variable
    arModel <- paste0(title, subtitle1, "ar1 by lx1@1;\n")
    for (w in wavesList[-1]) {
        arModel <- paste0(arModel, "ar", w, " by lx", w, "@1;\n")
    }
    ## Then add stability paths
    subtitle2 <- "\n!!! Regression Statements;\n"
    arModel <- paste0(arModel, subtitle2, "ar2 on ar1(a);\n")
    for (w in wavesList[-c(1:2)]) {
        arModel <- paste0(arModel,
                          paste0("ar", w, " on ar", (w-1), "(a);\n"))
    }
    ## Finally add variances (waves after the first are set to be equal)
    arModel <- paste0(arModel, "\n!!! Autoregressive Component Variance\nar1 (v);\n")
    for (w in wavesList[-1]) {
        arModel <- paste0(arModel,
                          paste0("ar", w, " (c);\n"))
    }
    return(arModel)
}

## Setup phantom variables for missing waves (if needed; returns comment if none)
buildPhantom <- function(waves, actualWaves) {
    ifelse(length(actualWaves) < length(c(1:waves)), {
        title <- "!!! Phantom Variables;\n"
        phantom <- title
        phantomWaves <- c(1:waves)[-actualWaves]
        for (w in phantomWaves) {
            phantom <- paste0(phantom, "lx", w, " by ;\n")
        }
        return(phantom)
    },
    return("!!! No Phantom Variables;\n"))
}

## Build observed part of model
buildObserved <- function(actualWaves) {
    title <- "!!! Observed Variables;\n"
    ## First do "BY" statements
    observedModel <- paste0(title, "lx", actualWaves[1], " by x", actualWaves[1], "@1;\n")
    for (w in actualWaves[-1]) {
        observedModel <- (paste0(observedModel,
                                 "lx", w, " by x", w, "@1;\n"))
    }
    ## Then constrain variance to 0 (we have a separate state variable)
    observedModel <- paste0(observedModel,
                            "\n!!! Residual Variance Constrained to 0\nx",
                            actualWaves[1], "@0;\n")
    for (w in actualWaves[-1]) {
        observedModel <- paste0(observedModel,
                                "x", w, "@0;\n")
    }
    return(observedModel)
}


## Build State variance part of model
buildState <- function(waves) {
    wavesList <- 1:waves
    title <- "\n!!! States;\n"
    ## First do "BY" statements
    stateModel <- paste0(title, "s1 by lx1@1;\n")
    for (w in wavesList[-1]) {
        stateModel <- paste0(stateModel,
                             paste0("s", w, " by lx", w, "@1;\n"))
    }
    ## Then add variances (constrained to be equal for stationarity)
    subtitle <- "\n!!! State Variance;\n"
    stateModel <- paste0(stateModel, subtitle, "s1 (s);\n")
    for (w in wavesList[-1]) {
        stateModel <- paste0(stateModel,
                             paste0("s", w, " (s);\n"))
    }
    return(stateModel)
}

## Build latent occasion variance section
buildLatentVar <- function(waves) {
    wavesList <- 1:waves
    ## All we need to do is set variances here, as we have defined these elsewhere
    title <- "\n!!! Constrain Latent Occasion Residuals to 0\n"
    latentVar <- paste0(title, "lx1@0;\n")
    for (w in wavesList[-1]) {
        latentVar <- paste0(latentVar,
                            paste0("lx", w, "@0;\n"))
    }
    return(latentVar)
}

## Build final MODEL statement
buildModel <- function(waves, actualWaves) {
    startsModel <- paste0(buildObserved(actualWaves), "\n",
                          buildPhantom(waves, actualWaves), "\n",
                          buildLatentVar(waves), "\n",
                          buildState(waves), "\n",
                          buildTrait(waves), "\n",
                          buildAR(waves), "\n" 
                          )
    cat(startsModel)
    return(startsModel)
}

run_starts_mplus <- function(data, waves, actualWaves, title="test", output="test") {
    ## Build Constraint Statement
    constraints <- paste("c = v - v*a*a;",
                         "NEW(traitVar, arVar, stateVar);",
                         "traitVar=t/(t+v+s);",
                         "arVar=v/(t+v+s);",
                         "stateVar=s/(t+v+s);", sep="\n")
    startsModel <- buildModel(waves, actualWaves)
    inp <- mplusObject(TITLE=title,
                       rdata=data,
                       usevariables = paste0("x", actualWaves),
                       ANALYSIS = "MODEL=NOCOVARIANCES;",
                       MODEL=startsModel,
                       MODELCONSTRAINT = constraints
                       )
    
    mplusModeler(inp, modelout = paste0(output, ".inp"), run=1)
}




#run_starts_mplus(data, 10, c(1:10), title="New Test", output="newTest")


#source("scripts/gen_starts.R")
#data <- gen_starts(n=10000, nwaves=10, yx=0, xy=0, xr=1, yr=1)



#mplusModeler(inp, modelout = "test.inp", run=1)

## ## Test with lavaan
## library(lavaan)
## lStarts <- lavaan(starts_uni, data=data)
## summary(lStarts)
