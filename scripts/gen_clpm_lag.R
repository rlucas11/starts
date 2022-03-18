library(mnormt)
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
source("scripts/clpm10_l2.R")

## Testing
n=10000
nwaves=10
x=1
y=1
stab_x1=.4
stab_x2=.2
stab_y1=.4
stab_y2=.2
yx=.1
xy=.1
yx2=.1
xy2=.1
x=1
y=1
cor_xy=.5
x_rel=.9
y_rel=.9
xr=0
yr=0

gen_clpm <- function(n=10000,
                     nwaves=10,
                     stab_x1=.5,
                     stab_x2=.25,
                     stab_y1=.5,
                     stab_y2=.25,
                     yx=0, #y predicted from x
                     xy=0, #x predicted from y
                     x=1,
                     y=1,
                     cor_xy=.5,
                     xr=0,
                     yr=0
                     ) {

    ## Stationarity Constraints (not accurate yet)
    ##ifelse(x==0, wxr <- 0, wxr <- (1-stab_x^2)*x - 2*stab_x*xy*cor_xy - xy^2*y)
    ##ifelse(y==0, wyr <- 0, wyr <- (1-stab_y^2)*y - 2*stab_y*yx*cor_xy - yx^2*x)
    wxr2 <- (1-stab_x1^2)*x -
        xy^2*y -
        2*stab_x1*xy*cor_xy 
    wxr3 <- (x * ((1 + stab_x2) * (1 - stab_x1 - stab_x2) * (1 + stab_x1 - stab_x2))) / (1 - stab_x2) -
        xy^2*y -
        2*stab_x1*xy*cor_xy -
        2*stab_x1*xy*cor_xy*stab_x1 -
        2*stab_x2*x*yx*y*xy
    
    wyr2 <- (1-stab_y1^2)*y -
        yx^2*x -
        2*stab_y1*yx*cor_xy 
    wyr3 <- (y * ((1 + stab_y2) * (1 - stab_y1 - stab_y2) * (1 + stab_y1 - stab_y2))) / (1 - stab_y2) -
        yx^2*x -
        2*stab_y1*yx*cor_xy -
        2*stab_y1*yx*cor_xy*stab_y1 -
        2*stab_y2*y*xy*x*yx
    wyr4 <- wyr3

##    alt_cor_xyr <- cor_xy * (1 - stab_x1*stab_y1 - yx*xy) - (x * stab_x1 * yx) - (y * stab_y1 * xy)
    cor_xyr <- (1-stab_x1*stab_y1-xy*yx)*cor_xy - stab_x1*yx*x - stab_y1*xy*y
    cor_xyr2 <- cor_xyr*1.8
    cor_xyr3 <- cor_xyr*.7

    ## Initialize and Calculate Matrices
    ## Lambda
    lambda <- matrix(0, nrow = 2 * nwaves, ncol = 2 * nwaves,
                     dimnames = list(c(paste0("x",1:nwaves),
                                       paste0("y",1:nwaves)),
                                     c(paste0("x",1:nwaves),
                                       paste0("y",1:nwaves))))

    ## lambda
    for (i in 1:(2*nwaves)) {
        lrow <- i
        lcol <- i
        lambda[lrow, lcol] <- 1
    }


    ## Theta
    theta <- matrix(0, nrow = 2 * nwaves, ncol = 2 * nwaves,
                    dimnames= list(c(paste0("x", 1:nwaves),
                                     paste0("y", 1:nwaves)),
                                   c(paste0("x", 1:nwaves),
                                     paste0("y", 1:nwaves))))
    theta[1:nwaves, 1:nwaves] <- diag(xr, nrow = nwaves)
    theta[(nwaves+1):(2*nwaves),(nwaves+1):(2*nwaves)] <- diag(yr, nrow = nwaves)


    ## Psi
    psi <- matrix(0, nrow = 2 * nwaves, ncol = 2 * nwaves,
                  dimnames = list(c(paste0("x",1:nwaves),
                                    paste0("y", 1:nwaves)),
                                  c(paste0("x",1:nwaves),
                                    paste0("y", 1:nwaves))))
diag(psi) <- c(x, wxr2, wxr3, wxr4, rep(wxr4, nwaves-4),
               y, wyr2, wyr3, wyr4, rep(wyr4, nwaves-4))
    psi[(nwaves+1), 1] <- cor_xy
    psi[1, (nwaves+1)] <- cor_xy
    psi[(nwaves+2), 2] <- cor_xyr
    psi[2, (nwaves+2)] <- cor_xyr
    psi[(nwaves+3), 3] <- cor_xyr2
    psi[2, (nwaves+3)] <- cor_xyr2
    for (i in 4:nwaves) {
        prow <- i + nwaves
        pcol <- i
        psi[prow, pcol] <- cor_xyr3
        psi[pcol, prow] <- cor_xyr3
    }

    ## Beta
    beta <- matrix(0, nrow = 2 * nwaves, ncol = 2 * nwaves,
                   dimnames = list(c(paste0("x",1:nwaves),
                                     paste0("y",1:nwaves)),
                                   c(paste0("x",1:nwaves),
                                     paste0("y",1:nwaves))))
    ## beta
    for (i in 1:(nwaves-1)) {
        ## x lag-1 stabilities
        xsrow <- i+1
        xscol <- i
        beta[xsrow, xscol] <- stab_x1
        ## y lag-1 stabilities
        ysrow <- i+1+nwaves
        yscol <- i+nwaves
        beta[ysrow, yscol] <- stab_y1
    }
    for (i in 1:(nwaves-2)) {
        ## x lag-1 stabilities
        xsrow <- i+2
        xscol <- i
        beta[xsrow, xscol] <- stab_x2
        ## y lag-1 stabilities
        ysrow <- i+2+nwaves
        yscol <- i+nwaves
        beta[ysrow, yscol] <- stab_y2
    }
    for (i in 1:(nwaves-3)) {
        ## x lag-3 stabilities
        xsrow <- i+3
        xscol <- i
        beta[xsrow, xscol] <- stab_x3
        ## y lag-3 stabilities
        ysrow <- i+3+nwaves
        yscol <- i+nwaves
        beta[ysrow, yscol] <- stab_y3
    }
    for (i in 1:(nwaves-1)) {
        ## y~x cross-lagged
        ycrow <- i+1+nwaves
        yccol <- i
        beta[ycrow, yccol] <- yx
        ## x~y cross-lagged
        xcrow <- i+1
        xccol <- i+(nwaves)
        beta[xcrow, xccol] <- xy
    }

    diag_length <- nwaves + nwaves
    ## Generate latent factor scores
    eta <- rmnorm(n, varcov = (solve(diag(diag_length)-beta) %*% psi %*% t(solve(diag(diag_length)-beta))))
    ## Generate residuals (all zero in ri-clpm)
    ifelse(xr==0,
           ex <- matrix(0, nrow = n, ncol = nwaves),
           ex <- rmnorm(n, varcov = theta[1:nwaves,1:nwaves]))
    ifelse(yr==0,
           ey <- matrix(0, nrow = n, ncol = nwaves),
           ey <- rmnorm(n, varcov = theta[(nwaves+1):(2*nwaves),(nwaves+1):(2*nwaves)]))
    e <- cbind(ex,ey)
    ## Compute observed scores
    obs <- tcrossprod(eta, lambda) + e
    ## Make it a dataframe
    data <- as.data.frame(obs)
}


data <- gen_clpm(n = 100000, nwaves=100, x = 1,y = 1, stab_x1 = .4, stab_x2 = .2, stab_x3 = 0, stab_y1 = .4, stab_y2 = .2, stab_y3 = 0, yx=.2, xy=.2, xr=0, yr=0, cor_xy = .5)



cov(data[,c(91:100,191:200)])

dataF <- data
data <- data[,c(91:100,191:200)]
names(data) <- paste0(c(rep("x",10),rep("y",10)),1:10)

m1 <- lavaan(clpm10_c, data)
summary(m1)

m2 <- lavaan(ri_clpm10_c, data)
summary(m2)

m3 <- lavaan(clpm3_c, data)
summary(m3)

m4 <- lavaan(ri_clpm3_c, data)
summary(m4)

m5 <- lavaan(clpm10_l2, data)
summary(m5)


data2 <- gen_starts(yx=0, xy=0, n=100000)

mb1 <- lavaan(ri_clpm10_c, data2)
summary(mb1)
mb2 <- lavaan(clpm10_l2, data2)
summary(mb2)





lag2 <- '
x10 ~ x9 + x8
x9 ~ x8 + x7
x8 ~ x7 + x6
x7 ~ x6 + x5
x6 ~ x5 + x4
x5 ~ x4 + x3
x4 ~ x3 + x2
x3 ~ x2 + x1
x2 ~ x1
x1 ~~ x1
x2 ~~ x2
x3 ~~ x3
x4 ~~ x4
x5 ~~ x5
x6 ~~ x6
x7 ~~ x7
x8 ~~ x8
x9 ~~ x9
x10 ~~ x10
'

testM <- lavaan(lag2, data)
summary(testM, rsquare=TRUE)


## x2

1 - (1*.5^2)

## x3

1 - (1*.5^2) - (1*.25^2) - 2*.25*1*.50*1*.5



################################################################################
## scratch
################################################################################




solve(diag(diag_length)-beta) %*% psi %*% t(solve(diag(diag_length)-beta))



varcovT <- matrix(c(1, .50, .50, .375, .3125, .25, .2031250, .1640625, .1328125, .1074219,
                    .50, 1, .50, .50, .375, .3125, .25, .2031250, .1640625, .1328125,
                    .50, .50, 1, .50, .50, .375, .3125, .25, .2031250, .1640625,
                    .375, .50, .50, 1, .50, .50, .375, .3125, .25, .2031250,
                    .3125, .375, .50, .50, 1, .50, .50, .375, .3125, .25,
                    .25, .3125, .375, .50, .50, 1, .50, .50, .375, .3125,
                    .2031250, .25, .3125, .375, .50, .50, 1, .50, .50, .375,
                    .1640625, .2031250, .25, .3125, .375, .50, .50, 1, .50, .50,
                    .1328125, .1640625, .2031250, .25, .3125, .375, .50, .50, 1, .50,
                    .1074219, .1328125, .1640625, .2031250, .25, .3125, .375, .50, .50, 1),
                  nrow=10, ncol=10, byrow = TRUE)

eta <- rmnorm(n, varcov = varcovT)
obs <- tcrossprod(eta, lambda[1:10, 1:10])

testM <- lavaan(lag2, obs)
summary(testM, rsquare=TRUE)


varx2 <- output[18,4] * output[17,4]^2 + output[19,4]
r2x2 <- output[18,4] * output[17,4]^2 / varx1
varx2
r2x2

varx3 <- varx2 * output[15,4]^2 +
    output[18,4] * output[16,4]^2 +
    2 * output[16,4] * output[18,4] * output[17,4] * varx2 * output[15,4] +
    output[19,4]
r2x3 <- varx2 * output[15,4]^2 +
    output[18,4] * output[16,4]^2 +
    2 * output[16,4] * output[18,4] * output[17,4] * varx2 * output[15,4] /
    varx3
varx3
r2x3

varx1 <- 1
b1 <- .5
b2 <- .25

resid <- (varx1 * ((1 + b2) * (1 - b1 - b2) * (1 + b1 - b2))) / (1 - b2)


vx <- 1
vy <- 1
bx1 <- .5
bx2 <- .2
by1 <- .5
by2 <- .2
yx <- .2
xy <- .2
xyr <- .5
n <- 100000

## Generate Initial Xs and Ys
data <- as.data.frame(rmnorm(n=n, varcov = matrix(c(vx, xyr, xyr, vy), nrow=2, ncol=2, byrow=TRUE)))
names(data) <- c("x1", "y1")
x2e <- (vx-(bx1^2*var(data$x1) + xy^2*var(data$y1) + 2*cor(data$x1, data$y1)*bx1*xy))
y2e <- (vy-(by1^2*var(data$y1) + yx^2*var(data$x1) + 2*cor(data$x1, data$y1)*by1*yx))
cor_xyr <- (1-bx1*by1-xy*yx)*xyr - bx1*yx*vx - by1*xy*vy
w2resid <- as.data.frame(rmnorm(n=n, varcov = matrix(c(x2e, cor_xyr, cor_xyr, y2e), nrow=2, ncol=2, byrow=TRUE)))
data$x2 <- bx1*data$x1 + xy*data$y1 + w2resid$V1
data$y2 <- by1*data$y1 + yx*data$x1 + w2resid$V2
X3temp <- bx1*data$x2 + bx2*data$x1 + xy*data$y2
Y3temp <- by1*data$y2 + by2*data$y1 + yx*data$x2
X3Var <- var(X3temp)
X3Resid <- var(data$x2) - X3Var
Y3Var <- var(Y3temp)
Y3Resid <- var(data$y2) - Y3Var
X3Y3Cov <- cov(X3temp,Y3temp)
newCov <- cov(data$x2, data$y2) - X3Y3Cov
w3resid <- as.data.frame(rmnorm(n=n,
                             varcov=matrix(c(X3Resid, newCov, newCov, Y3Resid), nrow = 2, ncol = 2, byrow = TRUE)))
data$x3 <- bx1*data$x2 + bx2*data$x1 + xy*data$y2 + w3resid$V1
data$y3 <- by1*data$y2 + by2*data$y1 + yx*data$x2 + w3resid$V2
cov(data)




x3e <- (vx- (bx1^2*var(data$x2) +
        xy^2*var(data$y2) +
        bx2^2*var(data$x1) +
        2*bx2*var(data$x1)*bx1*var(data$x1)*bx1 +
        2*bx2*cor(data$y1, data$x1)*by1*xy +
        2*bx2*cor(data$x1, data$y1)*xy*bx1 +
        2*bx2*var(data$x1)*yx*var(data$y2)*xy +
        2*bx1*var(data$x2)*bx1*cor(data$x1, data$y1)*by1*var(data$y2)*xy +
        2*bx1*var(data$x2)*bx1*cor(data$x1, data$y1)*xy*var(data$x2)*bx1 +
        2*bx1*var(data$x2)*bx1*var(data$x2)*yx*var(data$y2)*xy +
        2*bx1*cor_xyr*xy))
y3e <- (vy- (by1^2*var(data$y2) +
        yx^2*var(data$x2) +
        by2^2*var(data$y1) +
        2*by2*var(data$y1)*by1*var(data$y1)*by1 +
        2*by2*cor(data$x1, data$y1)*bx1*yx +
        2*by2*cor(data$y1, data$x1)*yx*by1 +
        2*by2*var(data$y1)*xy*var(data$x2)*yx +
        2*by1*var(data$y2)*by1*cor(data$y1, data$x1)*bx1*var(data$x2)*yx +
        2*by1*var(data$y2)*by1*cor(data$y1, data$x1)*yx*var(data$y2)*by1 +
        2*by1*var(data$y2)*by1*var(data$y2)*xy*var(data$x2)*yx +
        2*by1*cor_xyr*yx))
cor_xyr3 <- cor_xyr -
    bx2*vx*bx1*vx*yx -
    bx2*vx*yx*vy*by1 -
    by2*vy*by1*vy*xy -
    by2*vy*xy*vx*bx1
#cor_xyr3 <- 0
w3resid <- as.data.frame(rmnorm(n=n, varcov = matrix(c(x3e, cor_xyr3, cor_xyr3, y3e), nrow=2, ncol=2, byrow=TRUE)))




### Figure out variance

X3 <- bx1*data$x2 + bx2*data$x1 + xy*data$y2

var(X3)
bx1^2*var(data$x2) 
xy^2*var(data$y2) 
bx2^2*var(data$x1) 
2*bx2*var(data$x1)*bx1*var(data$x1)*bx1 
2*bx2*cor(data$y1, data$x1)*by1*xy 
2*bx2*cor(data$x1, data$y1)*xy*bx1 
2*bx2*var(data$x1)*yx*var(data$y2)*xy 
2*bx1*var(data$x2)*bx1*cor(data$x1, data$y1)*by1*var(data$y2)*xy 
2*bx1*var(data$x2)*bx1*cor(data$x1, data$y1)*xy*var(data$x2)*bx1 
2*bx1*var(data$x2)*bx1*var(data$x2)*yx*var(data$y2)*xy 
2*bx1*cor_xyr*xy


bx1^2*var(data$x2) +
    xy^2*var(data$y2) +
        bx2^2*var(data$x1) +
        2*bx2*var(data$x1)*bx1*var(data$x1)*bx1 +
        2*bx2*cor(data$y1, data$x1)*by1*xy +
        2*bx2*cor(data$x1, data$y1)*xy*bx1 +
        2*bx2*var(data$x1)*yx*var(data$y2)*xy +
        2*bx1*var(data$x2)*bx1*cor(data$x1, data$y1)*by1*var(data$y2)*xy +
        2*bx1*var(data$x2)*bx1*cor(data$x1, data$y1)*xy*var(data$x2)*bx1 +
        2*bx1*var(data$x2)*bx1*var(data$x2)*yx*var(data$y2)*xy +
        2*bx1*cor_xyr*xy







## Predicted Variance
b2^2*var(data$V1)
b1^2*var(data$V3)
2 * b2 * b1 * b1
2 * b2 * xr * b3 * b1

b2^2*var(data$V1) +
    b1^2*var(data$V3) +
        2 * b2 * b1 * b1 +
        2 * b2 * xr * b3 * b1

    



v1 <- rnorm(10000, mean=0, sd=1)
v2 <- rnorm(10000, mean=0, sd=sqrt(.75))

.5^2*var(data$V1) + .5^2*var(data$V2) + 2*cor(data$V1, data$V2)*.5*.5



################################################################################
## New approach
################################################################################

gen_clpm <- function(n=10000,
                     nwaves=10,
                     stab_x1=.5,
                     stab_x2=.25,
                     stab_y1=.5,
                     stab_y2=.25,
                     stab_y3=0,
                     yx=0,
                     yx2=0,
                     xy=0,
                     xy2=0,
                     x=1,
                     y=1,
                     cor_xy=.5,
                     xr=0,
                     yr=0
                     ) {


    ## Temporary Simulation to ensure stationarity
    ## Generate Initial Xs and Ys
    vx <- x
    vy <- y
    xyr <- cor_xy
    bx1 <- stab_x1
    bx2 <- stab_x2
    by1 <- stab_y1
    by2 <- stab_y2
    ##
    ## Create dataframe for initial variables
    data <- as.data.frame(rmnorm(n=n, varcov = matrix(c(vx, xyr, xyr, vy), nrow=2, ncol=2, byrow=TRUE)))
    names(data) <- c("x1", "y1")
    ## Figure out Residual Variance for Time 2:
    x2e <- (vx-(bx1^2*var(data$x1) + xy^2*var(data$y1) + 2*cor(data$x1, data$y1)*bx1*xy))
    wxr2 <- x2e
    y2e <- (vy-(by1^2*var(data$y1) + yx^2*var(data$x1) + 2*cor(data$x1, data$y1)*by1*yx))
    wyr2 <- y2e
    ## Figure out Residual Covariance for Time 2:
    cor_xyr <- (1-bx1*by1-xy*yx)*xyr - bx1*yx*vx - by1*xy*vy
    ## Calculate residuals
    w2resid <- as.data.frame(rmnorm(n=n, varcov = matrix(c(x2e, cor_xyr, cor_xyr, y2e), nrow=2, ncol=2, byrow=TRUE)))
    ## Create new Wave-2 Variables
    data$x2 <- bx1*data$x1 + xy*data$y1 + w2resid$V1
    data$y2 <- by1*data$y1 + yx*data$x1 + w2resid$V2
    ## Create temporary Wave 3 variables used to figure out variances and covariances
    X3temp <- bx1*data$x2 + bx2*data$x1 + xy*data$y2 + xy2*data$y1
    Y3temp <- by1*data$y2 + by2*data$y1 + yx*data$x2 + yx2*data$x1
    X3Var <- var(X3temp)
    wxr3 <- var(data$x2) - X3Var ## Residual variance for Wave 3 and beyond
    Y3Var <- var(Y3temp)
    wyr3 <- var(data$y2) - Y3Var ## Residual variance for Wave 3 and beyond
    X3Y3Cov <- cov(X3temp,Y3temp)
    cor_xyr2 <- cov(data$x2, data$y2) - X3Y3Cov ## Residual covariance for Wave 3 and beyond
    w3resid <- as.data.frame(rmnorm(n=n,
                                    varcov=matrix(c(wxr3, cor_xyr2, cor_xyr2, wyr3), nrow = 2, ncol = 2, byrow = TRUE)))
    data$x3 <- bx1*data$x2 + bx2*data$x1 + xy*data$y2 + xy2*data$y1 + w3resid$V1
    data$y3 <- by1*data$y2 + by2*data$y1 + yx*data$x2 + yx2*data$x1 + w3resid$V2
    X4temp <- bx1*data$x3 + bx2*data$x2 + xy*data$y3 + xy2*data$y2
    Y4temp <- by1*data$y3 + by2*data$y2 + yx*data$x3 + yx2*data$x2
    X4Var <- var(X4temp)
    wxr4 <- var(data$x3) - X4Var
    Y4Var <- var(Y4temp)
    wyr4 <- var(data$x3) - Y4Var
    X4Y4Cov <- cov(X4temp,Y4temp)
    cor_xyr3 <- cov(data$x3, data$y3) - X4Y4Cov ## Residual covariance for Wave 3 and beyond
    w4resid <- as.data.frame(rmnorm(n=n,
                                    varcov=matrix(c(wxr4, cor_xyr3, cor_xyr3, wyr4), nrow = 2, ncol = 2, byrow = TRUE)))
    data$x4 <- bx1*data$x3 + bx2*data$x2 + xy*data$y3 + xy2*data$y2 + w4resid$V1
    data$y4 <- by1*data$y3 + by2*data$y2 + yx*data$x3 + yx2*data$x2 + w4resid$V2





    var(X3temp)
    var(Y3temp)
    var(X4temp)
    var(Y4temp)

    ## Initialize and Calculate Matrices
    ## Lambda
    lambda <- matrix(0, nrow = 2 * nwaves, ncol = 2 * nwaves,
                     dimnames = list(c(paste0("x",1:nwaves),
                                       paste0("y",1:nwaves)),
                                     c(paste0("x",1:nwaves),
                                       paste0("y",1:nwaves))))

    ## lambda
    for (i in 1:(2*nwaves)) {
        lrow <- i
        lcol <- i
        lambda[lrow, lcol] <- 1
    }


    ## Theta
    theta <- matrix(0, nrow = 2 * nwaves, ncol = 2 * nwaves,
                    dimnames= list(c(paste0("x", 1:nwaves),
                                     paste0("y", 1:nwaves)),
                                   c(paste0("x", 1:nwaves),
                                     paste0("y", 1:nwaves))))
    theta[1:nwaves, 1:nwaves] <- diag(xr, nrow = nwaves)
    theta[(nwaves+1):(2*nwaves),(nwaves+1):(2*nwaves)] <- diag(yr, nrow = nwaves)


    ## Psi
    psi <- matrix(0, nrow = 2 * nwaves, ncol = 2 * nwaves,
                  dimnames = list(c(paste0("x",1:nwaves),
                                    paste0("y", 1:nwaves)),
                                  c(paste0("x",1:nwaves),
                                    paste0("y", 1:nwaves))))
    diag(psi) <- c(x, wxr2, wxr3, rep(wxr3, nwaves-3),
                   y, wyr2, wyr3, rep(wyr3, nwaves-3))
    psi[(nwaves+1), 1] <- cor_xy
    psi[1, (nwaves+1)] <- cor_xy
    psi[(nwaves+2), 2] <- cor_xyr
    psi[2, (nwaves+2)] <- cor_xyr
    for (i in 3:nwaves) {
        prow <- i + nwaves
        pcol <- i
        psi[prow, pcol] <- cor_xyr2
        psi[pcol, prow] <- cor_xyr2
    }

    ## Beta
    beta <- matrix(0, nrow = 2 * nwaves, ncol = 2 * nwaves,
                   dimnames = list(c(paste0("x",1:nwaves),
                                     paste0("y",1:nwaves)),
                                   c(paste0("x",1:nwaves),
                                     paste0("y",1:nwaves))))
    ## beta
    for (i in 1:(nwaves-1)) {
        ## x lag-1 stabilities
        xsrow <- i+1
        xscol <- i
        beta[xsrow, xscol] <- stab_x1
        ## y lag-1 stabilities
        ysrow <- i+1+nwaves
        yscol <- i+nwaves
        beta[ysrow, yscol] <- stab_y1
    }
    for (i in 1:(nwaves-2)) {
        ## x lag-1 stabilities
        xsrow <- i+2
        xscol <- i
        beta[xsrow, xscol] <- stab_x2
        ## y lag-1 stabilities
        ysrow <- i+2+nwaves
        yscol <- i+nwaves
        beta[ysrow, yscol] <- stab_y2
    }
    for (i in 1:(nwaves-1)) {
        ## y~x cross-lagged
        ycrow <- i+1+nwaves
        yccol <- i
        beta[ycrow, yccol] <- yx
        ## x~y cross-lagged
        xcrow <- i+1
        xccol <- i+(nwaves)
        beta[xcrow, xccol] <- xy
    }
    for (i in 1:(nwaves-2)) {
        ## y~x cross-lagged
        ycrow <- i+2+nwaves
        yccol <- i
        beta[ycrow, yccol] <- yx2
        ## x~y cross-lagged
        xcrow <- i+2
        xccol <- i+(nwaves)
        beta[xcrow, xccol] <- xy2
    }


    diag_length <- nwaves + nwaves
    ## Generate latent factor scores
    eta <- rmnorm(n, varcov = (solve(diag(diag_length)-beta) %*% psi %*% t(solve(diag(diag_length)-beta))))
    ## Generate residuals (all zero in ri-clpm)
    ifelse(xr==0,
           ex <- matrix(0, nrow = n, ncol = nwaves),
           ex <- rmnorm(n, varcov = theta[1:nwaves,1:nwaves]))
    ifelse(yr==0,
           ey <- matrix(0, nrow = n, ncol = nwaves),
           ey <- rmnorm(n, varcov = theta[(nwaves+1):(2*nwaves),(nwaves+1):(2*nwaves)]))
    e <- cbind(ex,ey)
    ## Compute observed scores
    obs <- tcrossprod(eta, lambda) + e
    ## Make it a dataframe
    data <- as.data.frame(obs)
}
