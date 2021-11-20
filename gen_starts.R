################################################################################
## Function to generate STARTS/CLPM/RI-CLPM Data
################################################################################

library(mnormt)

gen_starts <- function(n=500,      # N to generate
                       nwaves=5,   # Number of waves
                       ri_x=1,     # Random intercept variance for X
                       ri_y=1,     # Random intercept variance for Y
                       cor_i=.5,   # Correlation between intercepts (as correlation)
                       x=1,        # AR variance for X
                       y=1,        # AR variance for Y
                       stab_x=.5,  # Stability of X
                       stab_y=.5,  # Stability of Y
                       yx=.4,      # Cross lag (Y on X)
                       xy=.2,      # Cross lag (X on Y)
                       cor_xy=.5,  # Correlation between X and Y (as correlation)
                       cor_xyr=.2, # Correlation between X and Y residuals  (as correlation)
                       xr=0,       # Measurement error for X
                       yr=0        # Measurement error for Y
                       ) {
    ## Stationarity Constraints
    wxr <- 1 - ((x*stab_x^2 + y*xy^2)/(1-cor_xy^2)) #x residual
    wyr <- 1 - ((y*stab_y^2 + x*yx^2)/(1-cor_xy^2)) #y residual
    ## Transform correlations into covariances for matrices
    cor_i <- cor_i * (sqrt(ri_x) * sqrt(ri_y))
    cor_xyr <- cor_xyr * (sqrt(wxr) * sqrt(wyr))
    ##
    ## Initialize Matrices
    lambda <- matrix(0, nrow = 2 * nwaves, ncol = 2 + 2 * nwaves,
                     dimnames = list(c(paste0("x",1:nwaves),
                                       paste0("y",1:nwaves)),
                                     c("ri_x", "ri_y", paste0("x",1:nwaves),
                                       paste0("y",1:nwaves))))
    theta <- matrix(0, nrow = 2 * nwaves, ncol = 2 * nwaves,
                    dimnames= list(c(paste0("x", 1:nwaves),
                                    paste0("y", 1:nwaves)),
                                    c(paste0("x", 1:nwaves),
                                    paste0("y", 1:nwaves))))
    psi <- matrix(0, nrow = 2 + 2 * nwaves, ncol = 2 + 2 * nwaves,
                  dimnames = list(c("ri_x", "ri_y", paste0("x",1:nwaves),
                                    paste0("y", 1:nwaves)),
                                  c("ri_x", "ri_y", paste0("x",1:nwaves),
                                    paste0("y", 1:nwaves))))
    beta <- matrix(0, nrow = 2 + 2 * nwaves, ncol = 2 + 2 * nwaves,
                   dimnames = list(c("ri_x", "ri_y", paste0("x",1:nwaves),
                                     paste0("y",1:nwaves)),
                                   c("ri_x", "ri_y", paste0("x",1:nwaves),
                                     paste0("y",1:nwaves))))
    ##
    ## Fill in Matrices
    ## lambda
    lambda[1:nwaves, 1] <- 1 ## X loadings
    lambda[(nwaves+1):(2*nwaves), 2] <- 1  ## Y loadings
    for (i in 1:(2*nwaves)) {
        lrow <- i
        lcol <- i + 2
        lambda[lrow, lcol] <- 1
    }
    ## theta
    theta[1:nwaves, 1:nwaves] <- diag(xr, nrow = nwaves)
    theta[(nwaves+1):(2*nwaves),(nwaves+1):(2*nwaves)] <- diag(yr, nrow = nwaves)
    ## psi
    psi[1:2,1:2] <- c(ri_x, cor_i, cor_i, ri_y)
    diag(psi)[3:(2*nwaves+2)] <- c(x, rep(wxr, nwaves-1), y, rep(wyr, nwaves-1))
    psi[(nwaves+3), 3] <- cor_xy
    psi[3, (nwaves+3)] <- cor_xy
    for (i in 2:nwaves) {
        prow <- i + nwaves + 2
        pcol <- i + 2
        psi[prow, pcol] <- cor_xyr
        psi[pcol, prow] <- cor_xyr
    }
    ## beta
    for (i in 1:(nwaves-1)) {
        ## x stabilities
        xsrow <- i+3
        xscol <- i+2
        beta[xsrow, xscol] <- stab_x
        ## y stabilities
        ysrow <- i+3+nwaves
        yscol <- i+2+nwaves
        beta[ysrow, yscol] <- stab_y
    }
    for (i in 1:(nwaves-1)) {
        ## y~x cross-lagged
        ycrow <- i+3+nwaves
        yccol <- i+2
        beta[ycrow, yccol] <- yx
        ## x~y cross-lagged
        xcrow <- i+3
        xccol <- i+2+(nwaves)
        beta[xcrow, xccol] <- xy
    }
    ## Remove rows from matrices before generating data if no variance
    ifelse(ri_x==0 & ri_y==0,
           psi <- psi[-c(1:2),-c(1:2)],
    ifelse(ri_x==0 & ri_y>0,
           psi <- psi[-1,-1],
    ifelse(ri_x>0 & ri_y==0,
           psi <- psi[-2,-2],
           psi <- psi)))
    ifelse(ri_x==0 & ri_y==0,
           beta <- beta[-c(1:2),-c(1:2)],
    ifelse(ri_x==0 & ri_y>0,
           beta <- beta[-1,-1],
    ifelse(ri_x>0 & ri_y==0,
           beta <- beta[-2,-2],
           beta <- beta)))
    ifelse(ri_x==0 & ri_y==0,
           lambda <- lambda[,-c(1:2)],
    ifelse(ri_x==0 & ri_y>0,
           lambda <- lambda[,-1],
    ifelse(ri_x>0 & ri_y==0,
           lambda <- lambda[,-2],
           lambda <- lambda)))
    diag_length <- (2*nwaves) + sum(ri_x>0, ri_y>0) ## Dimensions of identity matrix
    ## Generate latent factor scores
    ##eta <- rmnorm(500, mean = 0, varcov = psi)
    eta <- rmnorm(n, varcov = (solve(diag(diag_length)-beta) %*%
                               psi %*% t(solve(diag(diag_length)-beta))))
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
