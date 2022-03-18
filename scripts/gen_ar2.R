## Testing
n <- 10000
nwaves <- 20
stab_x1 <- .5
stab_x2 <- .25
xr <- 0

gen_ar2 <- function(n=10000,
                     nwaves=20,
                     stab_x1=.5,
                     stab_x2=.25,
                     xr=0
                     ) {
    x <- 1
    ## Stationarity Constraints (not accurate yet)
    wxr2 <- (1-stab_x1^2)*x 
    wxr3 <- (x * ((1 + stab_x2) * (1 - stab_x1 - stab_x2) * (1 + stab_x1 - stab_x2))) / (1 - stab_x2) 

    ## Initialize and Calculate Matrices
    ## Lambda
    lambda <- matrix(0, nrow = nwaves, ncol = nwaves,
                     dimnames = list(c(paste0("x",1:nwaves)),
                                     c(paste0("x",1:nwaves))))

    ## lambda
    for (i in 1:(nwaves)) {
        lrow <- i
        lcol <- i
        lambda[lrow, lcol] <- 1
    }


    ## Theta
    theta <- matrix(0, nrow = nwaves, ncol = nwaves,
                    dimnames= list(c(paste0("x", 1:nwaves)),
                                   c(paste0("x", 1:nwaves))))
    theta[1:nwaves, 1:nwaves] <- diag(xr, nrow = nwaves)

    ## Psi
    psi <- matrix(0, nrow = nwaves, nwaves,
                  dimnames = list(c(paste0("x",1:nwaves)),
                                  c(paste0("x",1:nwaves))))
    diag(psi) <- c(x, wxr2, wxr3, rep(wxr3, nwaves-3))

    ## Beta
    beta <- matrix(0, nrow = nwaves, ncol = nwaves,
                   dimnames = list(c(paste0("x",1:nwaves)),
                                   c(paste0("x",1:nwaves))))
    ## beta
    for (i in 1:(nwaves-1)) {
        ## x lag-1 stabilities
        xsrow <- i+1
        xscol <- i
        beta[xsrow, xscol] <- stab_x1
    }
    for (i in 1:(nwaves-2)) {
        ## x lag-1 stabilities
        xsrow <- i+2
        xscol <- i
        beta[xsrow, xscol] <- stab_x2
    }

    diag_length <- nwaves
    ## Generate latent factor scores
    eta <- rmnorm(n, varcov = (solve(diag(diag_length)-beta) %*% psi %*% t(solve(diag(diag_length)-beta))))
    ## Generate residuals (all zero in ri-clpm)
    ifelse(xr==0,
           ex <- matrix(0, nrow = n, ncol = nwaves),
           ex <- rmnorm(n, varcov = theta[1:nwaves,1:nwaves]))
    e <- cbind(ex)
    ## Compute observed scores
    obs <- tcrossprod(eta, lambda) + e
    ## Make it a dataframe
    data <- as.data.frame(obs)
}


avg_cor <- function(cor_mat) {
    waves <- nrow(cor_mat)
    results <- matrix(nrow=waves-1)
    sum_cor <- 0
    for (i in 1:(waves-1)) {
        for (j in 1:(waves-i)) {
            sum_cor <- sum_cor + cor_mat[i+j, j]
        }
        results[i] <- sum_cor/(waves - i)
        sum_cor <- 0
    }
    results
}


arData <- gen_ar2(stab_x1=.6, stab_x2=.3)
arDataS <- arData[,5:20]
arAvgCor <- avg_cor(cor(arDataS))
arAvgCor


stData <- gen_starts(n=10000, nwaves = 15, ri_x = 1, cor_i=0, x=1, stab_x=.5, yx=0, xy=0)
stDataS <- stData[,1:15]
stAvgCor <- avg_cor(cor(stDataS))
