library(mnormt)

## Testing
n=10000
nwaves=10
x=1
y=1
stab_x1=.4
stab_x2=.4
stab_y1=.4
stab_y2=.4
yx=.1
xy=.0
yx2=.1
xy2=.0
cor_xy=.5
xr=0
yr=0
drop=TRUE
stationary=TRUE
burnin=10


gen_clpm <- function(n=10000,
                     nwaves=10,
                     stab_x1=.4,
                     stab_x2=.2,
                     stab_y1=.4,
                     stab_y2=.2,
                     yx=.1,
                     yx2=.1,
                     xy=.1,
                     xy2=.1,
                     x=1,
                     y=1,
                     cor_xy=.5,
                     xr=0,
                     yr=0,
                     xresid = NA, ## Specify residual variance if not stationary
                     yresid = NA, ## Specify residual variance if not stationary
                     drop=TRUE,  ## Drop burnin waves to get stationary pattern
                     stationary=TRUE, ## 
                     burnin=10
                     ) {


    ## Rename variables (because I'm reusing old code)
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

    ## Generate Data for Additional Waves

    for (w in 3:(nwaves+burnin)) {

        x1var <- paste0("x", (w-2))
        x2var <- paste0("x", (w-1))
        y1var <- paste0("y", (w-2))
        y2var <- paste0("y", (w-1))
        ## Create predicted scores to figure out residual variance and covariances
        predicted <- mutate(data, xpred=bx1*.data[[x2var]] +
                                      bx2*.data[[x1var]] +
                                      xy*.data[[y2var]] +
                                      xy2*.data[[y1var]],
                            ypred=by1*.data[[y2var]] +
                                      by2*.data[[y1var]] +
                                      yx*.data[[x2var]] +
                                      yx2*.data[[x1var]],
                            .keep="none")
        predicted_cov <- cov(predicted)
        new_x_var <- x-predicted_cov[[1,1]]
        new_y_var <- y-predicted_cov[[2,2]]
        new_xy_cov <- xyr-predicted_cov[[1,2]]
        ## Generate new x and y residual variables
        new_resid <- as.data.frame(rmnorm(n=n, varcov = matrix(c(new_x_var, new_xy_cov, new_xy_cov, new_y_var),
                                                               nrow=2, ncol=2, byrow=TRUE)))
        names(new_resid) <- c("xresid", "yresid")
        new_vars <- cbind(predicted, new_resid)
        xname <- paste0("x",w)
        yname <- paste0("y",w)
        new_vars <- mutate(new_vars, "{xname}" := xpred + xresid,
                           "{yname}" := ypred + yresid)
        data <- cbind(data, new_vars[,5:6])
        
    }

        
    measurement_error <- matrix(0, nrow=n, ncol=(nwaves+burnin)*2)
    measurement_error_x <- matrix(0, nrow=n, ncol=(nwaves+burnin))
    measurement_error_y <- matrix(0, nrow=n, ncol=(nwaves+burnin))
    if(xr>0) { measurement_error_cov_x <- matrix(0, nrow=(nwaves+burnin), ncol=(nwaves+burnin))
        diag(measurement_error_cov_x) <- rep(xr, (nwaves+burnin))
        measurement_error_x <- rmnorm(n=n, varcov = measurement_error_cov_x)
    }
    if(yr>0) { measurement_error_cov_y <- matrix(0, nrow=(nwaves+burnin), ncol=(nwaves+burnin))
        diag(measurement_error_cov_y) <- rep(yr, (nwaves+burnin))
        measurement_error_y <- rmnorm(n=n, varcov = measurement_error_cov_y)
    }
    measurement_error[,seq(1, (nwaves+burnin)*2, by=2)] <- measurement_error_x
    measurement_error[,seq(2, (nwaves+burnin)*2, by=2)] <- measurement_error_y
    data <- data+measurement_error

    ifelse(drop==TRUE,
    {data <- data[,-c(1:(2*burnin))]
        names(data) <- paste0(rep(c("x","y"), nwaves), rep(1:nwaves, each=2))},
    {data <- data
    names(data) <- c(paste0(rep(c("x","y"), burnin), rep(-(burnin-1):0, each=2)),
                     paste0(rep(c("x","y"), nwaves), rep(1:nwaves, each=2)))})

    return(data)
}

