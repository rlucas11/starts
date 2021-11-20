run_sim <- function() {
    data <- gen_starts(n=n,      # N to generate
                       nwaves=nwaves,   # Number of waves
                       ri_x=ri_x,     # Random intercept variance for X
                       ri_y=ri_y,     # Random intercept variance for Y
                       cor_i=cor_i,   # Correlation between intercepts
                       x=x,        # AR variance for X
                       y=y,        # AR variance for Y
                       stab_x=stab_x,  # Stability of X
                       stab_y=stab_y,  # Stability of Y
                       yx=yx,      # Cross lag (Y on X)
                       xy=xy,      # Cross lag (X on Y)
                       cor_xy=cor_xy,  # Correlation between X and Y
                       cor_xyr=cor_xyr, # Correlation between X and Y residuals
                       xr=xr,       # Measurement error for X
                       yr=yr        # Measurement error for Y
                       )
    ##
    fit_clpm2 <- lavaan(clpm2, data = data)
    fit_clpm5 <- lavaan(clpm5, data = data)
    fit_clpm10 <- lavaan(clpm10, data = data)
    fit_riclpm <- lavaan(ri_clpm10, data = data)
    output <- list(clpm2=parameterEstimates(fit_clpm2),
                   clpm5=parameterEstimates(fit_clpm5),
                   clpm10=parameterEstimates(fit_clpm10),
                   riclpm=parameterEstimates(fit_riclpm))
    output
}

run_sim2 <- function() {
    data <- gen_starts(n=n,      # N to generate
                       nwaves=nwaves,   # Number of waves
                       ri_x=ri_x,     # Random intercept variance for X
                       ri_y=ri_y,     # Random intercept variance for Y
                       cor_i=cor_i,   # Correlation between intercepts
                       x=x,        # AR variance for X
                       y=y,        # AR variance for Y
                       stab_x=stab_x,  # Stability of X
                       stab_y=stab_y,  # Stability of Y
                       yx=yx,      # Cross lag (Y on X)
                       xy=xy,      # Cross lag (X on Y)
                       cor_xy=cor_xy,  # Correlation between X and Y
                       cor_xyr=cor_xyr, # Correlation between X and Y residuals
                       xr=xr,       # Measurement error for X
                       yr=yr        # Measurement error for Y
                       )
    ##
    fit_clpm2 <- lavaan(clpm2, data = data)
    clpm2=parameterEstimates(fit_clpm2)
    clpm2
}

run_sim5 <- function() {
    data <- gen_starts(n=n,      # N to generate
                       nwaves=nwaves,   # Number of waves
                       ri_x=ri_x,     # Random intercept variance for X
                       ri_y=ri_y,     # Random intercept variance for Y
                       cor_i=cor_i,   # Correlation between intercepts
                       x=x,        # AR variance for X
                       y=y,        # AR variance for Y
                       stab_x=stab_x,  # Stability of X
                       stab_y=stab_y,  # Stability of Y
                       yx=yx,      # Cross lag (Y on X)
                       xy=xy,      # Cross lag (X on Y)
                       cor_xy=cor_xy,  # Correlation between X and Y
                       cor_xyr=cor_xyr, # Correlation between X and Y residuals
                       xr=xr,       # Measurement error for X
                       yr=yr        # Measurement error for Y
                       )
    ##
    fit_clpm5 <- lavaan(clpm5, data = data)
    clpm5=parameterEstimates(fit_clpm5)
    clpm5
}

run_sim10 <- function() {
    data <- gen_starts(n=n,      # N to generate
                       nwaves=nwaves,   # Number of waves
                       ri_x=ri_x,     # Random intercept variance for X
                       ri_y=ri_y,     # Random intercept variance for Y
                       cor_i=cor_i,   # Correlation between intercepts
                       x=x,        # AR variance for X
                       y=y,        # AR variance for Y
                       stab_x=stab_x,  # Stability of X
                       stab_y=stab_y,  # Stability of Y
                       yx=yx,      # Cross lag (Y on X)
                       xy=xy,      # Cross lag (X on Y)
                       cor_xy=cor_xy,  # Correlation between X and Y
                       cor_xyr=cor_xyr, # Correlation between X and Y residuals
                       xr=xr,       # Measurement error for X
                       yr=yr        # Measurement error for Y
                       )
    ##
    fit_clpm10 <- lavaan(clpm10, data = data)
    clpm10=parameterEstimates(fit_clpm10)
    clpm10
}

run_simri <- function() {
    data <- gen_starts(n=n,      # N to generate
                       nwaves=nwaves,   # Number of waves
                       ri_x=ri_x,     # Random intercept variance for X
                       ri_y=ri_y,     # Random intercept variance for Y
                       cor_i=cor_i,   # Correlation between intercepts
                       x=x,        # AR variance for X
                       y=y,        # AR variance for Y
                       stab_x=stab_x,  # Stability of X
                       stab_y=stab_y,  # Stability of Y
                       yx=yx,      # Cross lag (Y on X)
                       xy=xy,      # Cross lag (X on Y)
                       cor_xy=cor_xy,  # Correlation between X and Y
                       cor_xyr=cor_xyr, # Correlation between X and Y residuals
                       xr=xr,       # Measurement error for X
                       yr=yr        # Measurement error for Y
                       )
    ##
    fit_riclpm <- lavaan(ri_clpm10, data = data)
    clpmri=parameterEstimates(fit_riclpm)
    clpmri
}
