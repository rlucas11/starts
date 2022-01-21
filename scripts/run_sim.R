## Function to simulate clpm data
run_sim_clpm <- function(waves=10,
                         studyN,      
                         ri_x,     
                         ri_y,     
                         cor_i,   
                         x,        
                         y,        
                         stab_x, 
                         stab_y, 
                         yx,      
                         xy,      
                         cor_xy, 
                         reliability_x,       
                         reliability_y) {
    sx <- sum(ri_x, x)/reliability_x - sum(ri_x, x)
    sy <- sum(ri_y, y)/reliability_y - sum(ri_y, y)
    data <- gen_starts(n=studyN,    
                       nwaves=waves, 
                       ri_x=ri_x,    
                       ri_y=ri_y,    
                       cor_i=cor_i,  
                       x=x,       
                       y=y,       
                       stab_x=stab_x, 
                       stab_y=stab_y, 
                       yx=yx,     
                       xy=xy,     
                       cor_xy=cor_xy, 
                       xr=sx,      
                       yr=sy       
                       )
    ##
    if (waves==10) {
        clpm <- clpm10_c
    } else if (waves==5) {
        clpm <- clpm5_c
    } else if (waves==3) {
        clpm <- clpm3_c
    } else if (waves==2) {
        clpm <- clpm2_c
    } else {
        stop("No model defined for that many waves")
    }
    fit_clpm <- lavaan(clpm, data = data)
    output <- c(parameterEstimates(fit_clpm)[3,5],
                parameterEstimates(fit_clpm)[3,8],
                parameterEstimates(fit_clpm)[4,5],
                parameterEstimates(fit_clpm)[4,8])
    output
}

run_sim_riclpm <- function(waves=10,
                         studyN=10000,      
                         ri_x=1,     
                         ri_y=1,     
                         cor_i=.5,   
                         x=1,        
                         y=1,        
                         stab_x=.5, 
                         stab_y=.5, 
                         yx=0,      
                         xy=0,      
                         cor_xy=.5, 
                         reliability_x=.8,       
                         reliability_y=.8) {
    sx <- sum(ri_x, x)/reliability_x - sum(ri_x, x)
    sy <- sum(ri_y, y)/reliability_y - sum(ri_y, y)
    data <- gen_starts(n=studyN,    
                       nwaves=waves, 
                       ri_x=ri_x,    
                       ri_y=ri_y,    
                       cor_i=cor_i,  
                       x=x,       
                       y=y,       
                       stab_x=stab_x, 
                       stab_y=stab_y, 
                       yx=yx,     
                       xy=xy,     
                       cor_xy=cor_xy, 
                       xr=sx,      
                       yr=sy       
                       )
    ##
    fit_riclpm <- lavaan(ri_clpm10_c, data = data)
    output <- c(parameterEstimates(fit_riclpm)[59,5],
                parameterEstimates(fit_riclpm)[59,8],
                parameterEstimates(fit_riclpm)[68,5],
                parameterEstimates(fit_riclpm)[68,8])
    output
}


