summarizeR <- function(corMat, nvars=1) {

    averageRs <- matrix(nrow=(nrow(corMat)/nvars-1),
                        ncol=nvars)

    for (k in 1:nvars) {
        for (i in 1:((nrow(corMat)/nvars)-1)) {
            sumR <- 0
            nValid <- 0
            for (j in seq(1, (nrow(corMat)-nvars*i), by=nvars)) {
                if(!is.na(corMat[(j+(i*nvars)+(k-1)), j+(k-1)])) {
                    sumR <- sumR + corMat[(j+(i*nvars)+(k-1)), j+(k-1)]
                    nValid <- nValid + 1
                }
                
            }
            averageRs[i,k] <- sumR/(nValid)
        }
    }
    return(averageRs)
}


plotCors <- function(cors) {
    cors <- as.data.frame(cors)
    cors <- cors %>%
        mutate(lag=row_number()) %>%
        pivot_longer(cols=starts_with("V"))
    minCor <- min(cors$value)

    ggplot(aes(x=lag, y=value, group=name, color=name), data=cors) +
        geom_smooth(se=FALSE) +
        ylim(min(minCor,0), 1)
}


arCor <- function(waves, a) {
  cors <- matrix(nrow = (waves - 1), ncol = 1)
  for (i in 1:(waves - 1)) {
    cors[i, 1] <- a^i
  }
  return(cors)
}


