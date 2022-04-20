summarizeR <- function(corMat, nvars=1) {

    averageRs <- matrix(nrow=(nrow(corMat)/nvars-1),
                        ncol=nvars)

    for (k in 1:nvars) {
        for (i in 1:((nrow(corMat)/nvars)-1)) {
            sumR <- 0
            for (j in seq(1, (nrow(corMat)-nvars*i), by=nvars)) {
                sumR <- sumR + corMat[(j+(i*nvars)+(k-1)), j+(k-1)]
            }
            averageRs[i,k] <- sumR/((nrow(corMat)/nvars)-i)
        }
    }
    print(averageRs)
}


plotCors <- function(cors) {
    cors <- as.data.frame(cors)
    cors <- cors %>%
        mutate(lag=row_number()) %>%
        pivot_longer(cols=starts_with("V"))

    ggplot(aes(x=lag, y=value, group=name, color=name), data=cors) +
        geom_smooth(method="loess", se=FALSE, span=.4) +
        ylim(0, 1)
}



