summarizeR <- function(corMat, nvars) {

    averageRs <- matrix(nrow=(nrow(corMat)/nvars-1),
                        ncol=nvars)

    for (k in 1:nvars) {
        for (i in 1:((nrow(corMat)/nvars)-1)) {
            sumR <- 0
            for (j in seq(1, (nrow(corMat)-nvars*i), by=nvars)) {
                sumR <- sumR + corMat[(i+k+j+(i-1)), j+(k-1)]
            }
            averageRs[i,k] <- sumR/((nrow(corMat)/nvars)-i)
        }
    }
    print(averageRs)
}

corMat <- matrix(c(1, .5, .25,
                   .5, 1, .25,
                   .25, .5, 1),
                 nrow=3,
                 ncol=3)

summarizeR(corMat, nvars = 1)
