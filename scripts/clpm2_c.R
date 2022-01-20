clpm2_c <- '
## Regressions
##
x2 ~ a*x1
y2 ~ b*y1
x2 ~ c*y1
y2 ~ d*x1
##
## Variances
x1 ~~ x1
y1 ~~ y1
x2 ~~ x2
y2 ~~ y2
##
## Covariances
x1 ~~ y1
x2 ~~ y2
'
