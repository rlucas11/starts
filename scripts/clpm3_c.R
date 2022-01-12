clpm3 <- '
## Regressions
##
x2 ~ a*x1
y2 ~ b*y1
x2 ~ c*y1
y2 ~ d*x1
x3 ~ a*x2
y3 ~ b*y2
x3 ~ c*y2
y3 ~ d*x2
##
## Variances
x1 ~~ x1
y1 ~~ y1
x2 ~~ x2
y2 ~~ y2
x3 ~~ x3
y3 ~~ y3
##
## Covariances
x1 ~~ y1
x2 ~~ y2
x3 ~~ y3
'
