clpm5_c <- '
## Regressions
##
## Wave 2
x2 ~ a*x1
y2 ~ b*y1
x2 ~ c*y1
y2 ~ d*x1
## Wave 3
x3 ~ a*x2
y3 ~ b*y2
x3 ~ c*y2
y3 ~ d*x2
## Wave 4
x4 ~ a*x3
y4 ~ b*y3
x4 ~ c*y3
y4 ~ d*x3
## Wave 5
x5 ~ a*x4
y5 ~ b*y4
x5 ~ c*y4
y5 ~ d*x4
##
## Variances
x1 ~~ x1
y1 ~~ y1
x2 ~~ x2
y2 ~~ y2
x3 ~~ x3
y3 ~~ y3
x4 ~~ x4
y4 ~~ y4
x5 ~~ x5
y5 ~~ y5
##
## Covariances
x1 ~~ y1
x2 ~~ y2
x3 ~~ y3
x4 ~~ y4
x5 ~~ y5
'
