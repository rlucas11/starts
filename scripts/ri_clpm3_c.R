ri_clpm3_c <- '
## Random Intercepts
##
ri_x =~ 1*x1 + 1*x2 + 1*x3
ri_y =~ 1*y1 + 1*y2 + 1*y3
##
## Create within-person centered variables
  wx1 =~ 1*x1
  wx2 =~ 1*x2
  wx3 =~ 1*x3 
  wy1 =~ 1*y1
  wy2 =~ 1*y2
  wy3 =~ 1*y3
##
## Regressions
##
## Stabilities
wx2 ~ a*wx1
wx3 ~ a*wx2
wy2 ~ b*wy1
wy3 ~ b*wy2
##
##
## Cross-lags
wy2 ~ c*wx1
wy3 ~ c*wx2
wx2 ~ d*wy1
wx3 ~ d*wy2
##
## Variances
ri_x ~~ ri_x
ri_y ~~ ri_y
wx1 ~~ wx1
wy1 ~~ wy1
wx2 ~~ wx2
wy2 ~~ wy2
wx3 ~~ wx3
wy3 ~~ wy3
##
## Covariances
ri_x ~~ ri_y
wx1 ~~ wy1
wx2 ~~ wy2
wx3 ~~ wy3
'
