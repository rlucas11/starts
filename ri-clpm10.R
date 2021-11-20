ri_clpm10 <- '
## Random Intercepts
##
ri_x =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5 + 1*x6 + 1*x7 + 1*x8 + 1*x9 + 1*x10
ri_y =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 + 1*y5 + 1*y6 + 1*y7 + 1*y8 + 1*y9 + 1*y10
##
## Create within-person centered variables
  wx1 =~ 1*x1
  wx2 =~ 1*x2
  wx3 =~ 1*x3 
  wx4 =~ 1*x4
  wx5 =~ 1*x5
  wx6 =~ 1*x6
  wx7 =~ 1*x7
  wx8 =~ 1*x8 
  wx9 =~ 1*x9
  wx10 =~ 1*x10
  wy1 =~ 1*y1
  wy2 =~ 1*y2
  wy3 =~ 1*y3
  wy4 =~ 1*y4
  wy5 =~ 1*y5
  wy6 =~ 1*y6
  wy7 =~ 1*y7
  wy8 =~ 1*y8
  wy9 =~ 1*y9
  wy10 =~ 1*y10
##
## Regressions
##
## Stabilities
wx2 ~ wx1
wx3 ~ wx2
wx4 ~ wx3
wx5 ~ wx4
wx6 ~ wx5
wx7 ~ wx6
wx8 ~ wx7
wx9 ~ wx8
wx10 ~ wx9
wy2 ~ wy1
wy3 ~ wy2
wy4 ~ wy3
wy5 ~ wy4
wy6 ~ wy5
wy7 ~ wy6
wy8 ~ wy7
wy9 ~ wy8
wy10 ~ wy9
##
##
## Cross-lags
wy2 ~ wx1
wy3 ~ wx2
wy4 ~ wx3
wy5 ~ wx4
wy6 ~ wx5
wy7 ~ wx6
wy8 ~ wx7
wy9 ~ wx8
wy10 ~ wx9
wx2 ~ wy1
wx3 ~ wy2
wx4 ~ wy3
wx5 ~ wy4
wx6 ~ wy5
wx7 ~ wy6
wx8 ~ wy7
wx9 ~ wy8
wx10 ~ wy9
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
wx4 ~~ wx4
wy4 ~~ wy4
wx5 ~~ wx5
wy5 ~~ wy5
wx6 ~~ wx6
wy6 ~~ wy6
wx7 ~~ wx7
wy7 ~~ wy7
wx8 ~~ wx8
wy8 ~~ wy8
wx9 ~~ wx9
wy9 ~~ wy9
wx10 ~~ wx10
wy10 ~~ wy10
##
## Covariances
ri_x ~~ ri_y
wx1 ~~ wy1
wx2 ~~ wy2
wx3 ~~ wy3
wx4 ~~ wy4
wx5 ~~ wy5
wx6 ~~ wy6
wx7 ~~ wy7
wx8 ~~ wy8
wx9 ~~ wy9
wx10 ~~ wy10
'
