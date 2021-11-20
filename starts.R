starts_c <- '
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
wx1 ~~ arx*wx1
wy1 ~~ ary*wy1
wx2 ~~ arx2*wx2
wy2 ~~ ary2*wy2
wx3 ~~ arx2*wx3
wy3 ~~ ary2*wy3
wx4 ~~ arx2*wx4
wy4 ~~ ary2*wy4
wx5 ~~ arx2*wx5
wy5 ~~ ary2*wy5
wx6 ~~ arx2*wx6
wy6 ~~ ary2*wy6
wx7 ~~ arx2*wx7
wy7 ~~ ary2*wy7
wx8 ~~ arx2*wx8
wy8 ~~ ary2*wy8
wx9 ~~ arx2*wx9
wy9 ~~ ary2*wy9
wx10 ~~ arx2*wx10
wy10 ~~ ary2*wy10
x1 ~~ x1
x2 ~~ x2
x3 ~~ x3
x4 ~~ x4
x5 ~~ x5
x6 ~~ x6
x7 ~~ x7
x8 ~~ x8
x9 ~~ x9
x10 ~~ x10
y1 ~~ y1
y2 ~~ y2
y3 ~~ y3
y4 ~~ y4
y5 ~~ y5
y6 ~~ y6
y7 ~~ y7
y8 ~~ y8
y9 ~~ y9
y10 ~~ y10
##
## Covariances
ri_x ~~ rr*ri_y
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
##
'
