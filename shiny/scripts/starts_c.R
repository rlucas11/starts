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
wx2 ~ a*wx1
wx3 ~ a*wx2
wx4 ~ a*wx3
wx5 ~ a*wx4
wx6 ~ a*wx5
wx7 ~ a*wx6
wx8 ~ a*wx7
wx9 ~ a*wx8
wx10 ~ a*wx9
wy2 ~ b*wy1
wy3 ~ b*wy2
wy4 ~ b*wy3
wy5 ~ b*wy4
wy6 ~ b*wy5
wy7 ~ b*wy6
wy8 ~ b*wy7
wy9 ~ b*wy8
wy10 ~ b*wy9
##
##
## Cross-lags
wy2 ~ c*wx1
wy3 ~ c*wx2
wy4 ~ c*wx3
wy5 ~ c*wx4
wy6 ~ c*wx5
wy7 ~ c*wx6
wy8 ~ c*wx7
wy9 ~ c*wx8
wy10 ~ c*wx9
wx2 ~ d*wy1
wx3 ~ d*wy2
wx4 ~ d*wy3
wx5 ~ d*wy4
wx6 ~ d*wy5
wx7 ~ d*wy6
wx8 ~ d*wy7
wx9 ~ d*wy8
wx10 ~ d*wy9
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
x1 ~~ s*x1
x2 ~~ s*x2
x3 ~~ s*x3
x4 ~~ s*x4
x5 ~~ s*x5
x6 ~~ s*x6
x7 ~~ s*x7
x8 ~~ s*x8
x9 ~~ s*x9
x10 ~~ s*x10
y1 ~~ s*y1
y2 ~~ s*y2
y3 ~~ s*y3
y4 ~~ s*y4
y5 ~~ s*y5
y6 ~~ s*y6
y7 ~~ s*y7
y8 ~~ s*y8
y9 ~~ s*y9
y10 ~~ s*y10
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
## Constraints
arx2 == 1 - ((arx * a^2 + ary * d^2) / (1 - rr^2))
ary2 == 1 - ((ary * b^2 + arx * c^2) / (1 - rr^2))
'
