starts_uni <- '
## Random Intercepts
##
ri_x =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5 + 1*x6 + 1*x7 + 1*x8 + 1*x9 + 1*x10
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
##
##
##
## Variances
ri_x ~~ ri_x
wx1 ~~ arx*wx1
wx2 ~~ arx2*wx2
wx3 ~~ arx2*wx3
wx4 ~~ arx2*wx4
wx5 ~~ arx2*wx5
wx6 ~~ arx2*wx6
wx7 ~~ arx2*wx7
wx8 ~~ arx2*wx8
wx9 ~~ arx2*wx9
wx10 ~~ arx2*wx10
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
##
##
## Constraints
arx2 == 1 - (arx * a^2)
'
