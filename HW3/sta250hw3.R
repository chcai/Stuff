setwd("C:/Users/Christine/Google Drive/sta250/sta250hw3")

#########
### 1 ###
#########

# Bisection and Newton-Raphson are root finding algorithms
# so for maximization,
# find the roots of the first derivative
target.fun = function(x) 125/(2+x)-38/(1-x)+34/x
target.deriv.fun = function(x) -125/(2+x)^2-38/(1-x)^2-34/x^2

# input l and u st fun(l)fun(u) < 0
bisec.fun = function(fun, l, u, tol, max.iter, verbose) {
  c. = (l+u)/2
  fun.c = fun(c.)
  # iteration counter for debugging
  i = 1
  # loop until fun(c.) < tol (converge) or hit max iterations
  while(abs(fun.c) > tol && i < max.iter) {
    fun.lc = fun(l)*fun(c.)
    if(fun.lc < 0) u = c. else l = c.
    fun.c = fun(c.)
    if(verbose && i%%(max.iter/10) == 0)
      print(paste0('iter=', i, ',', 'fun(c)=', fun.c))
    i = i+1
  }
  c.
}

nr.fun = function(fun, fun.deriv, x, tol, max.iter, verbose) {
  fun.x = fun(x)
  # iteration counter for debugging
  i = 1
  # loop until fun(x) < tol (converge) or hit max iterations
  while(abs(fun.x) > tol && i < max.iter) {
    x = x-fun(x)/fun.deriv(x)
    fun.x = fun(x)
    if(verbose) print(paste0('iter=', i, ',', 'fun(x)=', fun.x))
    i = i+1
  }
  x
}

bisec.fun(target.fun, .62, .63, 1e-06, 1000, TRUE)

nr.fun(target.fun, target.deriv.fun, .5, 1e-06, 1000, TRUE)