### function for probit mcmc

pmcmc.fun = function(y, X, beta0, sig0inv, n.iter, burn, GPU) {
  n = length(y)
  p = ncol(X)
  n1 = sum(y) # number of successes
  n0 = n-n1 # number of failures
  
  # there are constants which will not update at each iteration
  # calculate them before iterating
  sig.b = solve(sig0inv+t(X)%*%X)
  C. = sig0inv%*%beta0
  
  # allocate vector/matrix
  z = numeric(n)
  beta = matrix(0, n.iter, p)
  beta[1,] = beta0
  
  # prep kernel inputs
  if(GPU) {
    a = numeric(n)
    a[y == 1] = 0
    a[y == 0] = -Inf
    b = numeric(n)
    b[y == 1] = Inf
    b[y == 0] = 0
  }
  
  # start iterations
  for(i in 2:n.iter) {
    # mean of z updates with beta
    mean.z = X%*%beta[i-1,]
    
    # sample z
    if(GPU) {
      z = tn.fun(n, z, mean.z, rep(1, n), a, b, i, TRUE)
    } else {
      z[y == 1] = rtnorm(n1, mean.z[y == 1], 1, 0, Inf)
      z[y == 0] = rtnorm(n0, mean.z[y == 0], 1, -Inf, 0)
    }
    
    # mean of beta updates with z
    mean.b = sig.b%*%(C.+t(X)%*%z)
    beta[i,] = rmvnorm(1, mean.b, sig.b)
  }
  # return beta
  beta
}