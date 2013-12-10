### runs tn.fun

library(msm) # rtnorm()
library(RCUDA)

source('tn_fun.R')

set.seed(421)

k = 4
mu = 0
sig = 1
a = -Inf # lower truncation
b = -10 # upper truncation
GPU = TRUE

tn = tn.fun(10^k, numeric(10^k), rep(mu, 10^k), rep(sig, 10^k), 
            rep(a, 10^k), rep(b, 10^k), 1, GPU)
# tn
# mean(tn)

# save system.time output as R object
time. = lapply(k, function(k) system.time({
  tn.fun(10^k, numeric(10^k), mu, sig, a, b, GPU)
}))
save(time., file = 'tnCtime.rda')