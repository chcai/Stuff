### runs pmcmc.fun

library(msm) # rtnorm()
library(mvtnorm) # rmvnorm()
library(coda) # mcmc

source('tn_fun.R')
source('pmcmc_fun.R')

file.ind = '00'
GPU = TRUE
# must have this if-condition for general code
# non-GPU's may not even be able to install RCUDA
if(GPU) library(RCUDA)

dat = as.matrix(read.table(paste0('data/data_', file.ind, '.txt'), 
                           header = TRUE))
y = dat[,1]
X = dat[,-1]
p = ncol(X)

beta0 = rep(0, p)
sig0inv = matrix(0, p, p)

burn = 5000
n.iter = 20000+burn

# assign beta matrix and system.time
time. = system.time({
  beta = pmcmc.fun(y, X, beta0, sig0inv, n.iter, burn, GPU)
})

# if want mcmc object, dump burn period
beta = as.mcmc(beta[burn:n.iter,])

if(!GPU) what = 'R'
if(GPU) what = 'C'

# save beta mcmc object and time as R objects
save(beta, file = paste0('results/beta', what, file.ind, '.rda'))
save(time., 
     file = paste0('results/pmcmcTime', what, file.ind, '.rda'))