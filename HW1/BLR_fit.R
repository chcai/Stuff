#############
# Problem 2 #
#############

##
#
# Logistic regression
# 
# Y_{i} | \beta \sim \textrm{Bin}\left(n_{i},e^{x_{i}^{T}\beta}/
# (1+e^{x_{i}^{T}\beta})\right)
# \beta \sim N\left(\beta_{0},\Sigma_{0}\right)
#
##

# library(mvtnorm)
# library(coda)
library('MASS')

#####################################################################
#####################################################################
## Handle batch job arguments:

# 1-indexed version is used now.
args <- commandArgs(TRUE)

cat(paste0("Command-line arguments:\n"))
print(args)

####
# sim_start ==> Lowest simulation number to be analyzed by this 
# particular batch job
###

#######################
sim_start <- 1000
length.datasets <- 200
#######################

if (length(args)==0){
  sinkit <- FALSE
  sim_num <- sim_start + 1
  set.seed(1330931)
} else {
  # Sink output to file?
  sinkit <- TRUE
  # Decide on the job number, usually start at 1000:
  sim_num <- sim_start + as.numeric(args[1])
  # Set a different random seed for every job number!!!
  set.seed(762*sim_num + 1330931)
}

# Simulation datasets numbered 1001-1200

#####################################################################
#####################################################################

# This function calculates the log posterior:
post.f=function(beta.t,n,y,X,beta.0,Sigma.0.inv,m){
  tmp=sapply(1:m,function(i) -n[i]*log(1+exp(X[i,]%*%beta.t)))
  tmp=sum(tmp)
  tmp2=sapply(1:m,function(i) y[i]*X[i,]%*%beta.t)
  tmp2=sum(tmp2)
  tmp+(tmp2-1/2*t(beta.t-beta.0)%*%Sigma.0.inv%*%(beta.t-beta.0))
}

# Will use optim() to initialize beta^(0). optim() minimizes.
# Minimize negative log posterior = maximize log posterior.
# This function just gives negative of log posterior:
post.f.neg=function(beta.t,n,y,X,beta.0,Sigma.0.inv,m){
  -post.f(beta.t,n,y,X,beta.0,Sigma.0.inv,m)
}

bayes.logreg=function(n,y,X,beta.0,Sigma.0.inv,niter,burnin,
                      print.every,retune,verbose,p,m){
  # Allocate beta matrix:
  beta=matrix(nrow=niter+burnin,ncol=p)
  # Initialize beta^(0):
  beta[1,]=optim(beta.0,post.f.neg,n=n,y=y,X=X,beta.0=beta.0,
                 Sigma.0.inv=Sigma.0.inv,m=m)$par
  # Initialize proposal variances:
  v=rep(1,p)
  # Counter for acceptance rates:
  i=rep(0,p)
  # Ignore! if(p==11) {v[c(8,11)]=1e-03; v[7]=1e-04; v[2]=1e-06}
  
  # Loop through each time state:
  for(t. in 1:(niter+burnin-1)){
    beta[t.+1,]=beta[t.,]; beta.star=beta[t.+1,]
    # Loop through each beta_j:
    for(j in 1:p){
      # Sample beta_j^* from proposal:
      tmp=rnorm(1,beta[t.,j],sqrt(v[j]))
      beta.star[j]=tmp
      
      # Calculate log(alpha), the Metropolis ratio:
      top=post.f(beta.star,n,y,X,beta.0,Sigma.0.inv,m)
      bot=post.f(beta[t.+1,],n,y,X,beta.0,Sigma.0.inv,m)
      # Take difference instead of quotient
      # since using log scale:
      alpha=top-bot
      
      # Sample from U(0,1):
      u=runif(1)
      
      # Accept or reject beta_j^*?
      if(log(u)<alpha){
        # Accept and update:
        beta[t.+1,j]=tmp
        # Add 1 to acceptance counter:
        i[j]=i[j]+1
      }else{
        # Reject:
        beta.star[j]=beta[t.,j]}
      
      # Retune proposal variance every 'retune' iterations
      # and only during burnin:
      if(t.%%retune==0&&t.<burnin){
        # If verbose = TRUE, print acceptance counter:
        if(verbose) print(
          sprintf('t = %0.f: acceptance rate (j = %0.f) = %0.f',
                  t.,j,i[j]))
        # Retune process:
        if(i[j]/retune<0.1) v[j]=v[j]/4
        if(i[j]/retune>=0.1&&i[j]/retune<0.2) v[j]=v[j]/3
        if(i[j]/retune>=0.2&&i[j]/retune<0.3) v[j]=v[j]/2
        if(i[j]/retune>=0.6&&i[j]/retune<0.7) v[j]=v[j]*2.1
        if(i[j]/retune>=0.7&&i[j]/retune<0.8) v[j]=v[j]*3.1
        if(i[j]/retune>=0.8&&i[j]/retune<0.9) v[j]=v[j]*4.1
        if(i[j]/retune>=0.9) v[j]=v[j]*5.1
        # Reset j-th acceptance counter for next retune period:
        i[j]=0
      }
    }
    if(verbose&&t.%%print.every==0) 
      print(list('t'=t.,'beta'=beta[t.,],'proposal variances'=v))
  }
  # Return beta matrix:
  beta
}

# Read data corresponding to appropriate sim_num:
dat=read.csv(sprintf('~/Stuff/HW1/BayesLogit/data/blr_data_%.0f.csv',
             sim_num))
# If testing on local machine:
# dat=read.csv('blr_data_1001.csv')

# Extract X and y:
y=dat$y
m=length(y)
n=dat$n
X=as.matrix(dat[,3:4])
p=ncol(X)

#################################################
# Set up the specifications:
beta.0=rep(0,p)
Sigma.0.inv=solve(diag(p))
niter=10000
burnin=1000
print.every=1000
retune=burnin/10
verbose=TRUE
#################################################

# Fit the Bayesian model:
beta=bayes.logreg(n,y,X,beta.0,Sigma.0.inv,niter,burnin,
                  print.every,retune,verbose,p,m)

# Extract posterior quantiles...
beta.q=
  sapply(1:p,function(j) quantile(beta[,j],seq(.01,.99,by=.01)))

# Write results to a (99 x p) csv file...
write.matrix(beta.q,file=
  sprintf('~/Stuff/HW1/BayesLogit/results/blr_res_%.0f.csv',
  sim_num),sep=',')

cat("done. :)\n")