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

#############
# Problem 3 #
#############

library('MASS')
library(MCMCpack)
library(xtable)

# Use post.f(), post.f.neg(), and bayes.logreg() defined
# in problem 2.

# Read data:
dat=read.table('breast_cancer.txt',header=TRUE)

# Extract X and y:
y=dat$diagnosis
# In y, set 'M'=1 and 'B'=0:
y=as.numeric(y=='M')
m=length(y)
n=rep(1,m)
X0=as.matrix(dat[,-11])
# Standardize X:
X=sapply(1:ncol(X0),function(j) (X0[,j]-mean(X0[,j]))/sd(X0[,j]))
# Add vector of ones to X matrix:
X=cbind(n,X)
p=ncol(X)

#################################################
# Set up the specifications:
beta.0=rep(0,p)
Sigma.0.inv=solve(diag(1000,p))
niter=15000
burnin=5000
print.every=1000
retune=burnin/10
verbose=TRUE
#################################################

# Fit the Bayesian model:
beta=bayes.logreg(n,y,X,beta.0,Sigma.0.inv,niter,burnin,
                  print.every,retune,verbose,p,m)

# Write/save beta matrix to a file:
write.matrix(beta,file='beta.csv',sep=',')

##########
### 3b ###
##########

set.seed(421)

beta=as.matrix(read.csv('beta.csv'))
# Remove burnin:
beta=beta[-c(1:burnin),]
n=nrow(beta)
# Convert to mcmc object:
beta=mcmc(beta)
# Number of predictive datasets (response vectors):
M=1000

# Ignore! plot(beta)
# Ignore! effectiveSize(beta)

# Autocorrelations:
acor=acf(beta,lag.max=1,plot=FALSE)$acf
# Subset lag 1 autocorrelations:
acor.mx=sapply(1:p,function(k) acor[,,k][2,])
# Latex:
xtable(acor.mx)

##########
### 3c ###
##########

# Posterior intervals:
perc=sapply(1:p,function(j) quantile(beta[,j],c(.025,.975)))
# Latex:
xtable(perc)

##########
### 3d ###
##########

# Sample M beta^(t)'s from posterior
# equivalent to sample M indices:
ind=sample(1:n,M)

# For each beta^(t), simulate a predictive response vector
# from p(y|beta^(t)):
# Remember probabilities are logit^(-1)(x_i^T beta).
prob=sapply(ind,function(ind) 
  sapply(1:nrow(X),function(i) 
    exp(X[i,]%*%beta[ind,])/(1+exp(X[i,]%*%beta[ind,]))))
# For posterior predictive check, only need statistic,
# not actual predictive datasets.
# Calculate statistic, number of successess,
# i.e. number of malignant:
suc=sapply(1:M,function(j) sum(rbinom(m,1,prob[,j])))

# Make histogram:
png('sta250hw1_3d.png')
hist(suc,xlab='Number of Malignant',main=
       paste('Histogram of Number of Malignant
        in',M,'Simulated Datasets')
)
abline(v=sum(y),lty=2,col='blue')
legend('topright',legend='Number of malignant\nfrom real dataset',
       lty=2,col='blue',bty='n')
dev.off()