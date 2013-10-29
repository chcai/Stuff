#############
# Problem 3 #
#############

library('MASS')

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