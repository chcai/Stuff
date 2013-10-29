##########
### 3b ###
##########

library(coda)
library(MCMCpack)
library(xtable)
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