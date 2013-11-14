#####################
# BLB_lin_reg_job.R #
#####################

#============================== Setup for running on Gauss... ======#

args <- commandArgs(TRUE)

cat("Command-line arguments:\n")
print(args)

####
# sim_start ==> Lowest possible dataset number
###

###################
sim_start <- 1000
###################

sim_num <- sim_start + as.numeric(args[1])

cat(paste("\nAnalyzing dataset number ",sim_num,"...\n\n",sep=""))

# Find r and s indices:
s=5
r=50
s_index=ceiling(as.numeric(args[1])/r)
r_index=as.numeric(args[1])%%r
if(r_index==0) r_index=r

#============================== Run the simulation study ===========#

library(bigmemory)

# Attach big.matrix :
dat=attach.big.matrix('/home/pdbaines/data/blb_lin_reg_data.desc')

# Remaining BLB specs:
n=nrow(dat)
gamma=.7
b=ceiling(n^gamma)

# Extract the subset:
# Each job with same s_index should sample same sub-dataset
set.seed(s_index)
inds=sample(1:n,b,replace=FALSE)
dat.sub=dat[inds,]
y=dat.sub[,ncol(dat.sub)]
X=dat.sub[,-ncol(dat.sub)]

# Bootstrap dataset:
set.seed(as.numeric(args[1]))
weights.=rmultinom(1,n,rep(1/b,b))

# Fit lm:
fit=lm(y~0+.,data=as.data.frame(X),weights=weights.)$coef

# Output file:
outfile=paste0("output/","coef_",sprintf("%02d",s_index),"_",sprintf("%02d",r_index),".txt")
# Save estimates to file:
write(fit,file=outfile,ncolumns=1)