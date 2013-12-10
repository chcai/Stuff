### this file is just for processing output
# i.e. making plots, etc

setwd("C:/Users/Christine/Google Drive/sta250/sta250hw4")

load('results/tnRtime.rda')
time.R = sapply(time., function(x) x[3])

load('results/tnCtime.rda')
time.C = sapply(time., function(x) x[3])

png('sta250hw4_1e.png')
plot(log(10^(1:8)), time.R, type = 'b', 
     xlab = 'log(n)', ylab = 'Seconds', 
     main = 'Runtimes')
points(log(10^(1:8)), time.C, type = 'b', lty = 2)
legend('topleft', legend = c('CPU', 'GPU'), lty = 1:2)
dev.off()

png('sta250hw4_1e2.png')
plot(log(10^(1:6)), time.R[1:6], type = 'b', 
     xlab = 'log(n)', ylab = 'Seconds', 
     main = 'Runtimes')
points(log(10^(1:6)), time.C[1:6], type = 'b', lty = 2)
legend('topleft', legend = c('CPU', 'GPU'), lty = 1:2)
dev.off()

alpha = (a-mu)/sig
beta = (b-mu)/sig
Z = pnorm(beta)-pnorm(alpha)
mu+(dnorm(alpha)-dnorm(beta))*sig/Z