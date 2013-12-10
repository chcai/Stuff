### this file is just for processing outputs
# i.e. making plots, etc

setwd("C:/Users/Christine/Google Drive/sta250/sta250hw4")

# library(coda)

file.ind = '00'
GPU = TRUE
if(!GPU) what = 'R'
if(GPU) what = 'C'

load(paste0('results/beta', what, file.ind, '.rda'))
png(paste0('sta250hw4_2d', what, file.ind, '.png'))
par(mfrow = c(3, 3), oma = c(0, 0, 3, 0))
sapply(1:ncol(beta), function(j) {
  plot(beta[,j], type = 'l', ylab = paste0('beta ', j))
})
if(file.ind == '00') {
  mtext(paste0('Mini Data ', what, ' Traceplots'), 3, 
        outer = TRUE)
} else {
  mtext(paste0('Data ', file.ind, ' ', what, ' Traceplots'), 3, 
        outer = TRUE)
}
dev.off()
round(apply(beta, 2, mean), 2)

load(paste0('results/pmcmcTime', what, file.ind, '.rda'))
time.

time.R = c(9.197, 61.492, 691.622, 5048.904, 47903.002)
time.C = c(121.385, 132.683, 250.912, 1440.813, 13863.011)

png('sta250hw4_2e.png')
plot(1:5, time.R/3600, type = 'b', 
     xlab = 'Dataset', ylab = 'Hours', 
     main = 'Runtimes')
points(1:5, time.C/3600, type = 'b', lty = 2)
legend('topleft', legend = c('CPU', 'GPU'), lty = 1:2)
dev.off()