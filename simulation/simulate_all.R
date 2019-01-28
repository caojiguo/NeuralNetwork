######################################################################
#
#
#
#
#
#
#
######################################################################


#Clear workspace
rm(list=ls())


#Load dependencies
library(R.matlab)
library(fda)
library(MASS)
library(parallel)
source('ode_functions.R')
###################
#Preparation
numberofinterior = c(6,3,3,3,3,3,3,2,5,4,5,3)
bs.list = list()
basislist = list()
derlist   = list()
secondlist = list()
time.start = 0 
time.end   = 20000
sporder = 4
timeseq = seq(time.start,time.end, length.out = 20001)
for (jj in 1:12){
  knots  = seq(time.start, time.end, length.out = numberofinterior[jj]+2)
  nbasis =  length(knots) + sporder - 2
  bs.list[[jj]] = create.bspline.basis(c(0,20000),nbasis, sporder,knots)
}

coeflist = list()
coeflist[[1]]    = c(-5.1307,-5.4168,-5.4352,-4.7943,-5.2277,-5.1568,-3.0955,-6.6271,-9.3371,-10.6941)
coeflist[[2]]    = c(  -8.4529 ,  -6.2820 ,  -6.2715 ,  -6.2126   ,-9.9787 ,  -4.1193 ,  -7.9767)
coeflist[[3]]    = c(-20.9127  ,-17.7111  , -9.8625  , -6.7971 , -10.1419 ,  -3.9452 , -11.4773)
coeflist[[4]]    = c(-6.0241   ,-8.6272   ,-5.5635  , -6.3896  , -6.3152  , -3.3511   ,-5.8368)
coeflist[[5]]    = c(-11.3194  ,-10.0080   ,-7.0482  , -7.0632  , -9.7340 ,  -6.2688  , -6.8309)
coeflist[[6]]    = c(-5.4294   ,-6.4731  , -7.5013  , -4.6526  , -5.1576   ,-4.5361   ,-3.8525)
coeflist[[7]]    = c(-7.8624   ,-6.6731   ,-4.3848   ,-8.1334  , -2.2517  , -7.6259, -109.9612)
coeflist[[8]]    = c(-16.4025  , -5.9799  , -7.4510   ,-5.9575 , -15.6495 , -26.8798)
coeflist[[9]]    = c(  -6.2949  ,  -6.8799 ,  -5.8076  , -5.1444 ,  -6.0723  , -3.1912  , -4.9163,  -10.8792 , -13.2856)
coeflist[[10]]   = c( -6.2077 ,  -7.1254  , -7.6184  , -6.4745 ,  -5.9687  , -5.3781 ,  -5.5484 ,  -3.9116)
coeflist[[11]]   = c( -6.2247 ,  -5.6359  , -5.6657  , -5.8015 ,  -6.1338 ,  -4.8073  , -5.4679 ,  -8.6507,  -10.0717)
coeflist[[12]]   = c( -15.3299,  -12.8632,   -6.2573  , -6.6603 , -13.1295  , -8.2579  , -5.9124)

neuron_num = 12
fd.list = list()
for(ii in 1:neuron_num){
  fd.list[[ii]] = fd(coef = coeflist[[ii]], basisobj = bs.list[[ii]])
}


true_intmat = matrix(NA, nrow = length(timeseq), ncol = neuron_num)
for(kk in 1:neuron_num){
  true_intmat[,kk] = exp(eval.fd(timeseq, fd.list[[kk]]))
}

#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
sim_result <- list()
spikingsum <- list()
temp       <- list()
for (index in 1:neuron_num){
  
  temp[[index]] <- simulate_NHPP(true_intensity = true_intmat[,index], time_start =0 , time_end = 20000, numofsim = 1000, seed = 1657945)
  sim_result[[index]] <- unlist(temp[[index]])
  spikingsum[[index]] <- unlist(lapply(temp[[index]], length))
  
}
unlist(lapply(spikingsum,mean))

true_spiking = readMat('truespiking.mat')$temp
index = 5
truespk = as.vector(unlist(true_spiking[[index]]))
plot(true_intmat[,index], type='l')
for (ii in 1:length(truespk)){
  segments(x0 = truespk[ii], x1 = truespk[ii], y0 = 0 , y1 = 0.001,col='blue')
}
for (ii in 1:length(temp[[index]])){
  segments(x0 = temp[[index]][[2]][ii], x1 = temp[[index]][[2]][ii], y0 = 0.001 , y1 = 0.002,col='red')
}



for (index in 1:12){
  name = paste('sim_neuron',index,'.mat',sep='')
  writeMat(name, sim_result = sim_result[[index]], spikingsum = spikingsum[[index]])
}

#------------------------------------------------------------------------------------------------------------------