######################################################################
# Run parameter cascade method for estimating regulation function
# and reconstruct time-varying network of neurons 
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
#Load basis objects and coefficients
#Express intensity functions as linear combinations
numberofinterior = c(6,3,3,3,3,3,3,2,5,4,5,3)
bs.list = list()
basislist = list()
derlist   = list()
secondlist = list()
time.start = 0 
time.end   = 20000
sporder = 4
timeseq = seq(time.start,time.end, length.out = 201)
neuron_num = 12

#load smoothing spline result from 
#estimating intenstiy function in terms of 
#basis coefficients and objects
load('smoothingresult.Rdata')

#####################
#Set up for regulation functions
samplesize = 101
support = c(0,1)
num.knots = 11

sp.order = 4
time_seq = seq(0,20000,length.out = samplesize)
knots = seq(support[1],support[2],len = num.knots)
basis_obj = create.bspline.basis(rangeval = support, nbasis = num.knots +sp.order - 2, norder = sp.order, breaks = knots)

phi.list = get_PHImat(fd.list = fd.list,time_seq = time_seq, basis_obj = basis_obj, neuron_num = neuron_num,all.phi=FALSE)
#####################
#Get intensity trajectories
X_t = matrix(NA, nrow = length(time_seq), ncol = neuron_num)
for(kk in 1:neuron_num){
  X_t[,kk] = exp(eval.fd(time_seq, fd.list[[kk]]))
}

####################
#Plot intensity functions 
library(ggplot2)
library(cowplot)

true_spiking = readMat('spike_train_data.mat')$temp
plot <- list()
y1_vec = c(0.0008,0.000125,0.00007,0.00035,0.00006,0.001,
           0.00035,0.000035,0.0008,0.0007,0.0002,0.00015)
for (ii in 1:12){
  data.mat <- cbind(time_seq, X_t[,ii])
  colnames(data.mat) <- c('time','Intensity')
  rownames(data.mat) <- NULL
  data.mat <- as.data.frame(data.mat)
  spk.data <- as.vector(unlist(true_spiking[[ii]]))
  temp <- rep(y1_vec[ii], length(spk.data))
  spk.data <- cbind(spk.data, temp)
  colnames(spk.data) <- c('spikes','height')
  spk.data <- as.data.frame(spk.data)
  plot[[ii]] <- ggplot(data = data.mat , aes(x = time, y = Intensity)) + geom_line() + 
                          geom_segment(data = spk.data,aes(x = spikes, y = 0 ,
                                        xend = spikes, yend = height)) + theme_minimal() +
                                scale_y_continuous(name="") + scale_x_continuous(name="") + theme(axis.ticks.y=element_blank(),
                                                                                                  axis.text.y = element_blank())
  
}
plot_grid(plot[[1]],plot[[2]],
          plot[[3]],plot[[4]],
          plot[[5]],plot[[6]],
          plot[[7]],plot[[8]],
          plot[[9]],plot[[10]],
          plot[[11]],plot[[12]],nrow = 3,labels=1:12)





#Normalize intensity 
Z_t = apply(X_t, 2, function(x) { (x - min(x)) / (max(x) - min(x))})

#Get derivative of intensity for each neuron stored in response matrix  
response = matrix(NA,nrow = samplesize, ncol = neuron_num)
for (target in 1:neuron_num){
 temp = eval.fd(time_seq,fd.list[[target]],1) * exp(eval.fd(time_seq,fd.list[[target]]))
 temp_mean = mean(temp)
 response[,target] = temp-temp_mean
}


###############################
target_cand = seq(1,neuron_num,1)########
amp_cand    <- c(1e6,1e7,1e7,1e6,2e7,1e6,1e6,1e7,1e6,1e6,5e6,1e7)
#sparsity grid
lambda_cand = seq(0.01,0.5,0.01)
#lambda_cand  <- c(0.01,seq(0.05,0.5,0.05))
#list to store parameter estimations
par.list    = list()
#list to store information criteria w.r.t. each sparsity penalty 
info.list   = list()
#list t
leveling.list = list()
#
par.res = matrix(NA, nrow = basis_obj$nbasis * neuron_num, ncol = neuron_num)
sparsity.res = rep(NA, neuron_num)
for (ii in 1:neuron_num){
  print(paste("Neuron ",ii,sep=''))
  info.list[[ii]] = matrix(NA, nrow = 4, ncol = length(lambda_cand))
  rownames(info.list[[ii]]) = c('penalty','AIC','AICc','BIC')
  #columns correspond to penalty parameters 
  par.list[[ii]]  = matrix(NA, nrow = neuron_num * basis_obj$nbasis, ncol = length(lambda_cand))
  for (jj in 1:length(lambda_cand)){
    print(paste('Running on tuning #',jj,' of ',length(lambda_cand),sep=''))
    temp = getbeta(target = ii, basis_list = bs.list, num.knots = 11,amp=amp_cand[ii],coef_list = coeflist,lambda = lambda_cand[jj], alpha = 3.7,lambda_I = 100,gamma = 100000,iter.max = 300)
    par.list[[ii]][,jj] = temp$coef_est
    info.list[[ii]]['penalty',jj] = lambda_cand[jj]
    info.list[[ii]]['AIC',jj]     = temp$AIC
    info.list[[ii]]['AICc',jj]     = temp$AICc
    info.list[[ii]]['BIC',jj]     = temp$BIC
  }
  leveling.list[[ii]] <- matrix(NA, nrow = 3, ncol = length(lambda_cand))
  rownames(leveling.list[[ii]])    <- c('penalty','AICc','change')
  
  leveling.list[[ii]]['penalty',]  <- rev(lambda_cand)
  leveling.list[[ii]]['AICc',]     <- rev(info.list[[ii]]['AICc',])
  
  leveling.list[[ii]]['change',1]  <- 0
  
  for (kk in 2:length(lambda_cand)){
    leveling.list[[ii]]['change',kk] = (leveling.list[[ii]]['AICc',kk-1] - leveling.list[[ii]]['AICc',kk])/abs(leveling.list[[ii]]['AICc',kk])
  }
  index <- which(leveling.list[[ii]]['change',] == max(leveling.list[[ii]]['change',]))
  par.res[,ii] <- par.list[[ii]][,index]
  sparsity.res[ii] <- leveling.list[[ii]]['penalty',index]
  plot_res(response_der = response, target = target_cand[ii],coef_vec = par.res[,ii], phi_list = phi.list, Z_t = Z_t, time_seq = time_seq, amp = amp_cand[ii], neuron_num = neuron_num, basis_num = basis_obj$nbasis,lambda =  sparsity.res[ii])
}
dev.off()

