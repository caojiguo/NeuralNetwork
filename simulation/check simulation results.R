rm(list=ls())
library(fda)


load('coef_realresult.Rdata')
load('sorted_curves_all.Rdata')



source('ode_functions.R')
neuron_num   <- 12
numberofinterior <- c(6,3,3,3,3,3,3,2,5,4,5,3)
time.start   <- 0 
time.end     <- 20000
sporder       <- 4
timeseq      <- seq(time.start, time.end, length.out = 101)
bs.list      <- list()
#Setting up basis objects for intensity curves
for (jj in 1:12){
  knots         <- seq(time.start, time.end, length.out = numberofinterior[jj] + 2)
  nbasis        <- length(knots) + sporder - 2
  bs.list[[jj]] <- create.bspline.basis(c(0,20000), nbasis, sporder, knots)
}

#coefficients for basis objects 
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
#-----------------------------------------------------------------------

fd.list = list()
for(ii in 1:neuron_num){
  fd.list[[ii]] = fd(coef = coeflist[[ii]], basisobj = bs.list[[ii]])
}
samplesize = 101
support = c(0,1)
num.knots = 11

sp.order = 4
time_seq = seq(0,20000,length.out = samplesize)
knots = seq(support[1],support[2],len = num.knots)
basis_obj = create.bspline.basis(rangeval = support, nbasis = num.knots +sp.order - 2, norder = sp.order, breaks = knots)

X_t = matrix(NA, nrow = length(time_seq), ncol = neuron_num)
for(kk in 1:neuron_num){
  X_t[,kk] = exp(eval.fd(time_seq, fd.list[[kk]]))
}
#Normalize intensity 
Z_t = apply(X_t, 2, function(x) { (x - min(x)) / (max(x) - min(x))})
#---
#plot all ci
for (neuron in 1:12){
  filename = paste('neuron',neuron,'_simresult.Rdata',sep='')
  load(filename)
  
  neuron.coef   <- temp
  phi_list = PHI_mat()
  regulation_mat = list()
  
  for (kk in 1:12){
    regulation_mat[[kk]] = matrix(NA, nrow = 101,ncol = 100)
  }
  
  
  for (iter in 1:100){
    
    sim_coef_vec <- neuron.coef[,iter]
    coef_mat = matrix(sim_coef_vec, nrow = 13, ncol = 12,byrow=FALSE)
    curves <- matrix(NA, nrow = 101, ncol = 12)
    
    par(mfrow=c(3,2))
    for (jj in 1:12){
      curves[,jj] = as.vector(coef_mat[,jj] %*% phi_list[[jj]])
      plot_data = curves[,jj][order(Z_t[,jj])]
      #plot(sort(X_t[,jj]), plot_data,type='l',lwd=2,main = paste('g = ',jj,', l = ',9),ylab=paste('f_',jj,',',9,sep=''),xlab=paste('X_',jj,sep=''))
      regulation_mat[[jj]][,iter] <- plot_data
    }
  }
  
  pointwise_info <- list()
  
  for (kk in 1:12){
    pointwise_info[[kk]] = matrix(NA, nrow = 3, ncol = 101)
    pointwise_info[[kk]][2,] = apply(regulation_mat[[kk]],1,mean)
    
    variance_vec = apply(regulation_mat[[kk]],1,var)
    pointwise_info[[kk]][1,] = pointwise_info[[kk]][2,] - 1.96*sqrt(variance_vec)
    pointwise_info[[kk]][3,] = pointwise_info[[kk]][2,] + 1.96*sqrt(variance_vec)
  }
  
  library(ggplot2)
  library(cowplot)
  plot <- list()
  for (k in 1:12){
    
    plot.data <- cbind(sort(X_t[,k]), true_reg[[neuron]][,k],
                       pointwise_info[[k]][2,], pointwise_info[[k]][1,],
                       pointwise_info[[k]][3,])
    colnames(plot.data) <- c('X','true','mean','lower','upper')
    temp_knots <- seq(min(X_t[,k]),max(X_t[,k]),length.out = 4)[2:4]
    temp_lab   <- scales::scientific(temp_knots,digits = 1)
    plot.data <- as.data.frame(plot.data)
    plot[[k]] <- ggplot(data = plot.data,aes(x = X)) + geom_line(aes(y = mean),size=1.1) +
      geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.45) +
      theme_light()+
      scale_y_continuous(name="",breaks=c(-4,-2,0,2,4),limits=c(-12,12)) +
      theme(axis.text = element_text(face = "bold",size = 10))+
      scale_x_continuous(name="Intensity",labels = c(0,temp_lab),
                         limits =c(0,NA),breaks = c(0,temp_knots) )
  }
  
  plotname <- paste('neuron',neuron,'_reg.png',sep='')
  png(file=plotname,height = 1320,width=810,units='px')
  plot_grid(plot[[1]],plot[[2]],
            plot[[3]],plot[[4]],
            plot[[5]],plot[[6]],
            plot[[7]],plot[[8]],
            plot[[9]],plot[[10]],
            plot[[11]],plot[[12]],nrow = 4,labels=1:12)
  ggsave(filename=plotname,width=8.13,height = 13.75,units='in')
  dev.off()
}

