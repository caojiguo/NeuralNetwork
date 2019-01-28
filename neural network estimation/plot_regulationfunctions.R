###########################################################
# Plot regulation functions for each neuron's differential equation
#
#
#
#
###########################################################

rm(list=ls())

source('ode_functions.R')
#Load data after running parameter cascade method 
load('leveloff.Rdata')
library(ggplot2)
library(cowplot)
#ii indicates the target neuron 
ii = 1;

  response_der = response
  target = target_cand[ii]
  coef_vec = par.res[,ii]
  phi_list = phi.list
  Z_t = Z_t
  time_seq = time_seq
  amp = amp_cand[ii]
  neuron_num = neuron_num
  basis_num = basis_obj$nbasis
  lambda =  sparsity.res[ii]
  
  
  
  nontrivial = length(coef_vec) - length(which(coef_vec == 0))
  phi_mat = do.call(rbind,phi_list)
  der_est = as.vector(coef_vec %*% phi_mat)
  
  # par(mfrow=c(1,1))
  # plot(time_seq,response_der[,target], col='red',type='l',lwd = 2, main='True and estimated derivative of intensity of neuron 9 along time',
  #      xlab='Time', ylab='Derivative')
  # lines(time_seq,der_est/amp, col='blue',lwd=2)
  # legend('topleft',lwd = c(2,2),col=c('red','blue'),c('True','Estimated'),bty='n')
  # 
  
  
  coef_mat = matrix(coef_vec, nrow = basis_num, ncol = neuron_num,byrow=FALSE)
  curves = matrix(NA, nrow = length(time_seq), ncol = neuron_num)
  sorted_curves = matrix(NA, nrow = length(time_seq), ncol = neuron_num)    
  
  # par(mfrow=c(3,2))
  for (jj in 1:neuron_num){
    curves[,jj] = as.vector(coef_mat[,jj] %*% phi_list[[jj]])
    #plot_data = curves[,jj][order(Z_t[,jj])]; sorted_curves[,jj] = plot_data
    # plot(sort(X_t[,jj]), plot_data,type='l',lwd=2,main = paste('g = ',jj,', l = ',target),ylab=paste('f_',jj,',',target,sep=''),xlab=paste('X_',jj,sep=''))
  }
  
  
  #-----------------------------
  #use ggplot 
  
  
  plot <- list()
  for (jj in 1:12){
    plot.data <- NULL
    curves[,jj] = as.vector(coef_mat[,jj] %*% phi_list[[jj]])
    temp = curves[,jj][order(Z_t[,jj])]; sorted_curves[,jj] = temp
    plot.data <- cbind(sort(X_t[,jj]),temp)
    colnames(plot.data) <- c('X','fX')
    rownames(plot.data) <- NULL
    plot.data <- as.data.frame(plot.data)
    plot[[jj]] <- ggplot(data = plot.data,aes(x = X,y = fX)) + geom_line(size = 1.15) + theme_bw() +  geom_line(aes(y=0),linetype="dashed") +
      scale_y_continuous(name="") + scale_x_continuous(name="Intensity") #+ theme(axis.ticks.x=element_blank(),
                                                                                               # axis.text.x = element_blank())
  }
  #png(paste('neuron',ii,'_reg.png',sep=''),width = 975,height = 613)
  plot_grid(plot[[1]],plot[[2]],
            plot[[3]],plot[[4]],
            plot[[5]],plot[[6]],
            plot[[7]],plot[[8]],
            plot[[9]],plot[[10]],
            plot[[11]],plot[[12]],nrow = 4,labels=1:12)
  #dev.off()
  
  
  
  
  
  
  
