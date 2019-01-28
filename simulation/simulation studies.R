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


#load dependencies
library(R.matlab)
library(fda)
library(MASS)
source('ode_functions.R')



truedata <-  get_truefd()
numberofinterior <- c(6,3,3,3,3,3,3,2,5,4,5,3)
sporder  <- 4
time.start <- 0
time.end   <- 20000
timeseq <- seq(time.start,time.end, length.out = 101)



#Settings for intensity curves 
amp_cand    <- c(1e6,1e7,1e7,1e6,2e7,1e6,1e6,1e7,1e6,1e6,5e6,1e7)
lambda_cand <- seq(0.01,0.3,0.01)
lambda_vec  <- c(0.14,0.31,0.44,0.35,0.39,0.29,0.37,0.45,0.16,0.4,0.21,0.29)
#lambda_cand  <- c(0.01,seq(0.05,0.5,0.05))

basisnum   <- 13
neuron_num <- 12


connection.list <- list()
connection.list[[1]] <- c(1,4,10)
connection.list[[2]] <- c(1,6,9,10,11)
connection.list[[3]] <- c(1,4,6,9,10,11,12)
connection.list[[4]] <- c(1,4,6,9,10,12)
connection.list[[5]] <- c(4,6,9,10,12)
connection.list[[6]] <- c(1,3,5,10,11,12)
connection.list[[7]] <- c(1,4,6,9,10,11)
connection.list[[8]] <- c(1,4,6,7,9,10)
connection.list[[9]] <- c(1,4,9,10,12)
connection.list[[10]] <- c(1,4,6,9,10,11,12)
connection.list[[11]] <- c(4,6,9,10)
connection.list[[12]] <- c(1,4,6,7,9,10,11,12)


simulationresult <- list()
penaltyresult    <- list()
selectionresult  <- matrix(NA, nrow = neuron_num, ncol = neuron_num)
truedata <-  get_truefd()
for (neuron in 1:neuron_num){
  print(paste('running on neuron #',neuron,sep=' '))
  #-----------------------------------------------------
  filename <- paste('coef_neuron',neuron,'.mat',sep='')
  simdata  <- readMat(filename)
  coef_simulated <- simdata$coef
  #-----------------------------------------------------
  
  true_coef <- truedata$true.list[[neuron]]$coef 
  #define knots, order of spline functions, and basis functions 
  knots      <- seq(time.start, time.end, length.out = numberofinterior[neuron] + 2)
  sporder   <- 4
  nbasis    <- length(knots) + sporder - 2
  #define basis 
  basis     <- create.bspline.basis(c(time.start, time.end), nbasis, sporder, knots)
  #-----------------------------------------------------
  deriv.mat = matrix(NA, nrow = length(timeseq), ncol = 1000)
  for (jj in 1:1000){
    fundata <- fd(coef = coef_simulated[,jj], basisobj = basis)
    deriv.mat[,jj] = eval.fd(timeseq, fundata,1) * exp(eval.fd(timeseq, fundata))
  }
  filter    <- apply(deriv.mat, 2, function(x) { mean((x - truedata$true.der[,neuron])^2)})
  
  
  #filter = apply(coef_simulated, 2,  function(x){sum(abs(x - true_coef))})
  
  sim.result <- coef_simulated[,order(filter)]
  sim.result <- sim.result[,1:100]
  
  sim.fdlist <- list()
  #combine simulated coefficients and basis object into functional data object 
  for (jj in 1:ncol(sim.result)){
    sim.fdlist[[jj]] <- fd(coef = sim.result[,jj], basisobj = basis)
  }
  
  plot(truedata$true.der[,neuron], type='l', lwd = 2, col = 'red')
  for( jj in 1:ncol(sim.result)){
    lines(eval.fd(timeseq, sim.fdlist[[jj]],1) * exp(eval.fd(timeseq, sim.fdlist[[jj]])), col='blue',lwd = 1.5)
  }
  
  
  #
  temp <- sim_estimation(whichneuron = neuron, coef_sim = sim.result, lambda_cand = lambda_vec[neuron],amp_vec=amp_cand,
                            connection_vec = connection.list[[neuron]])
  
  sim_filename = paste('neuron',neuron,'_simresult.Rdata',sep='')
  
  save(sim.result, temp, file = sim_filename)
  
  simulationresult[[neuron]] <- temp
  
  
  selection = matrix(NA, nrow = 12, ncol = ncol(sim.result))
  phi.list = PHI_mat()
  for (jj in 1:ncol(sim.result)){
    selection[,jj] = checkconnection(coef_vec = simulationresult[[neuron]][,jj], phi_list = phi.list)
  }
  #
  selectionresult[,neuron] <- apply(selection,1,sum)
}

simulationresult <- list()
for (neuron in 1:12){
  filename = paste('neuron',neuron,'_simresult.Rdata',sep='')
  load(filename)
  simulationresult[[neuron]] <- temp
  selection = matrix(NA, nrow = 12, ncol = ncol(sim.result))
  phi.list = PHI_mat()
  for (jj in 1:ncol(sim.result)){
    selection[,jj] = checkconnection(coef_vec = simulationresult[[neuron]][,jj], phi_list = phi.list)
  }
  #
  selectionresult[,neuron] <- apply(selection,1,sum)
  
}




phi.list = PHI_mat()
for (jk in 1:12){
  plotname <- paste('Neuron_',jk,'_simulation','.png',sep='')
  #point-wise mean and 95% limits 
  sim_temp = matrix(NA, nrow = 101, ncol = 100)
  for (qq in 1:100){
    sim_temp[,qq] = as.vector(simulationresult[[jk]][,qq] %*% do.call(rbind,phi.list))
  }
  mean_vec <- apply(sim_temp, 1, mean)
  limits   <- t(apply(sim_temp,1, function(x) { return(quantile(x, probs = c(0.025, 0.975)))}))
  response <- get_response()
  png(plotname)
  plot(response[,jk]*amp_cand[jk],type='l',lwd = 2, ylab='Derivative of intensity')
  lines(mean_vec, col = 'red',lwd = 1.5)
  lines(limits[,1], col = 'blue')
  lines(limits[,2], col = 'blue')
  legend('topleft',lwd=c(2,1.5,1),col=c('black','red','blue'), legend = c('True derivative','Point-wise average','Point-wise 95% limits'),bty='n')
  dev.off()
}
