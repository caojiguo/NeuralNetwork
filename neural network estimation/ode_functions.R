#####################################################################
#
#  Functions used in 'ode_model.R' for running parameter cascading
#
#
#
#
#
#####################################################################


getbeta <- function(target, basis_list, coef_list, amp, support = c(0,1), sp.order = 4, num.knots = 11, lambda, lambda_I = 100, alpha = 3.7, gamma = 1e5, iter.max = 500){
	#Input:
	#    target : index for indicating target neuron 
	# basis_list: list of basis objects for smoothe intensity trajectories
	#  coef_list: list of coefficients for the basis_list 
	#   sp.order: order of spline functions for f_gl 
	#  num.knots: number of knots including boundary and interior knots 
	#   lambda  : sparsity penalty in fSCAD penalty 
	# lambda_I  : Identifiability parameter
	#    alpha  : sparsity penalty in fSCAD penalty 
	#    gamma  : smoothness penalty for f_gl 
	#   iter.max: maximum number of iterations for computing GAM coefficients for f_gl 
	# 
	#Output: 
	#   coef_est : matrix with estimated coefficients stored in columns for each iteration 
	#     AIC    : AIC calculated with coefficients at last iteration 
	#    AICc    : AICc calculated with coefficients at last iteration 
	#     BIC    : BIC calculated with coefficients at last iteration 
	#
	#Comment: 

	#-----------------------------------------------------
	#Set up uniform basis object for all regulation functions f_gl 
	time_seq <- seq(0,20000,length.out = 101)
	knots <- seq(support[1], support[2], len = num.knots)
	basis_obj <- create.bspline.basis(rangeval = support, nbasis = num.knots + sp.order - 2, norder = sp.order, breaks = knots)

	#number of neurons 
	neuron_num <- length(basis_list)
	#transform basis lists and coef lists into functional data object 
	fd.list = list()
	for(ii in 1:neuron_num){
		fd.list[[ii]] = fd(coef = coef_list[[ii]], basisobj = basis_list[[ii]])
	}
	#-----------------------------------------------------

	#-----------------------------------------------------
	#Calculating derivative of intensity for target neuron 
	target_der = eval.fd(time_seq,fd.list[[target]],1) * exp(eval.fd(time_seq,fd.list[[target]]))
	#center the response 
	der_mean = mean(target_der)
	target_der = target_der - der_mean 
	#amplify the response by amp 
	target_der = target_der * amp
    #-----------------------------------------------------

    #sample size of intensity: n 
	numofpoints = length(time_seq )

	#-----------------------------------------------------
	#Start estimating coefficients 
	#store result each iteration in columns 
    coef_est = matrix(NA, nrow = basis_obj$nbasis * neuron_num, ncol = iter.max)
  
	for (iter in 1:iter.max){
	  	if (iter == 1){
		    phi_mat = get_PHImat(fd.list = fd.list, time_seq = time_seq, basis_obj = basis_obj, neuron_num = neuron_num,all.phi = TRUE)
			OMEGA = get_OMEGAmat(fd.list = fd.list,time_seq = time_seq, basis_obj = basis_obj, neuron_num = neuron_num)
			coef_est[,iter] = ginv(phi_mat %*% t(phi_mat) / numofpoints + lambda_I * OMEGA) %*% phi_mat %*% target_der / numofpoints
		}else{
     		coef_last = coef_est[,iter-1]
		    W_mat = get_Wmat(fd.list = fd.list,time_seq = time_seq, basis_obj = basis_obj, coef =  coef_last, neuron_num = neuron_num,lambda = lambda , alpha = alpha,knots = knots)$W_mat
		    (betanonzero = get_Wmat(fd.list = fd.list,time_seq = time_seq, basis_obj = basis_obj, coef =  coef_last, neuron_num = neuron_num,lambda = lambda , alpha = alpha,knots = knots)$betanonzero)
	
		  
		    W_mat = W_mat[betanonzero,betanonzero]
		    phi_mat = get_PHImat(fd.list = fd.list, time_seq = time_seq, basis_obj = basis_obj, neuron_num = neuron_num,all.phi = TRUE)
		    phi_mat = phi_mat[betanonzero,]
		  
		    V_mat = get_Vmat(fd.list = fd.list,time_seq = time_seq, basis_obj = basis_obj, neuron_num = neuron_num)
		    V_mat = V_mat[betanonzero,betanonzero]
			
			O_mat = get_OMEGAmat(fd.list = fd.list,time_seq = time_seq, basis_obj = basis_obj, neuron_num = neuron_num)
			O_mat = O_mat[betanonzero,betanonzero]
			
			coef_temp = ginv(phi_mat %*% t(phi_mat) / numofpoints + gamma * V_mat + W_mat + lambda_I * O_mat) %*% phi_mat %*% target_der / numofpoints
			coef_est[betanonzero,iter] = coef_temp
			coef_est[!betanonzero,iter] = 0
		
			if(iter == iter.max){
			  #Calculate infomration criterion given the coefficients at last iteration 
			  temp = sum(diag(t(phi_mat) %*% ginv(phi_mat %*% t(phi_mat) / numofpoints + gamma * V_mat + W_mat + lambda_I * O_mat) %*% phi_mat))/numofpoints
			  phi_origin = get_PHImat(fd.list = fd.list, time_seq = time_seq, basis_obj = basis_obj, neuron_num = neuron_num,all.phi = TRUE)
			  AIC  = length(time_seq) * log(mean((target_der - as.vector(coef_est[,iter] %*% phi_origin) )^2))  + 2 * temp
			  AICc = AIC + 2*temp*(temp+1)/(length(target_der) - temp- 1)
			  BIC  = length(time_seq) * log(  mean((target_der - as.vector(coef_est[,iter] %*% phi_origin) )^2)) + temp * log(length(time_seq))
			}
	
		}
	}


  return(list(coef_est = coef_est[,iter.max],AIC=AIC,AICc=AICc,BIC=BIC))

}

get_Vmat <- function(fd.list,time_seq,basis_obj,neuron_num){
  #Input:
  #     fd.list  : 
  #    time_seq  :
  #    basis_obj :
  #  neuron_num  : 
  #	
  #Output:
  #
  #     V_mat    : 
  #
  #Comment:
  #
  #
  #
  
  
  
	#Evaluate intensity trajectories
	X_t = matrix(NA, nrow = length(time_seq), ncol = neuron_num)
	for(kk in 1:neuron_num){
		X_t[,kk] = exp(eval.fd(time_seq, fd.list[[kk]]))
	}
	#Normalize intensity 
	Z_t = apply(X_t, 2, function(x) { (x - min(x)) / (max(x) - min(x))})
	distance = apply(X_t, 2, function(x){return(max(x) - min(x))})

	#Derivative 
	dX_dt = matrix(NA, nrow = length(time_seq), ncol = neuron_num)
	for (ii in 1:neuron_num){
		dX_dt[,ii] = eval.fd(time_seq, fd.list[[ii]],1) * exp(eval.fd(time_seq, fd.list[[ii]]))
	}
	dZ_dt = matrix(NA, nrow = length(time_seq), ncol = neuron_num)
	for (jj in 1:neuron_num){
		dZ_dt[,jj] = dX_dt[,jj] / distance[jj]
	}

	#Second derivative 
	d2X_dt = matrix(NA, nrow = length(time_seq), ncol = neuron_num) 
	for (ii in 1:neuron_num){
		d2X_dt[,ii] = eval.fd(time_seq, fd.list[[ii]],2) * exp(eval.fd(time_seq, fd.list[[ii]])) + (eval.fd(time_seq, fd.list[[ii]],1)^2) * exp(eval.fd(time_seq, fd.list[[ii]]))
	}
	d2Z_dt = matrix(NA, nrow = length(time_seq), ncol = neuron_num)
	for (kk in 1:neuron_num){
		d2Z_dt[,kk] = d2X_dt[,kk] / (distance[kk]^2)
	}
	
	
	#Evaluation of integral by quadrature
	weight = c(1,rep(c(4,2),(length(time_seq)-3)/2),4,1)*(length(time_seq)-1)/(length(time_seq)-1)/3
	#V_gl 
	V_gl = list()
	for(jj in 1:neuron_num){
		temp1 = matrix(rep(dZ_dt[,jj],basis_obj$nbasis), ncol =basis_obj$nbasis, byrow =FALSE)
		temp2 = matrix(rep(d2Z_dt[,jj],basis_obj$nbasis), ncol = basis_obj$nbasis, byrow = FALSE)
    
		der_mat = eval.basis(basis_obj, Z_t[,jj], 1)
		sec_mat = eval.basis(basis_obj, Z_t[,jj],2)
		
		d_mat = sec_mat * (temp1^2) + der_mat * temp2
		V_gl[[jj]] = t(d_mat)%*%(d_mat*(weight%*%t(rep(1,basis_obj$nbasis))))
	}
	V_mat = as.matrix(bdiag(V_gl))
	return(V_mat)
}	

get_PHImat <- function(fd.list,time_seq, basis_obj, neuron_num,all.phi=TRUE){
  #Input:
  #    fd.list  :
  #   time_seq  :
  #   basis_obj :
  #  neuron_num :
  #   all.phi   :
  #Output:
  #
  #		phi_mat
  #     phi_gl
  #
  #Comment:
  #
  #
  #
  #
  
  #
  #Evaluate intensity trajectories
  #Intensities are stored in columns of X_t 
	X_t = matrix(NA, nrow = length(time_seq), ncol = neuron_num)
	for(kk in 1:neuron_num){
		X_t[,kk] = exp(eval.fd(time_seq, fd.list[[kk]]))
	}

	#Normalize intensity to (0,1)
	Z_t = apply(X_t, 2, function(x) { (x - min(x)) / (max(x) - min(x))})

	#For each g-th normalized intensity ,evaluate each basis function and store in rows of matrix phi_gl 
	#phi_gl can be thought of as basis values given a basis object for each neuron intensity 
	phi_gl = list()
	for (ii in 1:neuron_num){
		phi_gl[[ii]] = t(eval.basis(basis_obj, Z_t[,ii]))
	}
	#Stack phi_gl matrices to get phi_l matrix 
	phi_mat = do.call(rbind, phi_gl)
	#given logical value all.phi to determin which form of basis value to return 
  if(all.phi){
    return(phi_mat)
  }else{
    return(phi_gl)
  }
}

get_OMEGAmat <- function(fd.list, time_seq, basis_obj, neuron_num){
  #Input:
  #    fd.list  :
  #   time_seq  :
  #   basis_obj :
  #  neuron_num :
  #
  #Output:
  #
  #   omega_mat
  #
  #Comment:
  #
  #
  #
  #
  
  
  
	#Evaluate intensity trajectories
	X_t = matrix(NA, nrow = length(time_seq), ncol = neuron_num)
	for(kk in 1:neuron_num){
		X_t[,kk] = eval.fd(time_seq, fd.list[[kk]])
	}
	#Normalize intensity 
	Z_t = apply(X_t, 2, function(x) { (x - min(x)) / (max(x) - min(x))})

	omega_gl = list()
	for( ii in 1:neuron_num){
		omega_glk = apply(eval.basis(Z_t[,ii], basis_obj),2,sum)
		omega_gl[[ii]] = omega_glk %*% t(omega_glk)
	}
	omega_mat = as.matrix(bdiag(omega_gl))
	return(omega_mat)
}

get_Wmat <- function(fd.list, time_seq, knots,basis_obj, coef, neuron_num,lambda, alpha,threshold = 1e-5){
  #Input:
  #
  #    fd.list:
  #   time_seq:
  #      knots: 
  #  basis_obj:
  #       coef:
  # neuron_num:
  #     lambda:
  #      alpha:
  #  threshold:
  #
  #Output:
  #
  #         W_mat:
  #  beta_nonzero:
  #
  #Comment:
  #
  #
  #
  #
  
  
	#Evaluate intensity trajectories
	X_t = matrix(NA, nrow = length(time_seq), ncol = neuron_num)
	for(kk in 1:neuron_num){
		X_t[,kk] = eval.fd(time_seq, fd.list[[kk]])
	}
	#Normalize intensity 
	Z_t = apply(X_t, 2, function(x) { (x - min(x)) / (max(x) - min(x))})

	#int.temp contains the M_glj matrix for j = 1,..., M_gl
	int.temp = eval.quad(knots = knots, numofpoints = 200, basis_obj = basis_obj)

	#
	M_gl = length(knots) - 1
	coef_mat = matrix(coef, nrow = basis_obj$nbasis, ncol = neuron_num, byrow = FALSE)
	
	#L-2 norm of basis functions over each sub-interval for each neuron 
    int.l2norm = matrix(NA, nrow = M_gl, ncol=neuron_num)
    for (jj in 1:neuron_num){
      for (ii in 1:M_gl){
        int.l2norm[ii,jj] = sqrt(coef_mat[,jj]%*%int.temp[[ii]]%*%coef_mat[,jj])
      }
    }
    
    
    #Derivativer of fSCAD penalty given function norm in each sub-interval (1 to M_gl) for each neuron (1 to neuron_num)
    pen.der = matrix(NA, nrow = M_gl, ncol = neuron_num)
    for (kk in 1:neuron_num){
      for (qq in 1:M_gl){
        pen.der[qq,kk] = fSCADpen(lambda = lambda, alpha = alpha, x = int.l2norm[qq,kk] * sqrt(M_gl), der.return=TRUE)
      }
    }
    
    order = basis_obj$nbasis - length(basis_obj$params)
    
    
    #finding 0 norm of functions and reducing dimension of W_gl and subsequentely W_l matrix 
    betanonzero = rep(TRUE, length(coef))
    W_gl  = list()
    for (jj in 1:neuron_num){
      W_gl[[jj]] = matrix(NA, nrow = basis_obj$nbasis, ncol=basis_obj$nbasis)
      temp = list()
      for(ii in 1:M_gl){
          if (int.l2norm[ii,jj] < threshold){
            betanonzero[ii:(ii+order-1) + (jj-1)*basis_obj$nbasis]=FALSE
            temp[[ii]] = matrix(0,nrow = basis_obj$nbasis,ncol=basis_obj$nbasis)
          }else{
              temp[[ii]] = pen.der[ii,jj] * sqrt(M_gl) * int.temp[[ii]] / int.l2norm[ii,jj]
          }
      }
      W_gl[[jj]] = 0.5 * Reduce('+',temp)
    }
    W_mat = as.matrix(bdiag(W_gl))
    
    


    # W_gl = list()
    # for (index in 1:neuron_num){
    # 	W_gl[[index]] = matrix(0, nrow = basis_obj$nbasis, ncol = basis_obj$nbasis)
    # 	temp = list()
    # 	for (iij in 1:M_gl){
    # 	  temp[[iij]] = pen.der[iij,index] * sqrt(M_gl) * int.temp[[iij]] / int.l2norm[iij,index]
    # 		#(W_gl[[index]] = W_gl[[index]] + pen.der[iij,index] * sqrt(M_gl) * int.temp[[iij]] / int.l2norm[iij,index])
    # 	}
    # 	W_gl[[index]] = 0.5 * Reduce('+',temp)
    # }
    # W_mat = as.matrix(bdiag(W_gl))
    return(list(W_mat = W_mat, betanonzero = betanonzero))
}

fSCADpen <- function(lambda,alpha = NULL, x, der.return = FALSE){
 
	#fSCADpen: Calculate the SCAD penalty or its derivative given tuning parameters
	#
	#Input:
	#       lambda   = regularization parameter
	#       alpha    = location parameter (default at 3.7)
	#          x     = norm of the coefficients for penalization
	#     der.return = logical value to determine if the function should return penalty or the derivative of penalty function
	#Output:
	#       penalty  = evaluation of SCAD penalty function for a given norm of coefficients
	#       pen.der  = evaluation of derivative of SCAD penalty function for a given norm of coefficients 
	#Comment: 

  if(is.null(alpha)){
    alpha = 3.7
  }
  indi = (x <= lambda)
  penalty = lambda * x * indi + (-0.5 * pmax(alpha * lambda - x,0)^2/(alpha - 1) + 0.5 * (alpha + 1) * lambda^2) * (1 - indi)
  
  pen.der = lambda * indi + (pmax(alpha * lambda - x, 0)/(alpha - 1)) * (1 - indi)
  
  
  if (der.return){
    #print('Return the derivative of SCAD penalty function')
    return(pen.der)
  }else{
    #print('Returning SCAD penalty')
    return(penalty)
  }
}



plot_res <-  function(response_der,target,coef_vec, phi_list,Z_t, time_seq,amp,neuron_num,basis_num, vs.time = TRUE, vs.Zg = TRUE,lambda,test=FALSE,sim=FALSE,whichsim=NULL){
	#plot_res:
	#
	#
	#
	#
	#Input:
	#	resposne_der:
	#	target      :
	#   coef_vec    :
	#   phi_list    :
	#   Z_t         :
	#   time_seq    :
	#   amp         :
	#   neuron_num  :
	#   basis_num   :
	#   lambda      :
	#   vs.time     :
	#   vs.Zg       :
	#   test        :
	#   sim         :
	#   whichsim    :
	#   vs.time     :
	#   vs.Zg       :
	#   test        :
	#   sim         :
	#Output:
	#
	#
	#comment:
	#

  nontrivial = length(coef_vec) - length(which(coef_vec == 0))
  if(vs.time){
    phi_mat = do.call(rbind,phi_list)
    der_est = as.vector(coef_vec %*% phi_mat)
    if(sim){
      name = paste('Neuron_',target,'sim_',whichsim,'_lambda','_',lambda,'_nontrivial_',nontrivial,'.png',sep='')
    }else{
      name = paste('Neuron_',target,'_lambda','_',lambda,'_nontrivial_',nontrivial,'.png',sep='') 
    }
    if(!test){
      png(filename = name)
      par(mfrow=c(1,1))
      plot(time_seq,response_der[,target] * amp, col='red',type='l',lwd = 2, main='True derivative of intensity vs. estimated along time')
      lines(time_seq,der_est, col='blue',lwd=2)
      dev.off()
    }else{
      par(mfrow=c(1,1))
      plot(time_seq,response_der[,target] * amp, col='red',type='l',lwd = 2, main='True derivative of intensity vs. estimated along time')
      lines(time_seq,der_est, col='blue',lwd=2)
    }
  }
  if(vs.Zg){
    if(!test){
      coef_mat = matrix(coef_vec, nrow = basis_num, ncol = neuron_num,byrow=FALSE)
      curves = matrix(NA, nrow = length(time_seq), ncol = neuron_num)
      for (jj in 1:neuron_num){
        if (jj == 1 || jj == 7){
          if(sim){
            name = paste('Neuron_',target,'sim_',whichsim,'curves_',jj,'_lambda','_',lambda,'_nontrivial_',nontrivial,'.png',sep='')
          }
          else{
            name = paste('Neuron_',target,'curves_',jj,'_lambda','_',lambda,'_nontrivial_',nontrivial,'.png',sep='') 
          }
          png(filename = name)
          par(mfrow=c(3,2))
        }
        curves[,jj] = as.vector(coef_mat[,jj] %*% phi_list[[jj]])
        plot_data = curves[,jj][order(Z_t[,jj])]
        plot(sort(Z_t[,jj]), plot_data,type='l',lwd=2,main = paste('f_gl, g = ',jj,', l = ',target))
        if(jj == 6||jj == 12){
          dev.off()
        }
      }
    }else{
      par(mfrow=c(3,2))
      coef_mat = matrix(coef_vec, nrow = basis_num, ncol = neuron_num,byrow=FALSE)
      curves = matrix(NA, nrow = length(time_seq), ncol = neuron_num)
      for (jj in 1:neuron_num){
        curves[,jj] = as.vector(coef_mat[,jj] %*% phi_list[[jj]])
        plot_data = curves[,jj][order(Z_t[,jj])]
        plot(sort(Z_t[,jj]), plot_data,type='l',lwd=2,main = paste('f_gl, g = ',jj,', l = ',target))
      }
    }

  }
  
}

simulate_NHPP <- function(true_intensity, time_start = 0 , time_end = 20000, numofsim = 300,seed){
  
	#Simulate_NHPP:
	#
	#
	#
	#Input:
	#
	#
	#Output:
	#
	#
	#Comment:
	#
	#
	#
	#

  set.seed(seed)
  #Upper limit of intensity function
  lambda_upper = max(true_intensity)
  result = list()
  for (jj in 1:numofsim){
  
    spktrain_sim = c()
    time = 0 
    while (time <= time_end){
      u1 <- runif(1)
      time <- time - log(u1) / lambda_upper 
      if(time < 1){ time = 1}
      u2 <- runif(1)
      if(time <= time_end){
        if( u2  <= (true_intensity[round(time)] / lambda_upper)){
          spktrain_sim = append(spktrain_sim,round(time)) 
        } 
      }
    }
    result[[jj]] = spktrain_sim
  }
  
  
  return(result)
}

plot_spkint <- function(true_int = true_intensity, true_obs = true_spiking, sim_spiking){
  plot(true_int,type='l')
  for (ii in 1:length(true_spiking)){
    segments(x0 = true_obs[ii], x1 = true_obs[ii], y0 = 0 , y1 = 0.001,col='blue')
  }
  for (ii in 1:length(true_spiking)){
    segments(x0 = sim_spiking[ii], x1 = sim_spiking[ii], y0 = 0.003 , y1 = 0.004,col='red')
  }
}

ode_sim  <- function(whichsim,fd_sim,lambda_cand, amp=1e6,target = 1){
  
  #Preparation
  numberofinterior = c(6,3,3,3,3,3,3,2,5,4,5,3)
  bs.list = list()
  basislist = list()
  derlist   = list()
  secondlist = list()
  time.start = 0 
  time.end   = 20000
  sporder = 4
  timeseq = seq(time.start,time.end, length.out = 201)
  for (jj in 1:12){
    knots  = seq(time.start, time.end, length.out = numberofinterior[jj]+2)
    nbasis =  length(knots) + sporder - 2
    bs.list[[jj]] = create.bspline.basis(c(0,20000),nbasis, sporder,knots)
    basislist[[jj]]  = eval.basis(bs.list[[jj]],timeseq)
    derlist[[jj]]    = eval.basis(bs.list[[jj]],timeseq,1)
    secondlist[[jj]] = eval.basis(bs.list[[jj]], timeseq,2)
  }
  
  coeflist = list()
  coeflist[[1]]    = fd_sim[[whichsim]]$coefs
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
  ###################
  neuron_num = 12
  fd.list = list()
  for(ii in 1:neuron_num){
    fd.list[[ii]] = fd(coef = coeflist[[ii]], basisobj = bs.list[[ii]])
  }
  
  #####################
  samplesize = 101
  support = c(0,1)
  num.knots = 11
  sp.order = 4
  time_seq = seq(0,20000,length.out = samplesize)
  knots = seq(support[1],support[2],len = num.knots)
  basis_obj = create.bspline.basis(rangeval = support, nbasis = num.knots +sp.order - 2, norder = sp.order, breaks = knots)
  #
  
  phi.list = get_PHImat(fd.list = fd.list,time_seq = time_seq, basis_obj = basis_obj, neuron_num = neuron_num,all.phi=FALSE)
  
  #####################
  
  X_t = matrix(NA, nrow = length(time_seq), ncol = neuron_num)
  for(kk in 1:neuron_num){
    X_t[,kk] = exp(eval.fd(time_seq, fd.list[[kk]]))
  }
  #Normalize intensity 
  Z_t = apply(X_t, 2, function(x) { (x - min(x)) / (max(x) - min(x))})
  
  ######################
  response = matrix(NA,nrow = samplesize, ncol = neuron_num)
  for (target in 1:neuron_num){
    temp = eval.fd(time_seq,fd.list[[target]],1) * exp(eval.fd(time_seq,fd.list[[target]]))
    temp_mean = mean(temp)
    response[,target] = temp-temp_mean
  }
  
  #Trying tuning with AICc
  try_result = matrix(NA, nrow = 12*basis_obj$nbasis ,ncol = length(lambda_cand))
  info_res   <- matrix(NA, nrow = 4, ncol = length(lambda_cand))
  rownames(info_res) <- c('penalty','AIC','AICc','BIC')

  for (jj in 1:length(lambda_cand)){
    print(paste('Running on tuning #',jj,' of ',length(lambda_cand),sep=''))
    temp = getbeta(target = 1, basis_list = bs.list, amp = 1e6, coef_list = coeflist, lambda = lambda_cand[jj],alpha = 3.7, lambda_I = 100, gamma =100000, iter.max = 100)
    try_result[,jj] = temp$coef_est
    info_res['penalty',jj]    <- lambda_cand[jj]
    info_res['AIC',jj]         = temp$AIC 
    info_res['AICc',jj]        = temp$AICc
    info_res['BIC',jj]         = temp$BIC
  }
  
  leveling <- matrix(NA, nrow = 3, ncol = length(lambda_cand))
  rownames(leveling) <- c('penalty','AICc','change')
  leveling['penalty',] <- rev(lambda_cand)
  leveling['AICc',]    <- rev(info_res['AICc',])
  leveling['change',1] <- 0

  for (jk in 2:length(lambda_cand)){
  	leveling['change',jk] = (leveling['AICc',jk-1] - leveling['AICc',jk]) /abs(leveling['AICc',jk])
  }
  index <- which(leveling['change',] == max(leveling['change',]))

  coef_result = try_result[,index]

  
  
  plot_res(response_der = response, target = 1,coef_vec = coef_result,phi_list = phi.list, Z_t = Z_t, time_seq = time_seq, amp = 1e6, neuron_num = neuron_num, basis_num = 13,lambda = lambda_cand[index],sim=TRUE, whichsim = whichsim)
  
  return(list(coef_result = coef_result, penalty = lambda_cand[index]))
}

checkconnection <- function(coef_vec,phi_list,basis_num = 13, neuron_num = 12,simnum = 100){
  coef_mat = matrix(coef_vec, nrow = basis_num, ncol = neuron_num,byrow=FALSE)
  curves = matrix(NA, nrow = 101, ncol = neuron_num)
  for (jj in 1:neuron_num){
    curves[,jj] = as.vector(coef_mat[,jj] %*% phi_list[[jj]])
  }
  
  result = rep(NA, neuron_num)
  for (ii in 1:neuron_num){
    if (all(curves[,ii] == 0)){
      result[ii] = 0 
    }else{
      result[ii] = 1
    }
  }
  
  return(result)

}

PHI_mat <- function (abc = NULL){
  #Preparation
  numberofinterior = c(6,3,3,3,3,3,3,2,5,4,5,3)
  bs.list = list()
  basislist = list()
  derlist   = list()
  secondlist = list()
  time.start = 0 
  time.end   = 20000
  sporder = 4
  
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
  ###################
  neuron_num = 12
  fd.list = list()
  for(ii in 1:neuron_num){
    fd.list[[ii]] = fd(coef = coeflist[[ii]], basisobj = bs.list[[ii]])
  }
  #####################
  samplesize = 101
  support = c(0,1)
  num.knots = 11
  sp.order = 4
  time_seq = seq(0,20000,length.out = samplesize)
  knots = seq(support[1],support[2],len = num.knots)
  basis_obj = create.bspline.basis(rangeval = support, nbasis = num.knots +sp.order - 2, norder = sp.order, breaks = knots)
  #
  
  phi.list = get_PHImat(fd.list = fd.list,time_seq = time_seq, basis_obj = basis_obj, neuron_num = neuron_num,all.phi=FALSE)
  
  return(phi.list)
}

#evaluation integral for each subintervals via quadrature (simpson)
#Input:
#     knots       = knots of spline basis including the boundary points and interior knots  
#     numofpoints = number of points to evaluate integral via quadrature 
#     basis_obj   = Bspline basis objects 
#Output:
#     int.res     = list object containing the evaluation of integral of square of basis functions over different 
#                   sub-intervals. 
eval.quad <- function(knots, numofpoints = 200, basis_obj){
  
  int.res = list()
  
  numofsubs  = length(knots) -1 
  numofbasis = basis_obj$nbasis
  
  
  for (ii in 1:numofsubs){
    support.temp = c(knots[ii],knots[ii+1])
    
    quadpts      <- seq(support.temp[1], support.temp[2], length.out = numofpoints+1)
    nquadpts     <- length(quadpts)
    quadwts      <- as.vector(c(1,rep(c(4,2),(nquadpts-2)/2),4,1),mode="any")
    quadwts      <- c(1,rep(c(4,2),(nquadpts-1)/2))
    quadwts[nquadpts] <- 1
    quadwts      <- quadwts/ (3 * numofpoints)
    
    basis.temp   <- eval.basis(basis_obj,quadpts)
    
    int.res[[ii]]<- t(basis.temp)%*%(basis.temp*(quadwts%*%t(rep(1,numofbasis))))
    
  }
  return(int.res)
}


get_response = function(abc = NULL){
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
  timeseq = seq(time.start,time.end, length.out = 201)
  for (jj in 1:12){
    knots  = seq(time.start, time.end, length.out = numberofinterior[jj]+2)
    nbasis =  length(knots) + sporder - 2
    bs.list[[jj]] = create.bspline.basis(c(0,20000),nbasis, sporder,knots)
    basislist[[jj]]  = eval.basis(bs.list[[jj]],timeseq)
    derlist[[jj]]    = eval.basis(bs.list[[jj]],timeseq,1)
    secondlist[[jj]] = eval.basis(bs.list[[jj]], timeseq,2)
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
  #Normalize intensity 
  Z_t = apply(X_t, 2, function(x) { (x - min(x)) / (max(x) - min(x))})
  
  #Get derivative of intensity for each neuron stored in response matrix  
  response = matrix(NA,nrow = samplesize, ncol = neuron_num)
  for (target in 1:neuron_num){
    temp = eval.fd(time_seq,fd.list[[target]],1) * exp(eval.fd(time_seq,fd.list[[target]]))
    temp_mean = mean(temp)
    response[,target] = temp-temp_mean
  }
  
  return(response)
}