rm(list=ls())


library(igraph)
library('RColorBrewer')
library(fda)
library(MASS)
source('ode_functions.R')



#load data
load('leveloff.Rdata')


#network data setup

from = c()
to   = c()
weight = c()

name = c()
for (jj in 1:neuron_num){
  name[jj] = paste('neuron_',jj,sep='')
}


time_network = seq(0,20000,length.out = 201)
network = list()
phi.temp = get_PHImat(fd.list = fd.list, time_seq = time_network, basis_obj = basis_obj, neuron_num = 12, all.phi = FALSE)
for( ii in 1:length(time_network)){
  network[[ii]] = matrix(NA, nrow = 12 * 12, ncol = 3)
  colnames(network[[ii]]) = c('to','from','weight')
  network[[ii]][,'to']  = rep(name,each=12)
  network[[ii]][,'from'] = rep(name,times = 12)
  
  weight = c()
  weight_temp <- matrix(NA,nrow = 12, ncol = 12)
  for( jj in 1:neuron_num){
    temp = c()
    for (kk in 1:neuron_num){
      
      coef_vec = par.res[seq(1,156,13)[kk]:(seq(1,156,13)[kk]+12),jj]
      phi_vec = phi.temp[[jj]][,ii]
      temp[kk] = coef_vec %*% phi_vec
      weight_temp[kk,jj] <-  coef_vec %*% phi_vec
    }
    weight = c(temp,weight)
  }
  weight <- as.vector(weight_temp)
  network[[ii]][,'weight'] = round(weight,digits=3)
  network[[ii]] = data.frame(network[[ii]])
  network[[ii]][network[[ii]]$weight==0 ,'from'] = NA
  network[[ii]] = na.omit(network[[ii]])
  network[[ii]][,'weight'] = abs(as.numeric(as.character(network[[ii]][,'weight'])))/2
}

description = data.frame(id = name,type=seq(1,12,1))
#color palette
color.pal = brewer.pal(neuron_num,'Set3')

#confirm layout of netowrk (consistency)
#net = graph_from_data_frame(d = network[[1]], vertices = description)
#tk.id = tkplot(net)
#coords =  tkplot.getcoords(tk.id) 
#coords.backup = as.vector(coords)
#matrix(coords.backup, ncol = 2, byrow=FALSE)
coords = matrix(c(56, 386 , 98, 204, 373, 219, 340, 345, 311 ,296 , 36  ,92 ,237, 256,   0  ,30 ,167, 369 , 11, 369 ,103, 300, 109, 358),ncol = 2 , byrow = FALSE)
for (index in 1:201){
  plot.name = paste('time at ', time_network[index],'.pdf',sep='')
  pdf(file = plot.name)
  #network data
  net = graph_from_data_frame(d = network[[index]][,c(2,1,3)], vertices = description)
  
  
  #node color 
  V(net)$color <- color.pal[V(net)$type]
  #edge (link) color
  edge.start <- ends(net, es=E(net), names=F)[,1]
  E(net)$color <- V(net)$color[edge.start]
  #edge width 
  E(net)$width <- E(net)$weight
  
  ###
  #Static network visualization
  ##
  plot(net,layout = coords,vertex.frame.color='white',vertex.label.color="black",main= paste('time at ', time_network[index],'ms',sep=''))
  
  ###
  #Dynamic network visualization
  ##
  
  dev.off()
}





