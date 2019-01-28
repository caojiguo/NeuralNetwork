# Estimating Time-varying Directed Neural Networks
We have included all the codes for estimating intensity functions with smoothing splines, estimating time-varying neural networks, and simulation studies in their folders respectively.

Codes and functions for estimating intensity function of neutrons are stored in folder ‘intensity function estimation’. For this folder,
1. ‘Estimate_intensity.m’  is the main code for estimating intensity function.
2. ‘cvcheck.m’ computes cross validation score to choose optimal smoothing parameter.
3. ‘penlik.m’ is the function corresponding to the penalized log likelihood function as the sum of log likelihood function and roughness penalty.
4. ‘spktime.mat’ is the matrix contains all spike times of 12 neurons.
5. ‘fdam.zip’. Please unzip this file first. This folder contains matlab functions for doing functional data analysis, it is equivalent to the ‘fda’ package in R.


Folder ‘neural network estimation’ contains codes and functions for estimating time-varying directed networks:
1. ‘ode_model.R’ is the main code for estimating networks.
2. ‘ode_functions.R’ contains the dependencies and functions for ‘ode_model.R’ .
3. ‘smoothingresult.Rdata’ stores the basis coefficients and objects from using smoothing spline method for estimating intensity functions.
4. ‘spike_train_data.mat’ contains all spike sequences of 12 neurons.
5. ‘animated_network.mp4’ is the visualized and animated dynamic network for 20000ms.
6. ‘animated_network.pdf’ contains time snap shots of the network.
7. ‘leveloff.Rdata’ stores the estimated regulation function for each differential equation after choosing an appropriate sparsity parameter lambda.
8. ‘plot_regulationfunctions.R’ can bed used for plotting regulation functions for each differential equations._
9. ‘network_visual.R’ is used to take time snap shots of the network._

Folder ‘simulation’ includes codes for simulating spike sequences and simulated data. Also, codes for estimating neural networks based on simulated data are included. For this folder,
1. ‘neuron_sim_result’ folder contains the plot of estimated intensity functions of 100 simulated spike sequences of all 12 neutrons.
2. ‘estimated intensity functions’ includes basis coefficients for all 100 simulated spike sequences of each neuron.
3. ‘simulated spike sequences’ includes all 100 simulated spike sequences used for simulation studies.
4. ‘fdam.zip’. Please unzip this file first. This folder contains matlab functions for doing functional data analysis, it is equivalent to the ‘fda’ package in R. 
5. ‘coefest_simulation.m’ contains the codes for obtaining basis coefficients stored in ‘estimated intensity functions’._
6.  ‘ penlik.m’ is the function corresponding to the penalized log likelihood function as the sum of log likelihood function and roughness penalty.
7. ‘simulate_all.R’ can simulate spike sequences based on thinning methodology of Lewis and Shedler (1979). Changing ‘seed = ‘ can generate a new random spike sequence.
8. ‘simulation studies.R’ contains codes for estimating neural networks based on simulated data. Selection table in the paper is generated frome this code.
9.  ‘ode_functions.R’ contains the dependencies and functions for R scripts.
10. ‘check simulation studies.R’ contains codes to check the results of simulation studies including calculation and plotting of confidence intervals for regulation functions of each differential equation by changing ‘neuron = ‘.
11. Series of ‘neuronx_simresult.Rdata’ corresponding to the estimation of regulation functions given simulated data of neuron x.
12. ‘coef_realresult.Rdata’ cotains the basis coefficients of all regulation function for all neurons based on real data.
13. ‘sorted_curves_all.Rdata’ contains all the estimated regulation functions for all neurons. This is for ploting with bootstrap confidence intervals.
