%--------------------------------------------------------------------------
%   Returns negative penalized log-likelihood of a point process
%   
%   Inputs: 
%
%   coef        :  coefficients of  B-splines 
%   spike_time  :  time of events (spikes)
%   pen_coef    :  penalty coefficients for smoothness
%   nbasis      :  number of basis function 
%   basis_value :  evaluation of basis functions along time
%   basis_2nd   : 2nd derivative of basis functions evaluated along time 
%--------------------------------------------------------------------------



function negloglik = penlik(coef,spike_time,pen_coef,nbasis,basis_2nd,basis_value)
    
    %Estimate intensity as linear combinations of basis function along time
    intensity_est = exp(basis_value * transpose(coef));
    %sum of log intensity at time of events 
    firstterm = sum(log(intensity_est(spike_time)));
    
    
    %Evaluate integral for the second term using Simpson's (quadrature) rule 
    m = 10000;
    quadwts = ones(2*m,1);
    quadwts(2*(1:(m-1))) = 4;
    quadwts(2*(1:m)+1) = 2;
    quadwts = quadwts(1:20001);
    quadwts(20001) = 1;
    %Evaluate integral 
    secondterm = (transpose(intensity_est) * quadwts)/3;
    
    %Evaluate the penalized second derivative integral 
    mat_temp = quadwts * ones(1,nbasis) / 3;
    not_temp =  times(basis_2nd, mat_temp);                                                                                                                                                                                                                                                                                                                                                                                                             mat_temp;
    Rmat       = transpose(basis_2nd) * not_temp;
    penalty = pen_coef * coef * Rmat * transpose(coef);
    
    %Calculate the negtive loglikelihood
    negloglik = -firstterm + secondterm + penalty;
end