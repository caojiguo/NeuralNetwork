clear;

simdata = cell(12,1000);

check_1st = NaN(1,12);
check_2nd = NaN(1,12);
for index = [5,8,12]
    name = strcat('sim_neuron',num2str(index),'.mat');
    load(name)
    checkpoints = double(reshape(spikingsum,1,1000));
    for ii = 1:1000
        if ii == 1
            simdata{index,ii} = sim_result(1:checkpoints(ii));
        else
            simdata{index,ii} = sim_result(cumsum(checkpoints(1:(ii-1))+1) : cumsum(checkpoints(1:(ii-1))+checkpoints(ii))); 
        end
        
    end
    
    %check simulation 
    len_vec = NaN(1,1000);
    for ii = 1:1000
        len_vec(ii) = size(simdata{index,ii},1);
    end
    
    check_1st(index) = all(len_vec == checkpoints);
    
    temp = 0;
    for jj = 1:1000
        temp = size(simdata{index,jj},1) + temp;
    end
   
    check_2nd(index) = (size(sim_result,1) == temp);
end

check_1st
check_2nd
%simulated data are stored in 'simdata'

numofinterior = [6,3,3,3,3,3,3,2,5,4,5,3];
penalty = 10.^[9,4,8,5,9,3,6,6,9,9,9,9];

time_start = 0;
time_end   = 20000;
time_range  = [time_start,time_end];
time = time_start:1:time_end;

%Define B-spline basis 
%Order of spline
sp_order = 4; 
%config for parameter estimation 
options = optimset('MaxIter',50000,'MaxFunEvals',50000,'TolX',1e-3,'display','off','TolFun',1e-3);

result = cell(1,12);


for index = 1:12
    fprintf('Start on neuron #%d\n',index);
    interior_knots  = numofinterior(index);
    smoothpenalty   = penalty(index);
    fprintf('interior knots of %d\n',interior_knots);
    fprintf('smooth penalty of %d\n',smoothpenalty);
    
    knots     = linspace(time_start, time_end, interior_knots+2);
    nbasis    = length(knots) + sp_order - 2;
    bspline_basis = create_bspline_basis(time_range, nbasis, sp_order, knots); 
    basis_value   = eval_basis(time, bspline_basis);    
    basis_1st = eval_basis(time,bspline_basis,1);
    basis_2nd = eval_basis(time,bspline_basis,2);
    
    result_temp = NaN(nbasis, 1000);
    
    for jj =1:1000
        fprintf('Just start on  simuated train #%d\n', jj);
        starting = fminsearch(@(coef) penlik(coef,simdata{index,jj},0,nbasis,basis_2nd,basis_value),...
                                rand(1,nbasis),options);  
                            
        result_temp(:,jj) = fminsearch(@(coef) penlik(coef,simdata{index,jj},smoothpenalty,nbasis,basis_2nd,basis_value),...
                                    starting, options);     
    end
    

    result{1,index} = result_temp;
    
end

%save('backup.mat');
%save('smoothall.mat','result');

for some = [5,8,12]
     name = strcat('coef_neuron',num2str(some),'.mat'); 
     coef = result{1,some};
     save(name,'coef');
end




