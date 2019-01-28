%----------------------------------------------------
% Estimate intensity functions of neurons by
# smoothing spline method
%
%
%
%
%
%
%----------------------------------------------------





clear;

%load spike train data
load('spktime.mat');

global time;           %time axis  %starting points for optimization, can be changed
       %options for fminsearch



time_start = 0;
time_end   = 20000;
time_range  = [time_start,time_end];
time = time_start:1:time_end;

%Define B-spline basis
%Order of spline
sp_order = 4;
options = optimset('MaxIter',50000,'MaxFunEvals',50000,'TolX',1e-3,'display','off','TolFun',1e-3);


knots_temp = [2,3,4,5,6,7,8,9,10];
penalty_range = 10.^[1,2,3,4,5,6,7,8,9];

%column for each train
coef_train = cell(1,12);
coef_basis = cell(1,12);
coef_1st   = cell(1,12);
coef_2nd   = cell(1,12);




for train = 1:12
   fprintf('Train #%d\n', train);
    coef_train{train} = cell(4,8);
    coef_basis{train} = cell(4,8);
    coef_1st{train}   = cell(4,8);
    coef_2nd{train}   = cell(4,8);

    for it = 1:9
        fprintf('Knots #%d\n', knots_temp(it));
        num_knots = knots_temp(it);
        knots     = linspace(time_start, time_end, num_knots+2);
        nbasis    = length(knots) + sp_order - 2;
        bspline_basis = create_bspline_basis(time_range, nbasis, sp_order, knots);
        basis_value   = eval_basis(time, bspline_basis);
        basis_1st = eval_basis(time,bspline_basis,1);
        basis_2nd = eval_basis(time,bspline_basis,2);

        for jt = 1:9
            fprintf('Penalty #%d\n', penalty_range(jt));
            coef_basis{train}{it,jt} = basis_value;
            coef_1st{train}{it,jt}   = basis_1st;
            coef_2nd{train}{it,jt}   = basis_2nd;

            starting = fminsearch(@(coef) penlik(coef,temp{train},0,nbasis,basis_2nd,basis_value),...
                                rand(1,nbasis),options);

            coef_train{train}{it,jt} = fminsearch(@(coef) penlik(coef,temp{train},penalty_range(jt),nbasis,basis_2nd,basis_value),...
                                    starting, options);

        end
    end

    save('refine_3.mat');
end

%Store intensity functions for each neuron
%with 1st and 2nd derivatives
intensity_mat(1,:) = exp(coef_train{1}{5,9} * transpose(coef_basis{1}{5,9}));
intensity_mat(2,:) = exp(coef_train{2}{2,4} * transpose(coef_basis{2}{2,4}));
intensity_mat(3,:) = exp(coef_train{3}{2,8} * transpose(coef_basis{3}{2,8}));
intensity_mat(4,:) = exp(coef_train{4}{2,5} * transpose(coef_basis{4}{2,5}));
intensity_mat(5,:) = exp(coef_train{5}{2,9} * transpose(coef_basis{5}{2,9}));
intensity_mat(6,:) = exp(coef_train{6}{2,3} * transpose(coef_basis{6}{2,3}));
intensity_mat(7,:) = exp(coef_train{7}{2,6} * transpose(coef_basis{7}{2,6}));
intensity_mat(8,:) = exp(coef_train{8}{1,6} * transpose(coef_basis{8}{1,6}));
intensity_mat(9,:) = exp(coef_train{9}{4,9} * transpose(coef_basis{9}{4,9}));
intensity_mat(10,:) = exp(coef_train{10}{3,9} * transpose(coef_basis{10}{3,9}));
intensity_mat(11,:) = exp(coef_train{11}{4,9} * transpose(coef_basis{11}{4,9}));
intensity_mat(12,:) = exp(coef_train{12}{2,9} * transpose(coef_basis{12}{2,9}));

deriv_mat(1,:) = (coef_train{1}{5,9} * transpose(coef_1st{1}{5,9})) .* intensity_mat(1,:);
deriv_mat(2,:) = (coef_train{2}{2,4} * transpose(coef_1st{2}{2,4})) .* intensity_mat(2,:);
deriv_mat(3,:) = (coef_train{3}{2,8} * transpose(coef_1st{3}{2,8})) .* intensity_mat(3,:);
deriv_mat(4,:) = (coef_train{4}{2,5} * transpose(coef_1st{4}{2,5})) .* intensity_mat(4,:);
deriv_mat(5,:) = (coef_train{5}{2,9} * transpose(coef_1st{5}{2,9})) .* intensity_mat(5,:);
deriv_mat(6,:) = (coef_train{6}{2,3} * transpose(coef_1st{6}{2,3})) .* intensity_mat(6,:);
deriv_mat(7,:) = (coef_train{7}{2,6} * transpose(coef_1st{7}{2,6})) .* intensity_mat(7,:);
deriv_mat(8,:) = (coef_train{8}{1,6} * transpose(coef_1st{8}{1,6})) .* intensity_mat(8,:);
deriv_mat(9,:) = (coef_train{9}{4,9} * transpose(coef_1st{9}{4,9})) .* intensity_mat(9,:);
deriv_mat(10,:) = (coef_train{10}{3,9} * transpose(coef_1st{10}{3,9})) .* intensity_mat(10,:);
deriv_mat(11,:) = (coef_train{11}{4,9} * transpose(coef_1st{11}{4,9})) .* intensity_mat(11,:);
deriv_mat(12,:) = (coef_train{12}{2,9} * transpose(coef_1st{12}{2,9})) .* intensity_mat(12,:);


save('finalresult.mat','intensity_mat', 'deriv_mat','temp');





clear;
load('finalresult.mat');
time_start = 0;
time_end   = 20000;
time_range  = [time_start,time_end];
time = time_start:1:time_end;
for train = 1:12
    figure

    subplot(2,1,1);
    plot(time(50:19000),intensity_mat(train,50:19000));
    ylim([0 max(intensity_mat(train,50:19000))*1.2])
               hold on;
           for spike = 1:length(temp{train})
                plot([temp{train}(spike) temp{train}(spike)], [0 max(intensity_mat(train,50:19000))/15])
           end
           hold off;

    subplot(2,1,2);
    plot(time(50:19000),deriv_mat(train,50:19000));
    temp1 = num2str(train);
               pngfilename = ['train_' temp1 '.png'];

           saveas(gcf, pngfilename);
end
