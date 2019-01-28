%--------------------------------------------------------------------------
%   Returns CV score for one penalty coefficient
%
%   pen_coef   :  penalty coefficient
%   spike_time  :  time of events
%   nbasis     :  number of basis functions
%   basis_2nd  :  2nd derivative matrix of basis functions evluated along
%                 time
%   basis_value:  evaluation of basis functions along time
%
%e
%--------------------------------------------------------------------------



function cvscore = cvcheck(pen_coef,spike_time,nbasis,basis_2nd,basis_value)
  %retriving globa variables
  global time;           %time axis
  global bspline_coef;   %starting points for optimization, can be changed
  global options;        %options for fminsearch

  %pre-allocate space for estimated coefficients and intensity for each
  %leave-one-out data
  coef_mat = zeros(length(spike_time), nbasis);
  intensity_mat = zeros(length(spike_time),length(time));
%--------------------------------------------------------------------------
%calculate cross validation score;

  %Estimate intensity for each leave-one-out point process
  for ii = 1:length(spike_time)
    leaveone_spike = spike_time;
    leaveone_spike(ii) = [];
    coef_mat(ii,:) = fminsearch(@(coef) penlik(coef,leaveone_spike,pen_coef,nbasis,basis_2nd,basis_value), bspline_coef, options);
    intensity_mat(ii,:) = exp(coef_mat(ii,:) * transpose(basis_value));
        fprintf('Just finished leaveoneout #%d\n', ii);
  end

  crossval = zeros(1,length(spike_time));

  for jj = 1:length(spike_time)
      if (jj == 1)
          %estimated log intensity at removal spike
          temp1 = log(intensity_mat(jj,spike_time(jj)));
          %
          if(spike_time(jj) - time(jj) <2 )
              temp2 = intensity_mat(jj,spike_time(jj));
          else
          time1 = 1;
          time2 = spike_time(jj);
          npoints = length(time1:1:time2);
          pointwts = ones(1,npoints);
          pointwts(2:2:npoints) = 4;
          pointwts(3:2:npoints-1) = 2;
          pointwts(npoints) = 1;
          temp2 = sum(times(intensity_mat(jj,1:spike_time(jj)),pointwts))/3;

          end
           crossval(jj) = temp1 - temp2;
      else
            temp1 = log(intensity_mat(jj,spike_time(jj)));
            if (spike_time(jj) - spike_time(jj-1)<2)
                 temp2 = intensity_mat(jj,spike_time(jj));
            else

             %
             time1 = spike_time(jj-1);
             time2 = spike_time(jj);
              npoints = length(time1:1:time2);
              pointwts = ones(1,npoints);
              pointwts(2:2:npoints) = 4;
              pointwts(3:2:npoints-1) = 2;
              pointwts(npoints) = 1;
              temp2 = sum(times(intensity_mat(jj,spike_time(jj-1):spike_time(jj)),pointwts))/3;

            end
              crossval(jj) = temp1 - temp2;
       end

       fprintf('Just finished leaveoneout #%d\n', jj);
  end

cvscore = -sum(crossval);

end
