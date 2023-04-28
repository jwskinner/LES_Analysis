function [out] = two_point_cor(variable, nam)

%   J.W. Skinner function
%   Computes the two point correlation function of a 2D lattice in the
%   zonal (x) direction for LES for a fixed width and height.
%
%   x - list of x coordinates of points
%   y - list of y coordinates of points
%   dr - binning distance for successive circles
%   blksize - optional, default 1000, number of points to be considered
%   in a single step.
%   verbose - if true, will print which step is currently processed
%
%   coorfun - two point correlation function
%   r - r-coordinates for the coordfun
%   rw - number of particles that participated in particular r value
%   corrfun computation. Low rw means corrfun is unreliable at that r.

data = variable.Data; % Extract data from the LES variable structure
% Define the number of data points
N = length(data);

% Compute the two-point correlation in the zonal direction
correlation = zeros(N, 1);
for i = 1:N
    for j = 1:N
        correlation(i) = correlation(i) + data(j)*data(j+i);
    end
    correlation(i) = correlation(i)/N;
end

out = correlation; 
% % Plot the results
% figure()
% plot(correlation)
% xlabel('Lag Distance')
% ylabel('Two-Point Correlation')





end