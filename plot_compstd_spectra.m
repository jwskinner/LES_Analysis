% Script for plotting the compensated 2D kinetic energy spectrum from
% Matheou and Teixiera 2019 figure 14
%
% Calls in the function: comp_1D_spectra.m to compute spectrum of data and
% plots the spectrum with length scale on xaxis and normalised by k^-{5/3}
% on the yaxis to isolate length scales with kolmogorov turbulence. 
%
% J.W.Skinner
% last modified: 7/03/2023

clear all

index = 10; % file index to plot
z = 1; % vertical level to plot

folders = ["./data/large_domain/CP_OUT/", "./data/small_domain/NOCP_OUT/"];


% Import the file
files_all = dir(strcat(folders(1), 'wrfout*'));
file = files_all(index).name;
fname=strcat(folders(1),file);
fname_1=strcat(folders(2),file);
data = ncread(fname,'U');

% Benchmarking test of isotropic turbulence for testing the KE spectrum
% function
%  b_data = "./benchmark_data/rayleightaylor-turbulence.dat"
%  data = h5read(b_data, '/Velocity_0001'); % Use the h5read function to import the dataset
%  data = data(1, :, :, :);

[E_avg, kx, ky, nbins] = comp_1D_spectra(fname, z, data);


% Convert the wavenumbers to length scales
lx = (2*pi./kx).*1000; % length scales in x direction [m]
ly = (2*pi./ky).*1000; % length scales in y direction [m]

% Compensate the energy spectrum by dividing by k^(5/3)
k_comp = kx.^(5/3);
E_comp = E_avg./k_comp';

%% Plot the kinetic energy spectrum
%loglog((1:nbins)*min(kx(2), ky(2)), E_avg, 'linewidth', 2); hold on
loglog(kx, E_comp, 'linewidth', 2); hold on
loglog(kx, 1e10.*kx.^(-5/3), '--k', 'linewidth', 1)
xlabel('Length Scale l [m]')
ylabel('k_x^{5/3} E_u')
title('2D Kinetic Energy Spectrum')
legend('CP', 'NOCP', 'k^{-5/3} scaling','k^{-9/3} scaling' )
%set(gca, 'XDir','reverse')
%xlim([1, 10^4])