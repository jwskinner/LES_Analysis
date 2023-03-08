% Script for plotting a simple 2D kinetic energy spectrum from the simulations
% Uses the function KE_2D_spectra to compute kinetic energy spectrum. 
% J.W.Skinner
% last modified: 7/03/2023
clear all 

index = 10; % file index to plot 
z = 1; % vertical level to plot 

folders = ["./data/large_domain/CP_OUT/", "./data/large_domain/NOCP_OUT/"];

% Import the file 
files_all = dir(strcat(folders(1), 'wrfout*'));
file = files_all(index).name;

fname_1=strcat(folders(1),file);
fname_2=strcat(folders(2),file);

% Compute the 2D KE Spectra for each file 
[E_avg, kx, ky, nbins] = KE_2D_spectra(fname_1, z);
[E_avg_2, kx, ky, nbins] = KE_2D_spectra(fname_2, z);

%% Plot the kinetic energy spectrum
loglog((1:nbins)*min(kx(2), ky(2)), E_avg, 'linewidth', 2); hold on
loglog((1:nbins)*min(kx(2), ky(2)), E_avg_2, 'linewidth', 2); hold on
loglog(kx, 1e12*kx.^(-5/3), '--k', 'linewidth', 1)
loglog(kx, 1e14*kx.^(-9/3), ':k', 'linewidth', 1)
xlabel('Wavenumber k')
ylabel('Kinetic Energy Spectrum')
title('2D Kinetic Energy Spectrum')
legend('NOCP', 'CP', 'k^{-5/3} scaling', 'k^{-9/3} scaling')
xlim([1, 10^4])