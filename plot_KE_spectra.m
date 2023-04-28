% Script for plotting a simple 2D kinetic energy spectrum from the simulations
% Uses the function KE_2D_spectra to compute kinetic energy spectrum. 
% J.W.Skinner
% last modified: 7/03/2023
clear all 

index = 97; % file number to plot 
z = 1; % vertical level to plot 

%folders = ["./data/large_domain/CP_OUT/", "./data/large_domain/NOCP_OUT/"];
folders = ["/scratch/05999/mkurowsk/GATE_NOCP_CONSTFLX/", "/scratch/05999/mkurowsk/GATE_CP_CONSTFLX/"]

% Import the file 
files_all = dir(strcat(folders(1), 'wrfout*'));

% Loop over all the files for creating a movie
% for i = 1:120
% index = i 

file = files_all(index).name;
fname_1=strcat(folders(1),file)
fname_2=strcat(folders(2),file);

% Compute the 2D KE Spectra for each file 
[E_avg, kx, ky, nbins] = KE_spectra(fname_1, z);
[E_avg_2, kx, ky, nbins] = KE_spectra(fname_2, z);

%% Plot the kinetic energy spectrum
f = gcf
loglog((1:nbins)*min(kx(2), ky(2)), E_avg, 'linewidth', 2); hold on
loglog((1:nbins)*min(kx(2), ky(2)), E_avg_2, 'linewidth', 2); hold on
loglog(kx, 1e12*kx.^(-5/3), '--k', 'linewidth', 1)
loglog(kx, 1e14*kx.^(-9/3), ':k', 'linewidth', 1)
xlabel('Wavenumber k')
ylabel('Kinetic Energy Spectrum')
title(strcat('t = ', num2str(i*0.5), ' hours'))
legend('NOCP', 'CP', 'k^{-5/3} scaling', 'k^{-9/3} scaling'); hold off 
xlim([1, 10^4])

%% Turn on export frames to gif
% exportgraphics(f,strcat('./plots/KE_Spec/', 'GATE', '.gif'),'Resolution',150, 'Append',true)
% 
% end 