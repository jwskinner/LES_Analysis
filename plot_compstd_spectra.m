% Script for plotting the compensated 2D kinetic energy spectrum from
% Matheou and Teixiera 2019 figure 14
%
%
% J.W.Skinner
% last modified: 7/03/2023

clear all

index = 10; % file index to plot
z = 1; % vertical level to plot

folders = ["./data/large_domain/CP_OUT/", "./data/large_domain/NOCP_OUT/"];

% Import the file
files_all = dir(strcat(folders(1), 'wrfout*'));
file = files_all(index).name;
fname=strcat(folders(1),file);
fname_1=strcat(folders(2),file);

[E_avg, kx, ky, nbins] = KE_2D_spectra(fname, z);
[E_avg_1, kx, ky, nbins] = KE_2D_spectra(fname_1, z);

% % Setup the simulation grid from the LES output
% nx = size(ncread(fname,'U'), 1);                % number of grid points in x direction
% ny = size(ncread(fname,'U'), 2);                % number of grid points in y direction
% dx = double(ncreadatt(fname,'/','DX'));         % grid spacing in x direction [m]
% dy = double(ncreadatt(fname,'/','DY'));         % grid spacing in y direction [m]
% Lx = (nx*dx)/1000;                              % length of domain in x direction [km]
% Ly = (ny*dy)/1000;                              % length of domain in y direction [km]
% x = linspace(0, Lx-dx, nx);                     % x coordinate vector
% y = linspace(0, Ly-dy, ny);                     % y coordinate vector
% [X, Y] = meshgrid(x, y);                        % meshgrid of x and y coordinates
% 
% % Read the velocity field
% [u, v, w] = loadNetCDF(fname, 'U', 'V', 'W');
% u = u.Data; v = v.Data; w = w.Data;
% 
% u=0.5*(u(1:end-1,:,:,:)+u(2:end,:,:,:)); %C->A grid
% v=0.5*(v(:,1:end-1,:,:)+v(:,2:end,:,:));
% w=0.5*(w(:,:,1:end-1,:)+w(:,:,2:end,:));
% 
% % restrict to height, z
% u = u(:,:,z); v = v(:, :, z);
% 
% % Compute the Fourier transform of the velocity field
% uhat = fft2(u);
% vhat = fft2(v);
% 
% % Compute the kinetic energy spectrum
% kx = [0:(nx/2-1) (-nx/2):-1];  % wavenumbers in x direction
% ky = [0:(ny/2-1) (-ny/2):-1];  % wavenumbers in y direction
% [KX, KY] = meshgrid(kx, ky);
% 
% K = sqrt(KX.^2 + KY.^2);  % magnitude of the wavevector
% E = 0.5*(abs(uhat).^2 + abs(vhat).^2);  % kinetic energy spectrum
% 
% % Average energy over annular bands in Fourier space
% nbins = size(kx, 2);  % number of annular bins
% E_avg = zeros(nbins, 1);  % average energy in each bin
% count = zeros(nbins, 1);  % number of wavevectors in each bin
% for i = 1:nbins
%     kmin = (i-1)*min(kx(2), ky(2));
%     kmax = i*min(kx(2), ky(2));
%     ind = find((K >= kmin) & (K < kmax));
%     count(i) = length(ind);
%     E_avg(i) = sum(E(ind))/count(i);
% end

% Convert the wavenumbers to length scales
lx = (2*pi./kx).*1000; % length scales in x direction [m]
ly = (2*pi./ky).*1000; % length scales in y direction [m]

E_comp = E_avg ./ (kx'.^(-5/3))
E_comp_1 = E_avg_1 ./ (kx'.^(-5/3))

%% Plot the kinetic energy spectrum
%loglog((1:nbins)*min(kx(2), ky(2)), E_avg, 'linewidth', 2); hold on
loglog(lx, E_comp, 'linewidth', 2); hold on
loglog(lx, E_comp_1, 'linewidth', 2); hold on
xlabel('Length Scale l [m]')
ylabel('k_x^{5/3} E_u')
title('2D Kinetic Energy Spectrum')
legend('CP', 'NOCP', 'k^{-5/3} scaling','k^{-9/3} scaling' )
set(gca, 'XDir','reverse')
xlim([1, 10^4])