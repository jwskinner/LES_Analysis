% J.W.Skinner -- 10/03/2023
%
% This code loads wind data from a Weather Research and Forecasting (WRF)
% model output file and computes the kinetic energy spectrum of the wind
% field at a specified height index.
% Specifically, the code loads the U-component wind data from the file,
% computes its 2D Fourier transform in the x and y directions, and then
% computes the 1D kinetic energy spectrum by averaging the squared
% magnitude of the Fourier coefficients along the y direction.
% The resulting TKE spectrum is plotted on a log-log scale versus
% wavenumber in units of 1/m, along with a line representing the kolmogorov
% k^-5/3 scaling often observed in 3D isotropic turbulence.
%
% [This is still a rough/unfinished script!]

loc_cp = "./data/GATE_CP_int_domain/wrfout_d01_2020-01-03_00:30:00";
% loc_cp = "./data/large_domain/CP_OUT/wrfout_d01_2020-01-01_18:00:00";
% loc_cp = "./data/small_domain/CP_OUT/wrfout_d01_0001-01-01_15:00:00";

ucp = ncread(loc_cp, 'U'); 
vcp = ncread(loc_cp, 'V'); 

U_CP=0.5*(ucp(1:end-1,:,:,:)+ucp(2:end,:,:,:)); %C->A grid
V_CP=0.5*(vcp(:,1:end-1,:,:)+vcp(:,2:end,:,:)); 

% Define the fixed height index for the TKE spectrum
z_idx = 1; % First layer above the boundary (correponds to ~60m)

% Define the grid spacing
dx = double(ncreadatt(loc_cp,'/','DX')) / 1000; % grid spacing in km 

nx = size(ucp, 1); 

% Compute Spectrum of CP and NOCP solutions
[TKE_spectrum_CP, KX] = compute_TKE_spectrum(U_CP, V_CP, z_idx, dx); 

%% Plot the TKE spectrum versus wavenumber in units of 1/m or Length in m
%figure;
loglog(pi./KX, TKE_spectrum_CP, 'LineWidth',1.5); hold on
loglog(2*pi./KX(:), 1e7*(2*pi./KX(:)).^(5/3), 'k--', 'LineWidth',1.2); % Add line for k^-5/3 scaling
xlabel('Length, [km]');
set(gca, 'XDir','reverse')
ylabel('E (k)');
% title(sprintf('E_{uu} spectrum at z=%.1f m', z_idx*dz));
legend('CP', 'k^{-5/3}')
xlim([0.01, 3*10^3])
ylim([100, 10^10])

function [Spectrum, KX] = compute_1D_spectrum(U, z_idx) % Still in progress

% Compute the size of U
[nx, ny] = size(U(:,:,z_idx));

% Compute the 2D Fourier transform of U in the x and y directions
% fft2 function returns a normalized 2D Fourier transform, where each 
% coefficient is divided by the total number of grid points (nx * ny) 
% before being centered with the fftshift function
U_spectral_2d = fftshift(fft2(U(:,:,z_idx))); 

% Define the wavenumber array in the x direction
KX = (-nx/2:nx/2-1) / (nx*dx);

% Compute the 1D kinetic energy spectrum
U_spectral_mag2 = abs(U_spectral_2d).^2;
U_spectral_mag2_mean = mean(U_spectral_mag2, 2);
TKE_spectrum = 0.5 * U_spectral_mag2_mean';

% Normalise the spectrum if desired 
% TKE_spectrum = TKE_spectrum / sum(TKE_spectrum);

end

function [TKE_spectrum, KX] = compute_TKE_spectrum(U, V, z_idx, dx)

% Compute the size of U and V
[nx, ny] = size(U(:,:,z_idx))

size(V(:,:,z_idx))

% Compute the 2D Fourier transform of U and V in the x and y directions
% fft2 function returns a normalized 2D Fourier transform, where each 
% coefficient is divided by the total number of grid points (nx * ny) 
% before being centered with the fftshift function
U_spectral_2d = fftshift(fft2(U(:,:,z_idx))); 
V_spectral_2d = fftshift(fft2(V(:,:,z_idx))); 

% Define the wavenumber array in the x direction
KX = 2*pi*(-nx/2:nx/2-1) / (nx*dx);

% Compute the 2D kinetic energy spectrum
U_spectral_mag2 = abs(U_spectral_2d).^2;
V_spectral_mag2 = abs(V_spectral_2d).^2;

TKE_spectrum_2D = 0.5 * (U_spectral_mag2 + V_spectral_mag2);

% Compute the 1D kinetic energy spectrum by averaging over y direction
TKE_spectrum_1D = mean(TKE_spectrum_2D, 2);

% Return the 1D TKE spectrum
TKE_spectrum = TKE_spectrum_1D';

% Normalise the spectrum if desired 
% TKE_spectrum = TKE_spectrum / sum(TKE_spectrum);

end
