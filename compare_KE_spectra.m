% J.W.Skinner -- 10/03/2023
%
% Loads wind data from WRF
% Outputs kinetic energy spectrum of the wind field at a specified height.
% 
% Steps:
% 1) Load U-component wind data from the file. 
% 2) Compute vertical structure in [m] to help select height in [m]
% 3) Compute 2D Fourier transform in the x and y directions
% 4) Compute 1D kinetic energy spectrum by averaging the squared
% magnitude of the Fourier coefficients along the y direction.
% 5) plot TKE spectrum in log-log scale vs wavenumber in units of 1/m, 
% and line representing the k^-5/3 turbulence scaling. 

%% Load the data from files (two here for comparing CP and NOCP)
loc_cp = "./data/large_domain/CP_OUT/wrfout_d01_2020-01-02_06:00:00";
loc_nocp = "./data/large_domain/NOCP_OUT/wrfout_d01_2020-01-02_06:00:00";

% Define the grid spacing
dx = 100; % meters
dy = 100; % meters
dz = 100; % meters

ucp = ncread(loc_cp, 'U');
vcp = ncread(loc_cp, 'V');
U_CP=0.5*(ucp(1:end-1,:,:,:)+ucp(2:end,:,:,:)); %C->A grid
V_CP=0.5*(vcp(:,1:end-1,:,:)+vcp(:,2:end,:,:));

unocp = ncread(loc_nocp, 'U');
vnocp = ncread(loc_nocp, 'V');
U_NOCP=0.5*(unocp(1:end-1,:,:,:)+unocp(2:end,:,:,:)); %C->A grid
V_NOCP=0.5*(vnocp(:,1:end-1,:,:)+vnocp(:,2:end,:,:));

%% Define the fixed height index for the spectrum in [m]
height_m = 400;

% Get the vertical structure so we can select a height in [m]
ph =ncread(loc_cp,'PH' );                                                  % geopotential perturbation [m2/s2]
phb=ncread(loc_cp,'PHB');                                                  % base geopotential [m2/s2)
p  =ncread(loc_cp,'P'  );                                                  % pressure perturbation [Pa]
pb =ncread(loc_cp,'PB' );                                                  % base pressure [Pa]
th =ncread(loc_cp,'T'  )+300;                                              % Potential temperature [K]

s=size(phb);
n=s(1); m=s(2); l=s(3); nm=n*m;

HS=mean(reshape(ph+phb,nm,l));
ZS=HS./9.81;                                                               % height at w-levels [m]

% Find the index for the height by minimising the difference of the list
% and height_m
diffs = abs(ZS - height_m);
[~, z_idx] = min(diffs);

%% Compute Spectrum of CP and NOCP solutions calling function defined below 

[TKE_spectrum_CP, KX] = compute_TKE_spectrum(U_CP, V_CP, z_idx);
[TKE_spectrum_NOCP, KX] = compute_TKE_spectrum(U_NOCP, V_NOCP, z_idx);

%% Plot the TKE spectrum versus wavenumber in units of 1/m
figure;
loglog(KX, TKE_spectrum_CP, 'LineWidth',1.5); hold on
loglog(KX, TKE_spectrum_NOCP, 'LineWidth',1.5); hold on
loglog(KX(1020:end), 0.5*1e11*(KX(1020:end)).^(-5/3), 'k--', 'LineWidth',1.2); % Add line for k^-5/3 scaling
xlabel('Wavenumber, k');
ylabel('TKE (m^2/s^2)');
% title(sprintf('E_{uu} spectrum at z=%.1f m', z_idx*dz));
legend('CP', 'NOCP', 'k^{-5/3}')
xlim([1, 3*10^3])
ylim([100, 10^10])


%% Functions used are defined below [will split these into seperate files later]
function [Spectrum, KX] = compute_1D_spectrum(U, z_idx)

% Compute the size of U
[nx, ny] = size(U(:,:,z_idx));

% Compute the 2D Fourier transform of U in the x and y directions
% fft2 function returns a normalized 2D Fourier transform, where each
% coefficient is divided by the total number of grid points (nx * ny)
% before being centered with the fftshift function
U_spectral_2d = fftshift(fft2(U(:,:,z_idx)));

% Define the wavenumber array in the x direction
KX = (-nx/2:nx/2-1) %/ (nx*dx);

% Compute the 1D kinetic energy spectrum
U_spectral_mag2 = abs(U_spectral_2d).^2;
U_spectral_mag2_mean = mean(U_spectral_mag2, 2);
TKE_spectrum = 0.5 * U_spectral_mag2_mean';

% Normalise the spectrum if desired
% TKE_spectrum = TKE_spectrum / sum(TKE_spectrum);

end

function [TKE_spectrum, KX] = compute_TKE_spectrum(U, V, z_idx)

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
KX = (-nx/2:nx/2-1); %/ (nx*dx);

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
