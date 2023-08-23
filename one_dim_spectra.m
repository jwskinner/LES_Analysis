% J.W.Skinner -- 14/03/2023
% This code loads wind data from a Weather Research and Forecasting (WRF)
% model output file and computes the kinetic energy spectrum of the wind
% field at a specified height index.
% Specifically, the code loads the U-component wind data from the file,
% computes its 2D Fourier transform in the x and y directions, and then
% computes the 1D kinetic energy spectrum by averaging the squared
% magnitude of the Fourier coefficients along the y direction.
% The resulting TKE spectrum is plotted on a log-log scale versus
% wavenumber in units of 1/m, along with a line representing the k^-5/3
% scaling often observed in atmospheric turbulence.

fname = "./data/small_domain/NOCP_OUT/wrfout_d01_0001-01-01_15:30:00"; 

% Define the fixed height index for the TKE spectrum
z_idx = 2;

% Setup physical and numerical constants
nam.R=287.04;                                                              % [J/kg/K]
nam.cp=1004.67;                                                            % [J/kg/K]
nam.g=9.81;                                                                % [m/s2]
nam.Ll=2.50e6;                                                             % latent heat of evaporation (vapor:liquid) at 0C [J/kg]
nam.Li=2.83e6;                                                             % latent heat of sublimation (vapor:solid) at 0C [J/kg]
nam.T0=300;                                                                % ncread(fname,'T00'); % base state temperature [K]
nam.P0=1.e5;                                                               % ncread(fname,'P00'); % base state pressure [Pa]
nam.dt = 0.5;                                                              % Output frequency [hours]
nam.levs = size(ncread(fname,'U'), 3);                                     % Number of vertical levels in the simulation
nam.nx = size(ncread(fname,'U'), 1);                                       % Number of x grid points in simulation
nam.ny = size(ncread(fname,'U'), 2);                                       % Number of y grid points in simulation

dx=double(ncreadatt(fname,'/','DX'));                                      % [m]
dy=double(ncreadatt(fname,'/','DY'));                                      % [m]

%% Import all the variables needed
u = ncread(fname, 'U'); 
v = ncread(fname, 'V');
w = ncread(fname, 'W');

u=0.5*(u(1:end-1,:,:,:)+u(2:end,:,:,:)); %C->A grid
v=0.5*(v(:,1:end-1,:,:)+v(:,2:end,:,:));

ph =ncread(fname,'PH' );                                                   % geopotential perturbation [m2/s2]
phb=ncread(fname,'PHB');                                                   % base geopotential [m2/s2)
p  =ncread(fname,'P'  );                                                   % pressure perturbation [Pa]
pb =ncread(fname,'PB' );                                                   % base pressure [Pa]
th =ncread(fname,'T'  )+nam.T0;                                             % Potential temperature [K]

qv=ncread(fname,'QVAPOR');                                                 % Water vapor mixing ratio [kg/kg]
qc=ncread(fname,'QCLOUD');                                                 % Cloud water mixing ratio [kg/kg]
qi=ncread(fname,'QICE');                                                   % Ice mixing ratio [kg/kg]

qt=qv+qc+qi;                                                               % total water mixing ratio [kg/kg] (no precipitating elements)

s=size(u);
n=s(1); 
m=s(2); 
l=s(3); 
nm=n*m;

HS=mean(reshape(ph+phb,nm,l+1));                                           
H=0.5*(HS(:,1:end-1)+HS(:,2:end));                                         

ZS=HS./nam.g;                                                              % height at w-levels [m]
Z=H./nam.g';                                                               % height at mass-levels [m]
p=p+pb;                                                                    % pressure
exn=(p/nam.P0).^(nam.R/nam.cp);                                            % exner function
qt=qv+qc+qi;                                                               % total water mixing ratio [kg/kg] (no precipitating elements)

TH=mean(reshape(th,nm,l));
t=th.*exn;                                                                 % temperature
tv=t.*(1+0.608*qv);                                                        % virtual temperature, bouyancy is tv - ql (eq. 1 Marcin)
thil=th-(nam.Ll*qc+nam.Li*qi)./(nam.cp*exn);                               % liquid water potential temperature

Nx = size(u, 1); % number of grid points along the x-direction
k = 2*pi*(-Nx/2:Nx/2-1)/(Nx*dx); % wavenumber vector

len = 2*pi ./ abs(k); %converts wavenumbers to length 

u = u(:,:, z_idx); v = v(:,:,z_idx); w = w(:,:,z_idx);
thil = thil(:, :, z_idx); qt = qt(:,:,z_idx);

%% Compute the one-dimensional spectra along the x-direction
Fu = fftshift(fft(u, [], 1), 1); % Fourier transform of u along the x-direction
Fw = fftshift(fft(w, [], 1), 1); % Fourier transform of w along the x-direction
Ful = fftshift(fft(thil, [], 1), 1); % Fourier transform of ul along the x-direction
Fqt = fftshift(fft(qt, [], 1), 1); % Fourier transform of qt along the x-direction

Pu = mean(abs(Fu).^2,2); % power spectrum of u
Pw = mean(abs(Fw).^2,2); % power spectrum of w
Pthil = mean(abs(Ful).^2,2); % power spectrum of ul
Pqt = mean(abs(Fqt).^2,2); % power spectrum of qt

% Compute the normalized spectral function F(k)
Fu_norm = Pu / var(u, [], 1)'; % normalized spectral function of u
Fw_norm = Pw / var(w, [], 1)'; % normalized spectral function of w
Fthil_norm = Pthil / var(thil, [], 1)'; % normalized spectral function of thil
Fqt_norm = Pqt / var(qt, [], 1)'; % normalized spectral function of qt

% Compute the premultiplied spectrum k*F(k)
kFu_norm = k.'.*Fu_norm; % premultiplied spectrum of u
kFw_norm = k.'.*Fw_norm; % premultiplied spectrum of w
kFthil_norm = k.'.*Fthil_norm; % premultiplied spectrum of thil
kFqt_norm = k.'.*Fqt_norm; % premultiplied spectrum of qt

%% Plot the spectra
figure(1);
subplot(1, 4, 1);
loglog(len, Fu_norm); hold on
loglog(len, 10^-2*len.^(5/3), 'k--');
xlabel('l_x (m)');
ylabel('\Phi_{uu} u');
xlim([10^0, 10^6])
set(gca, 'XDir', 'reverse')

subplot(1, 4, 2);
loglog(len, Fthil_norm); hold on 
loglog(len, 10^-2*len.^(5/3), 'k--');
xlabel('l_x (m)');
ylabel('\Phi_{\theta\theta} \theta_l');
xlim([10^0, 10^6])
set(gca, 'XDir', 'reverse')

subplot(1, 4, 3);
loglog(len, Fw_norm); hold on 
loglog(len, 10^-2*len.^(5/3), 'k--');
xlabel('l_x (m)');
ylabel('\Phi_{ww} w');
xlim([10^0, 10^6])
set(gca, 'XDir', 'reverse')

subplot(1, 4, 4);
loglog(len, Fqt_norm); hold on 
loglog(len, 10^-2*len.^(5/3), 'k--');
xlabel('l_x (m)');
ylabel('\Phi_{qq} q_t');
xlim([10^0, 10^6])
set(gca, 'XDir', 'reverse')

%% Plot the spectra
figure(2);
subplot(1, 4, 1);
loglog(len, kFu_norm); hold on
xlabel('l_x (m)');
ylabel('\Phi u');
xlim([10^0, 10^6])
set(gca, 'XDir', 'reverse')

subplot(1, 4, 2);
loglog(len, kFthil_norm); hold on
xlabel('l_x (m)');
ylabel('\Phi \theta_l');
xlim([10^0, 10^6])
set(gca, 'XDir', 'reverse')

subplot(1, 4, 3);
loglog(len, kFw_norm)
xlabel('l_x (m)');
ylabel('\Phi w');
xlim([10^0, 10^6])
set(gca, 'XDir', 'reverse')

subplot(1, 4, 4);
loglog(len, kFqt_norm)
xlabel('l_x (m)');
ylabel('\Phi q_t');
xlim([10^0, 10^6])
set(gca, 'XDir', 'reverse')