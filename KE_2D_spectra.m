%% Function for returning the 2D Kinetic energy spectrum 
% J.W.Skinner modified 7/03/2023

function [E_avg, kx, ky, nbins] = KE_2D_spectra(fname, z)

    % Setup the simulation grid from the LES output
    nx = size(ncread(fname,'U'), 1);                % number of grid points in x direction
    ny = size(ncread(fname,'U'), 2);                % number of grid points in y direction
    dx = double(ncreadatt(fname,'/','DX'));         % grid spacing in x direction [m]
    dy = double(ncreadatt(fname,'/','DY'));         % grid spacing in y direction [m]
    Lx = (nx*dx)/1000;                              % length of domain in x direction [km]
    Ly = (ny*dy)/1000;                              % length of domain in y direction [km]
    x = linspace(0, Lx-dx, nx);                     % x coordinate vector
    y = linspace(0, Ly-dy, ny);                     % y coordinate vector
    [X, Y] = meshgrid(x, y);                        % meshgrid of x and y coordinates

    % Read the velocity field
    [u, v, w] = loadNetCDF(fname, 'U', 'V', 'W');
    u = u.Data; v = v.Data; w = w.Data;

    u=0.5*(u(1:end-1,:,:,:)+u(2:end,:,:,:)); %C->A grid
    v=0.5*(v(:,1:end-1,:,:)+v(:,2:end,:,:));
    w=0.5*(w(:,:,1:end-1,:)+w(:,:,2:end,:));

    % restrict to height, z
    u = u(:,:,z); v = v(:, :, z);

    % Compute the Fourier transform of the velocity field
    uhat = fft2(u);
    vhat = fft2(v);

    % Compute the kinetic energy spectrum
%     kx = (2*pi/Lx)*[0:(nx/2-1) (-nx/2):-1];  % wavenumbers in x direction
%     ky = (2*pi/Ly)*[0:(ny/2-1) (-ny/2):-1];  % wavenumbers in y direction
    kx = [0:(nx/2-1) (-nx/2):-1];  % wavenumbers in x direction
    ky = [0:(ny/2-1) (-ny/2):-1];  % wavenumbers in y direction
    [KX, KY] = meshgrid(kx, ky);

    K = sqrt(KX.^2 + KY.^2);  % magnitude of the wavevector
    E = 0.5*(abs(uhat).^2 + abs(vhat).^2);  % kinetic energy spectrum

    % Average energy over annular bands in Fourier space
    nbins = size(kx, 2);  % number of annular bins
    E_avg = zeros(nbins, 1);  % average energy in each bin
    count = zeros(nbins, 1);  % number of wavevectors in each bin
    for i = 1:nbins
        kmin = (i-1)*min(kx(2), ky(2));
        kmax = i*min(kx(2), ky(2));
        ind = find((K >= kmin) & (K < kmax));
        count(i) = length(ind);
        E_avg(i) = sum(E(ind))/count(i);
    end
end