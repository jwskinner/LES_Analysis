% J.W.Skinner modification of M. Chinita's (2/14/2023) script for WRF data
% This code is a template for the auto-correlation function. It assumes
% the following variables have been loaded:
% z (vertical grid), u (3D zonal wind component), v, and w 

clear all

% nz, nx, and ny are the number of grid points in the 3D domain, e.g., 
% from Chinita et al 2022, run B = 150 x 256 x 256 grid points, so nz =
% 150, nx = 256, ny = 256 grid points.

fname = "./data/GATE_NOCP_int_domain/wrfout_d01_2020-01-03_00:00:00"; 

% Define the fixed height index for the TKE spectrum
% Vertical level at the top of the SBL
% here, using vertical TKE profiles the the top of the SBL is defined as 
% the height at which TKE reaches its minimum value (z = 12km from the 
% vertical profile plots).
hz = 109; % Corresponds to height of 11.95 [km] above the surface.

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

ph =ncread(fname,'PH' );                                                   % geopotential perturbation [m2/s2]
phb=ncread(fname,'PHB');                                                   % base geopotential [m2/s2)
p  =ncread(fname,'P'  );                                                   % pressure perturbation [Pa]
pb =ncread(fname,'PB' );                                                   % base pressure [Pa]
th =ncread(fname,'T'  )+nam.T0;                                            % Potential temperature [K]

s=size(th);
n=s(1); 
m=s(2); 
l=s(3); 
nm=n*m;

HS=mean(reshape(ph+phb,nm,l+1));                                           
H=0.5*(HS(:,1:end-1)+HS(:,2:end));                                         

ZS=HS./nam.g;                                                              % height at w-levels [m]
Z=H./nam.g';                                                               % height at mass-levels [m]

nx = nam.nx; 
ny = nam.ny; 
nz = nam.levs; % yes yes, I know this is a really messy way to do this but I will fix it later

% Instantaneous u-,v-, and w-components 
u2=zeros(nx,ny,nz);
v2=zeros(nx,ny,nz);
w2=zeros(nx,ny,nz);

% WRF model is in grid u = [nx + 1, ny, nz]
u2=0.5*(u(1:end-1,:,:,:)+u(2:end,:,:,:)); %C->A grid
v2=0.5*(v(:,1:end-1,:,:)+v(:,2:end,:,:));

% % The instantaneous fields of UConn's LES model are in a Arakawa C (staggered) grid, 
% % e.g., u = nz x nx+1 x ny, so here it converts to nz x nx x ny
% for i=1:nx
%    for j=1:ny
%        for k=1:nz
%            u2(i,j,k)=0.5*(u(k,i,j)+u(k,i+1,j));
%            v2(i,j,k)=0.5*(v(k,i,j)+v(k,i,j+1));
%            w2(i,j,k)=0.5*(w(k,i,j)+w(k+1,i,j));
%        end
%    end 
% end

%% Compute turbulent part of the flow (i.e., fluctuations)
u2_mean= squeeze(mean(mean(u2)));
v2_mean= squeeze(mean(mean(v2)));
w2_mean= squeeze(mean(mean(w2)));

nx = size(u2,1); 
ny = size(u2,2); 
nz = size(u2,3); 

for i=1:nx
   for j=1:ny
       for k=1:nz
           uf(i,j,k)=u2(i,j,k)-u2_mean(k);
           vf(i,j,k)=v2(i,j,k)-v2_mean(k);
           wf(i,j,k)=w2(i,j,k)-w2_mean(k);
       end
   end 
end

% Calculate variance
uf2 = uf.*uf;
variance_uf2 = squeeze(mean(mean(uf2)));

vf2 = vf.*vf;
variance_vf2 = squeeze(mean(mean(vf2)));

wf2 = wf.*wf;
variance_wf2 = squeeze(mean(mean(wf2)));

% Two-point autocorrelation function using rx by ry grid points 
box_x_m = 5; % Box size in [km] 
box_y_m = 5; % box size in [km] 

rx = round(box_x_m * 1000 / dx); % convert box size in km to grid points
ry = round(box_y_m * 1000 / dx);

% Let's create a frame of zeros around the turbulent quantities
uf_framed = zeros(nx+rx*2,ny+rx*2,hz); 
uf_framed(rx+1:end-rx,ry+1:end-ry,:) = uf(:,:,1:hz); % center

% The domain is doubly periodic so let's replace the zeros with the values
% of the opposite size
% uf_framed(1:rx,ry+1:end-ry,:) = uf(end-rx+1:end,:,1:hz); % left side uses the right side of uf
% uf_framed(end-rx+1,ry+1:end-ry,:) = uf(1:rx,:,1:hz); % right side uses the left side of uf [fix matrix sizing here]
% uf_framed(rx+1:end-rx,1:ry,:) = uf(:,end-ry+1:end,1:hz); % top uses the bottom of uf
% uf_framed(rx+1:end-rx,end-ry+1,:) = uf(:,1:ry,1:hz); % top uses the bottom of uf

% There was an array size bug above so testing a slightly new version
uf_framed(1:rx,ry+1:end-ry,:) = uf(end-rx+1:end,:,1:hz); % left side uses the right side of uf
uf_framed(end-rx+1:end,ry+1:end-ry,:) = uf(1:rx,:,1:hz); % right side uses the left side of uf [fix matrix sizing here]
uf_framed(rx+1:end-rx,1:ry,:) = uf(:,end-ry+1:end,1:hz); % top uses the bottom of uf
uf_framed(rx+1:end-rx,end-ry+1:end,:) = uf(:,1:ry,1:hz); % top uses the bottom of uf

vecx = (-rx:rx); 
vecy = (-rx:rx);

for izz=1:hz
    iz = izz
    cc = 0;
    for i = 1:size(vecx,2)
        for j = 1:size(vecy,2)
            for ii=rx+1:nx+rx
                for jj = ry+1:ny+ry
                    cc = cc+1;
                    helpR(cc) = uf_framed(ii,jj,iz).*uf_framed(ii+vecx(i),jj+vecy(j),iz);
                end
            end
            Ruu(i,j) = sum(helpR)./cc;
            Ruu_norm(i,j) = Ruu(i,j)./variance_uf2(iz);
            cc = 0;
        end
    end
        
    Ruu_norm_z(izz,:,:) = Ruu_norm;
end

% Integral length scale;  dx is the grid resolution
for izz = 1:hz
    L1(izz) = sum(Ruu_norm_z(izz,:,rx+1).*dx);
end

%% Figure
hsbl = Z(hz); % normalize by PBL height [[jack compute this and get scalings right]
dz = mean(diff(Z(1:hz))); % get the mean dz in the specified height range 

[X,Y] = meshgrid((vecy.*dz)./hsbl,(vecx.*dz)./hsbl);
hfig=figure; 
[C,h]=contour(Y,X,squeeze(Ruu_norm_z(1,:,:)),'LineColor','k','ShowText','on'); hold on
clabel(C,h,'Interpreter','Latex','FontSize',12)
% If there is overturning motions, plot those using dashed lines
[C,h]=contour(Y,X,squeeze(Ruu_norm_z(1,:,:)),'--','LevelList',(-0.1),'LineColor','k','ShowText','on')
clabel(C,h,'Interpreter','Latex','FontSize',12)
xlim([-0.3 0.3]); ylim([-0.3 0.3]);
% xlim([-1 1]); ylim([-1 1]);
ylabel('$r_y$','Interpreter','Latex','Fontsize',18); 
xlabel('$r_x$','Interpreter','Latex','Fontsize',18); 
text(4.5,5,'$R_{11}$', 'Interpreter', 'latex','Fontsize',12);
title('z/h = 0.20','Interpreter', 'latex','Fontsize',12)
set(gca,'FontSize',12,'TickLabelInterpreter','Latex')
pbaspect([1 1 1])

%% plot the mean u field 
figure; 
imagesc(squeeze(mean(u, 3)))
colorbar;
caxis([-1, -0.4]);

% Repeat for vv and ww 