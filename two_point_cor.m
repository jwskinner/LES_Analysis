function [outputArg1,outputArg2] = two_point_cor(file, folder, nam, output, i)

%   Computes the two point correlation function of a 2D lattice
%   of a fixed width and height
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


fname=strcat(folder,file);

%% ========================================================================
% read all relevant 3-D variables

u=ncread(fname,'U'); u=0.5*(u(1:end-1,:,:,:)+u(2:end,:,:,:)); %C->A grid
v=ncread(fname,'V'); v=0.5*(v(:,1:end-1,:,:)+v(:,2:end,:,:)); 
w=ncread(fname,'W'); w=0.5*(w(:,:,1:end-1,:)+w(:,:,2:end,:));

ph =ncread(fname,'PH' );                                                   % geopotential perturbation [m2/s2]
phb=ncread(fname,'PHB');                                                   % base geopotential [m2/s2)
tke=ncread(fname,'TKE');                                                   % TKE [m2/s2]
p  =ncread(fname,'P'  );                                                   % pressure perturbation [Pa]
pb =ncread(fname,'PB' );                                                   % base pressure [Pa]
th =ncread(fname,'T'  )+nam.T0;                                            % Potential temperature [K]


%% 

%% ========================================================================
% initial calculations

s=size(w);n=s(1); m=s(2); l=s(3); nm=n*m;

HS=mean(reshape(ph+phb,nm,l+1));                                           % Jack added +1 here to fix reshape [check] 
H=0.5*(HS(:,1:end-1)+HS(:,2:end));                                         % Jack flipped the indexing to get the half values here 

ZS=HS./nam.g;                                                              % height at w-levels [m]
Z=H./nam.g';                                                               % height at mass-levels [m]
p=p+pb;                                                                    % pressure
exn=(p/nam.P0).^(nam.R/nam.cp);                                            % exner function                                                        

TH=mean(reshape(th,nm,l));
t=th.*exn;                                                                 % temperature

%% ========================================================================
% Correlation function





end