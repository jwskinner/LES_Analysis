% J.W.Skinner -- 14/03/2023
% This code loads wind data from a Weather Research and Forecasting (WRF)
% model output file and formats it for paraview for 3D visualisation
% There is now a time loop for bringing in the time evolving data 

folder = "./data/small_domain/CP_OUT/"; 

% Returns a list of all files in the folder
files_all = dir(strcat(folder, 'wrfout*')); 
file = files_all(1).name; % take first file in the series to setup the global parameters
fname=strcat(folder,file);

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
nam.dx=double(ncreadatt(fname,'/','DX'));                                  % [m]
nam.dy=double(ncreadatt(fname,'/','DY'));                                  % [m]

time = 3; 
qt_out = zeros(nam.nx-1, nam.ny, nam.levs, time);
thil_out = zeros(nam.nx-1, nam.ny, nam.levs, time);
u_out = zeros(nam.nx-1, nam.ny, nam.levs, time);
time_out = zeros(time, 1); 

for i = 1:3 %length(files_all)

fname=strcat(folder,files_all(i).name)

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

%% Setup the length grid rather than using grid points 
% Simple linear grid
x = 1:nam.nx; 
y = 1:nam.ny;  

% Convert grid coordinates to meters
lx = ((x-1)*nam.dx)/1000; % Lx in [km] (subtract 1 to account for 0-based indexing)
ly = ((y-1)*nam.dy)/1000; % Ly in [km] (subtract 1 to account for 0-based indexing)

% write the variables out of the time loop 
u_out(:,:,:,i) = u; 
thil_out(:,:,:,i) = thil; 
qt_out(:,:,:,i) = qt; 
time_out(i) = (i - 1)*nam.dt; 

end

%% Create the netcdf for paraview of the variables above
%create the netcdf file
p_out = './data/paraview_variables_cp.nc'; % Put the paraview variables in data so they don't sync to git
ncid = netcdf.create(p_out,'NC_WRITE'); 

% define dimensions
nx_dimid = netcdf.defDim(ncid,'nx',size(u_out, 1));
ny_dimid = netcdf.defDim(ncid,'ny',size(u_out, 2));
levs_dimid = netcdf.defDim(ncid,'levs',size(u_out, 3));
time_dimid = netcdf.defDim(ncid,'time',size(u_out, 4));

% define variables
u_varid = netcdf.defVar(ncid,'u','double',[nx_dimid,ny_dimid,levs_dimid, time_dimid]);
thil_varid = netcdf.defVar(ncid,'thil','double',[nx_dimid,ny_dimid,levs_dimid, time_dimid]);
qt_varid = netcdf.defVar(ncid,'qt','double',[nx_dimid,ny_dimid,levs_dimid, time_dimid]);

% add attribute metadata
netcdf.putAtt(ncid,u_varid,'units','m/s');
netcdf.putAtt(ncid,thil_varid,'units','K');
netcdf.putAtt(ncid,qt_varid,'units','kg/kg');

% Define time variable and write time out
time_varid = netcdf.defVar(ncid,'time','double',time_dimid);
netcdf.putAtt(ncid, time_varid, 'units', 'hours');

%  end definitions and close file
netcdf.endDef(ncid);

% put variables in the file 
netcdf.putVar(ncid, u_varid, u_out);
netcdf.putVar(ncid, thil_varid, thil_out);
netcdf.putVar(ncid, qt_varid, qt_out);
netcdf.putVar(ncid, time_varid, time_out);

% close netcdf file
netcdf.close(ncid);