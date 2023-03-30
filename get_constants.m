function [nam] = get_constants(file)

fname = file; 

% Setup physical and numerical constants
nam.R=287.04;                                                              % [J/kg/K]
nam.cp=1004.67;                                                            % [J/kg/K]
nam.g=9.81;                                                                % [m/s2]
nam.Ll=2.50e6;                                                             % latent heat of evaporation (vapor:liquid) at 0C [J/kg]
nam.Li=2.83e6;                                                             % latent heat of sublimation (vapor:solid) at 0C [J/kg]
nam.T0=300;                                                                % ncread(fname,'T00'); % base state temperature [K]
nam.P0=1.e5;                                                               % ncread(fname,'P00'); % base state pressure [Pa]
nam.dx=double(ncreadatt(fname,'/','DX'));                                   % [m]
nam.dy=double(ncreadatt(fname,'/','DY'));                                   % [m]
nam.dt = 0.5;                                                              % Output frequency [hours]
nam.levs = size(ncread(fname,'U'), 3);                                      % Number of vertical levels in the simulation
nam.nx = size(ncread(fname,'U'), 1);                                        % Number of x grid points in simulation
nam.ny = size(ncread(fname,'U'), 2);   

end