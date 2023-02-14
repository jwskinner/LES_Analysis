% Prints out the fluxes from WRF for debugging the simulations. 

clear variables

txt = 'NOCP';

folder = './data/large_domain/CP_OUT/'
%folder = '/scratch/05999/mkurowsk/ocean_nocp/';

output = './plots/';

% Get a list of all files in the current directory
files_all = dir(strcat(folder, 'wrfout*'));

nam.R=287.04;                                                              % [J/kg/K]
nam.cp=1004.67;                                                            % [J/kg/K]
nam.g=9.81;                                                                % [m/s2]
nam.Ll=2.50e6;                                                             % latent heat of evaporation (vapor:liquid) at 0C [J/kg]
nam.Li=2.83e6;                                                             % latent heat of sublimation (vapor:solid) at 0C [J/kg]
nam.T0=300;                                                                % ncread(fname,'T00'); % base state temperature [K]
nam.P0=1.e5;                                                               % ncread(fname,'P00'); % base state pressure [Pa]
nam.dx=double(ncreadatt(strcat(folder,files_all(1).name),'/','DX'));       % [m]
nam.dy=double(ncreadatt(strcat(folder,files_all(1).name),'/','DY'));       % [m]
nam.txt = txt;                                                             % A label for the data
nam.dt = 0.5;                                                              % Output frequency [hours]

fprintf('Filename: %s\n', files_all(1).name);
file = files_all(3).name;
fname=strcat(folder,file);

[hfx_for, lh_for, tsk, hfx_tend, lh_tend, grxflx, swdown, glw, swnorm, ...
    hfx, qfx, lh] = loadNetCDF(fname, 'HFX_FORCE', 'LH_FORCE', 'TSK_FORCE', ...
    'HFX_FORCE_TEND', 'LH_FORCE_TEND', 'GRDFLX', 'SWDOWN', 'GLW', ...
    'SWNORM', 'HFX', 'QFX', 'LH');

col_names = {'Parameter', 'Description', 'Value', 'Units'};
% Create a formatted string for the column names
col_str = sprintf('%12s', col_names{:});

names = ['HFX_FORCE ', 'LH_FORCE ', 'TSK_FORCE ', ...
    'HFX_FORCE_TEND ', 'LH_FORCE_TEND ', 'GRDFLX ', 'SWDOWN ', 'GLW ', ...
    'SWNORM ', 'HFX ', 'QFX ', 'LH ']

data = [hfx_for.Data, lh_for.Data, tsk.Data, hfx_tend.Data, lh_tend.Data, ...
    mean(grxflx.Data,'all'), mean(swdown.Data, 'all'), mean(glw.Data, 'all'),...
    mean(swnorm.Data, 'all'), mean(hfx.Data, 'all'), mean(qfx.Data, 'all'),...
    mean(lh.Data, 'all')]