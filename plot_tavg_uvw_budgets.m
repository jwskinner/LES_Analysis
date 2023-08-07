%% J.W.Skinner 07/27/2023
% This script reads in data from specified location and time averages
% the momentum flux budgets from equation 12 of K&R. It loads data from
% either the atmosphere cluster or a local data folder,
% depending on the chosen data location. The script performs time averaging
% over specified time indices and calculates various budget terms related
% to momentum flux. It then creates various plots to visualize the results,
% including time series plots of the time-varying budgets, a main plot of
% time-averaged budgets, and subplots of individual time-averaged budget
% terms. 

clear variables 

addpath('../cmocean-main/')
addpath('./functs/')

%data_loc = '/data1/jwskinner/';                                           % Data location on atmosphere cluster
data_loc = './data/small_domain/';                                         % Data location on jacks laptop

% Folders for atmsophere cluster
folder = 'GATE_CP_CONSTFLX/';                                              % CP case 100km
%folder = 'GATE_NOEVP1.3km_CONSTFLX_100km/';                                % NOCP no aggregation case 100km

% Folders for Jacks laptop
folder = 'CP_OUT/';                                                        % CP case 100km
% folder = 'NOCP_OUT/';                                                    % NOCP no aggregation case 100km

% index of file to start and finish the averaging over
i_start = 10; 
i_end = 11 %97; % 48 hours

% Number of budget terms
budgets = {'Transport', 'Shear', 'Pressure redis.', 'Dissipation', 'Residual'};
title = '$0.5 \rho \langle \overline{\rm u^{\prime \, 2} + v^{\prime \, 2} } \rangle$'; 
n_budgets = sum(cellfun(@ischar, budgets)); 

% Takes the first folder for loading in params.
f_in = strcat(data_loc, folder); 

% Returns a list of all files in the folder
files_all = dir(strcat(f_in, 'wrfout*'));

% Load relevant physical constants into nam structure 
[nam] = get_constants(strcat(f_in, files_all(1).name));

% Load the vertical structure 
[Z, P, H] = vert_struct(strcat(f_in, files_all(1).name), nam);

% Pre-allocate arrays for the time averaging
mom_t = zeros((i_start - i_end), n_budgets, nam.levs);                     % dimensions:(time, # budget terms, nlevs)

j = 1;
for i = i_start: i_end

    file = files_all(i).name
    fname=strcat(f_in,file);

    % Call function to compute the budgets 
    [uv_var] = comp_uvw_budgets(fname, nam, i); 
    
    % Assign to time varying arrays
    uv_t(j, 1, :) = uv_var.Tterm;
    uv_t(j, 2, :) = uv_var.Sterm;
    uv_t(j, 3, :) = uv_var.Rterm;
    uv_t(j, 4, :) = uv_var.Dterm;
    uv_t(j, 5, :) = uv_var.Res;     
    
    j = j + 1; 
end

%% PLOT ALL THE BUDGET TERMS
plot_t_av_budgets(uv_t, Z, i_start, i_end, n_budgets, budgets, title)

  
