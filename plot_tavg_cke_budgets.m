%% J.W.Skinner 07/27/2023 
% This script reads in data from specified location and time averages 
% the CKE budgets from the Kardunov and Randall paper
% 
% The script has three main parts: 1) Import all the data; 2) compute
% budgets and time average them; 3) plot the budgets. The functions for
% these are in the ./functs directory. 
%
% Todos: - speed up code; chec tke budgets 

clear variables 

addpath('../cmocean-main/')
addpath('./functs/')

%data_loc = '/data1/jwskinner/';                                           % Data location on atmosphere cluster
data_loc = './data/small_domain/';                                         % Data location on jacks laptop

% Folders for atmsophere cluster
%folder = 'GATE_CP_CONSTFLX/';                                              % CP case 100km
folder = 'GATE_NOEVP1.3km_CONSTFLX_100km/';                                % NOCP no aggregation case 100km

% Folders for Jacks laptop
folder = 'CP_OUT/';                                                        % CP case 100km
% folder = 'NOCP_OUT/';                                                    % NOCP no aggregation case 100km

% index of file to start and finish the averaging over
i_start = 10; 
i_end = 11; %90; % 48 hours

% Number of budget terms
budgets = {'Bouyancy','Transport','Shear','Pressure corr.','Dissipation', 'Residual'};
title_str = '$0.5 \rho \langle \overline{\rm u^{\prime \, 2} + v^{\prime \, 2} } \rangle$'; 
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
cke_t = zeros((i_start - i_end), n_budgets, nam.levs);                     % dimensions:(time, # budget terms, nlevs)

j = 1;
for i = i_start: i_end

    file = files_all(i).name
    fname=strcat(f_in,file);

    % Sets up the paramters for the plotting
    [cke_var, scaling] = comp_cke_budgets(fname, nam, i); 
    
    % Assign to time varying arrays
    cke_t(j, 1, :) = cke_var.Bterm; 
    cke_t(j, 2, :) = cke_var.Tterm;
    cke_t(j, 3, :) = cke_var.Sterm;
    cke_t(j, 4, :) = cke_var.Pterm;
    cke_t(j, 5, :) = cke_var.Dterm;
    cke_t(j, 6, :) = cke_var.Res;    %residual 
    
    j = j + 1; 
end



%% Average over the time dimension of the arrays 
cke_tmean = mean(cke_t, 1);  

%% Plot the budget terms 
plot_t_av_budgets(cke_t, scaling, Z, i_start, i_end, n_budgets, budgets, title_str)

