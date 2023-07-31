%% J.W.Skinner 07/27/2023 
% This script reads in data from specified location and time averages 
% the CKE budgets from the Kardunov and Randall paper
%
% Todos: - speed up code; chec tke budgets 

clear variables 

addpath('../cmocean-main/')
addpath('./functs/')

data_loc = '/data1/jwskinner/';
folder = 'GATE_CP_CONSTFLX/'; % CP case 100km
folder = 'GATE_NOEVP1.3km_CONSTFLX_100km/'; %NOCP no aggregation case 100km

% index of file to start and finish the averaging over
i_start = 90; 
i_end = 90 %97; % 48 hours

% Takes the first folder for loading in params.
f_in = strcat(data_loc, folder); 

% Returns a list of all files in the folder
files_all = dir(strcat(f_in, 'wrfout*'));

% Load relevant physical constants into nam structure 
[nam] = get_constants(strcat(f_in, files_all(1).name));

% Load the vertical structure 
[Z, P, H] = vert_struct(strcat(f_in, files_all(1).name), nam);

% Pre-allocate arrays for the time averaging
cke_t = zeros((i_start - i_end), 6, nam.levs);  % dimensions:(time, # budget terms, nlevs)

j = 1;
for i = i_start: i_end

    file = files_all(i).name
    fname=strcat(f_in,file);

    % Sets up the paramters for the plotting
    [cke_var] = comp_cke_budgets(fname, nam, i); 
    
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

% Non dimensionalisation of the budget terms by 

% Plot the budgets 
%% first plot the budgets over time
figure('Renderer', 'painters', 'Position', [10 10 900 600]) % makes paper format figure

lw = 1.5
colors = {'#0072BD', '#A2142F', '#EDB120', '#77AC30', '#800080', '#000000'}; 
style = {'-', '-', '-', '-', '-', '--'};

% DIAGNOSTIC PLOT OF TIME VARYING BUDGETS 
alpha = 1.0; 
rgbC = [];
for i = 1:6 
    rgb = sscanf(colors{i}(2:end), '%2x')/255;
    rgbC = [rgbC, rgb]; 
end

subplot(1, 1, 1)
for i = 1:(i_end-i_start)
    plot(squeeze(cke_t(i, 1, :)),Z/10^3,'LineWidth',lw, 'Color', [rgbC(:,1)', alpha]); hold on; % B term
    plot(squeeze(cke_t(i, 2, :)),Z/10^3,'LineWidth',lw, 'Color', [rgbC(:,2)', alpha]); hold on; % T term
    plot(squeeze(cke_t(i, 3, :)),Z/10^3,'LineWidth',lw, 'Color', [rgbC(:,3)', alpha]); hold on; % S term
    plot(squeeze(cke_t(i, 4, :)),Z/10^3,'LineWidth',lw, 'Color',  [rgbC(:,4)', alpha]); hold on; %P term
    plot(squeeze(cke_t(i, 5, :)),Z/10^3,'LineWidth',lw, 'Color',  [rgbC(:,5)', alpha]); hold on; %D term
    plot(squeeze(cke_t(i, 6, :)),Z/10^3,'--','LineWidth',lw, 'Color',  [rgbC(:,6)', alpha]); hold on; %R term
end
ylabel('$z$ [km]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
legend('Bouyancy','Transport','Shear','Pressure corr.','Dissipation', 'Residual');
title('$\langle \overline{\rm CKE} \rangle$ variance budget', 'FontSize',15, 'Interpreter', 'Latex'); hold off;
xlim([-1e-3 1e-3])
ylim([0 18])



% MAIN PLOT OF TIME AVERAGED BUDGETS 
figure('Renderer', 'painters', 'Position', [10 10 900 600]) % makes paper format figure
subplot(1, 1, 1)
for i = 1:6
    plot(squeeze(cke_tmean(:, i, :)),Z/10^3,style{i},'LineWidth',lw,'Color',colors{i}); hold on; % Dissipation terms
end 
ylabel('$z$ [km]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
legend('Bouyancy','Transport','Shear','Pressure corr.','Dissipation', 'Residual');
title('$t$ average $\langle \overline{\rm CKE} \rangle$ variance budget', 'FontSize',15, 'Interpreter', 'Latex'); hold off;
xlim([-1e-3 1e-3])
ylim([0 18])





