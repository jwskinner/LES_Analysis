%% J.W.Skinner 2023 
% This script reads in data from specified location and time averages 
% the moisture variance and liquid water potential temperature variance 
% budgets then plots them. 
clear variables 

addpath('../cmocean-main/')
addpath('./functs/')

data_loc = '/data1/jwskinner/';
%folder = 'GATE_CP_CONSTFLX/'; 
folder = 'GATE_NOEVP1.3km_CONSTFLX_100km/';

% index of file to start and finish the averaging over
i_start = 89; 
i_end = 97; % 48 hours

% Takes the first folder for loading in params.
f_in = strcat(data_loc, folder); 

% Returns a list of all files in the folder
files_all = dir(strcat(f_in, 'wrfout*'));

% Load relevant physical constants into nam structure 
[nam] = get_constants(strcat(f_in, files_all(1).name));

% Load the vertical structure 
[Z, P, H] = vert_struct(strcat(f_in, files_all(1).name), nam);

% Pre-allocate arrays for the time averaging
qtvar_t = zeros((i_start - i_end), 5, nam.levs);  % dimensions:(time, # budget terms, nlevs)
thilvar_t = zeros((i_start - i_end), 5, nam.levs);

j = 1;
for i = i_start: i_end

    file = files_all(i).name
    fname=strcat(f_in,file);

    % Sets up the paramters for the plotting
    [qtvar, thilvar] = comp_moist_budgets(fname, nam, i); 
    
    % Assign to time varying arrays
    qtvar_t(j, 1, :) = qtvar.prod; 
    qtvar_t(j, 2, :) = qtvar.trns;
    qtvar_t(j, 3, :) = qtvar.srcs;
    qtvar_t(j, 4, :) = qtvar.diss;
    qtvar_t(j, 5, :) = qtvar.prod + qtvar.trns +qtvar.srcs +qtvar.diss; %residual 
    
    thilvar_t(j, 1, :) = thilvar.prod; 
    thilvar_t(j, 2, :) = thilvar.trns;
    thilvar_t(j, 3, :) = thilvar.srcs;
    thilvar_t(j, 4, :) = thilvar.diss;
    thilvar_t(j, 5, :) = thilvar.prod + thilvar.trns + thilvar.srcs + thilvar.diss;
    
    j = j + 1; 
end

% Average over the time dimension of the arrays 
qtvar_tmean = mean(qtvar_t, 1); 
thilvar_tmean = mean(thilvar_t, 1); 

%% Standard deviation over the time dimension 
% qtvar_std = std(qtvar_t, 0, 1);
% thilvar_std = std(thilvar, 0, 1);

% Plot the budgets 
%% first plot the budgets over time
figure('Renderer', 'painters', 'Position', [10 10 900 600]) % makes paper format figure

lw = 1.5
colors = {'#0072BD', '#A2142F', '#EDB120', '#77AC30', '#000000'}; 
style = {'-', '-', '-', '-', '--'};

% DIAGNOSTIC PLOT OF TIME VARYING BUDGETS 
alpha = 1.0; 
rgbC = []
for i = 1:5 
    rgb = sscanf(colors{i}(2:end), '%2x')/255;
    rgbC = [rgbC, rgb]; 
end

subplot(1, 2, 1)
for i = 1:(i_end-i_start)
    plot(squeeze(qtvar_t(i, 1, :)),Z/10^3,'LineWidth',lw, 'Color', [rgbC(:,1)', alpha]); hold on; % Production terms
    plot(squeeze(qtvar_t(i, 2, :)),Z/10^3,'LineWidth',lw, 'Color', [rgbC(:,2)', alpha]); hold on; % Transport terms
    plot(squeeze(qtvar_t(i, 3, :)),Z/10^3,'LineWidth',lw, 'Color', [rgbC(:,3)', alpha]); hold on; % Source/Sink terms
    plot(squeeze(qtvar_t(i, 4, :)),Z/10^3,'LineWidth',lw, 'Color',  [rgbC(:,4)', alpha]); hold on; % Dissipation terms
end
ylabel('$z$ [km]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
legend('PROD','TURB','SRC','DISS','RES');
title('$\langle \overline{q_t} \rangle$ variance budget', 'FontSize',15, 'Interpreter', 'Latex'); hold off;
xlim([-1e-3 1e-3])
ylim([0 18])

subplot(1, 2, 2)
for i = 1:(i_end-i_start)
    plot(squeeze(thilvar_t(i, 1, :)),Z/10^3,'LineWidth',lw, 'Color', [rgbC(:,1)', alpha]); hold on; % Production terms
    plot(squeeze(thilvar_t(i, 2, :)),Z/10^3,'LineWidth',lw, 'Color', [rgbC(:,2)', alpha]); hold on; % Transport terms
    plot(squeeze(thilvar_t(i, 3, :)),Z/10^3,'LineWidth',lw, 'Color', [rgbC(:,3)', alpha]); hold on; % Source/Sink terms
    plot(squeeze(thilvar_t(i, 4, :)),Z/10^3,'LineWidth',lw, 'Color', [rgbC(:,4)', alpha]); hold on; % Dissipation terms
end
ylabel('$z$ [km]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
legend('PROD','TURB','SRC','DISS','RES');
title('$\langle \overline{\theta_l} \rangle$ variance budget', 'FontSize',15, 'Interpreter', 'Latex'); hold off;
xlim([-1e-3 1e-3])
ylim([0 18])


% MAIN PLOT OF TIME AVERAGED BUDGETS 
figure('Renderer', 'painters', 'Position', [10 10 900 600]) % makes paper format figure
subplot(1, 2, 1)
for i = 1:5
    plot(squeeze(qtvar_tmean(:, i, :)),Z/10^3,style{i},'LineWidth',lw,'Color',colors{i}); hold on; % Dissipation terms
end 
ylabel('$z$ [km]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
% legend('Production','Turbulent transport','Microphysical sources/sinks','Dissipation','Residual');
legend('PROD','TURB','SRC','DISS','RES');
title('$t$ average $\langle \overline{q_t} \rangle$ variance budget', 'FontSize',15, 'Interpreter', 'Latex'); hold off;
xlim([-1e-3 1e-3])
ylim([0 18])

subplot(1, 2, 2)
for i = 1:5
    plot(squeeze(thilvar_tmean(:, i, :)),Z/10^3,style{i},'LineWidth',lw,'Color', colors{i}); hold on; % Dissipation terms
end 

ylabel('$z$ [km]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
% legend('Production','Turbulent transport','Microphysical sources/sinks','Dissipation','Residual');
legend('PROD','TURB','SRC','DISS','RES');
title('$t$ average $\langle \overline{\theta_l} \rangle$ variance budget', 'FontSize',15, 'Interpreter', 'Latex'); hold off;
xlim([-1e-3 1e-3])
ylim([0 18])



