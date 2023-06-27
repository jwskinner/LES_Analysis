%% J.W.Skinner 2023 
% This script reads in data from specified location and creates time averaged 
% histograms for distribution of each of the specified variables U, V, W, Theta and
% TKE. 

clear variables 

data_loc = '/data1/jwskinner/';
folder = 'GATE_NOEVP1.3km_CONSTFLX_100km/' %'GATE_NOCP_CONSTFLX/';

% index of file to start and finish the averaging over; set to same for
% instantanious histograms 
i_start = 90; 
i_end = 97; % 48 hours

% Takes the first folder for loading in params.
f_in = strcat(data_loc, folder); 

% Returns a list of all files in the folder
files_all = dir(strcat(f_in, 'wrfout*'));

% Load relevant physical constants into nam structure 
[nam] = get_constants(strcat(f_in, files_all(1).name));

% Load the vertical structure from a sample file
[Z, P, H] = vert_struct(strcat(f_in, files_all(1).name), nam);

% Pre-allocate arrays for the time averaging
u_t = zeros((i_end - i_start), nam.nx-1, nam.ny, nam.levs); 
v_t = zeros((i_end - i_start), nam.nx-1, nam.ny, nam.levs);
w_t = zeros((i_end - i_start), nam.nx-1, nam.ny, nam.levs);
th_t = zeros((i_end - i_start), nam.nx-1, nam.ny, nam.levs);
tke_t = zeros((i_end - i_start), nam.nx-1, nam.ny, nam.levs);

j = 1;
for i = i_start: i_end

    file = files_all(i).name
    fname=strcat(f_in,file);

    [u, v, w, th, tke, qc, qr] =  ...
        loadNetCDF(fname, 'U', 'V', 'W', 'T', 'TKE', 'QCLOUD', 'QRAIN');

    % Extract data component and convert to grid
    u = u.Data; v = v.Data; w = w.Data; th = th.Data + nam.T0; ...
        tke = tke.Data; qc = qc.Data; qr = qr.Data;  
    
    %C->A grid
    u=0.5*(u(1:end-1,:,:,:)+u(2:end,:,:,:)); 
    v=0.5*(v(:,1:end-1,:,:)+v(:,2:end,:,:)); 
    w=0.5*(w(:,:,1:end-1,:)+w(:,:,2:end,:));
    
    % Build the time arrays 
    u_t(j,:,:,:) = u; v_t(j,:,:,:) = v; w_t(j,:,:,:) = w; ... 
        th_t(j,:,:,:) = th; tke_t(j,:,:,:) = tke;
    
    j = j + 1; 
end

% Average over the time dimension of the arrays 
u_tmean = mean(u_t, 1); 
v_tmean = mean(v_t, 1);
w_tmean = mean(w_t, 1);
th_tmean = mean(th_t, 1);
tke_tmean = mean(tke_t, 1);

%% Make the histogram plots of velocity fields 
% Define variables and data array
variables = {'u', 'v', 'w'};
data = [u_tmean; v_tmean; w_tmean]; 

% Set histogram properties
binEdges = linspace(min(data,[], 'all'), max(data,[],'all'), 1000);         % Define custom bin edges
faceColor = [0, 0, 0.7]; % Custom face color
edgeColor = 'none'; % Black edge color

% Plot histograms of the velocities
figure('Renderer', 'painters', 'Position', [10 10 900 600])                % makes paper format figure

for i = 1:numel(variables)
    subplot(1, 3, i)
    % Create the histogram plot
    histogram(data(i,:, :, :), binEdges, 'FaceColor', faceColor, 'EdgeColor', edgeColor);
    % Add labels and title
    xlabel(['$', variables{i}, '$', ' [m/s]'], 'Interpreter', 'Latex');
    ylabel('Frequency');
%     xlim([-5, 5]);
    ylim([0, 4*10^6]);
    % Adjust appearance
    grid on;
    set(gca, 'LineWidth', 1.25, 'FontName', 'Arial', 'FontSize', 12, ...
        'XGrid', 'on', 'YGrid', 'off', 'GridLineStyle', '--', ...
        'MinorGridLineStyle', 'none');
    set(gca, 'GridAlpha', 0.25, 'XMinorTick', 'off');
end

set(gcf, 'Color', 'w');

%% Make the histogram plots of potential temperature field 
% Define variables and data array
data = [th_tmean]; 

% Set histogram properties
binEdges = linspace(min(data,[], 'all'), max(data,[],'all'), 200);         % Define custom bin edges
faceColor = [0.4, 0.6, 0.8]; % Custom face color
edgeColor = 'k'; % Black edge color

% Plot histograms of the velocities
figure('Renderer', 'painters', 'Position', [10 10 900 600])                % makes paper format figure

% Create the histogram plot
histogram(squeeze(data), binEdges, 'FaceColor', faceColor, 'EdgeColor', edgeColor);
% Add labels and title
xlabel(['$\theta$', ' [K]'], 'Interpreter', 'Latex');
ylabel('Frequency');
xlim([300, 400]);
ylim([0, 5000000]);
% Adjust appearance
grid on;
set(gca, 'LineWidth', 1.25, 'FontName', 'Arial', 'FontSize', 12, ...
    'XGrid', 'off', 'YGrid', 'off', 'GridLineStyle', '--', ...
    'MinorGridLineStyle', 'none');
set(gca, 'GridAlpha', 0.5, 'XMinorTick', 'off');
set(gcf, 'Color', 'w');




