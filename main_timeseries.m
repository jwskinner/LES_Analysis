%% 
% Script for computing timeseries properties of the convection data 
% Updated by J. W. Skinner (2022-12-14)
%
% convention: 
%   lowercase variables - 3-D or 2-D variables
%   uppercase variables - 1-D profiles
%
% INPUT: fname - filename pointing to WRF-LES output
% OUTPUT: (qtvar,thilvar) - structures with variance terms

clear variables 

cold_pools = ["CP", "NOCP"]; 

output = './plots/timeseries/';

%folder = strcat('./', case(i), '/');
folders = ["./new/TEST/CP_OUT/", "./new/TEST/NOCP_OUT/"];

% Get a list of all files in the current directory                         
files_all = dir(strcat(folders(2), 'wrf*'));

nam.R=287.04;                                                              % [J/kg/K]
nam.cp=1004.67;                                                            % [J/kg/K]
nam.g=9.81;                                                                % [m/s2]
nam.Ll=2.50e6;                                                             % latent heat of evaporation (vapor:liquid) at 0C [J/kg]
nam.Li=2.83e6;                                                             % latent heat of sublimation (vapor:solid) at 0C [J/kg]
nam.T0=300;                                                                % ncread(fname,'T00'); % base state temperature [K]
nam.P0=1.e5;                                                               % ncread(fname,'P00'); % base state pressure [Pa]
nam.dx=double(ncreadatt(strcat(folders(1),files_all(1).name),'/','DX'));    % [m]
nam.dy=double(ncreadatt(strcat(folders(1),files_all(1).name),'/','DY'));    % [m]
nam.dt = 0.5;                                                              % Output frequency in Hours

% Get a list of all files in the current directory                         
files_all = dir(strcat(folders(2), 'wrf*')); 

num_files = length(files_all);

% Preallocate all the empty arrays for filling 
TKE_CP = zeros(num_files);
TKE_NOCP = zeros(num_files);
LWP_CP = zeros(num_files);
LWP_NOCP = zeros(num_files);
RAINNC_CP = zeros(num_files);
RAINNC_NOCP = zeros(num_files); 
PRECR_CP = zeros(num_files);
PRECR_NOCP = zeros(num_files);

PRECG_CP = zeros(num_files);
PRECG_NOCP = zeros(num_files);
time_hours = zeros(num_files);


[Z, p, H] = vert_struct(strcat(folders(1),files_all(1).name), nam);        % Import vertical structure variables

% Loop over the files and cases 
for i = 1:num_files
    
    for j = 1:2
        cp = cold_pools(j);
        %     try
        fprintf('Filename: %s\n', files_all(i).name);
    
        file = files_all(i).name; 
        fname=strcat(folders(j),file); 

        [rainnc, tke, precr, precg] = loadNetCDF(fname, 'RAINNC', 'TKE', 'PRECR', 'PRECG');
        
        s=size(tke.Data); n=s(1); m=s(2); l=s(3); nm=n*m;
        TKE = mean(reshape(tke.Data,nm,l));                                % TKE
        LWP = mean_LWP(fname, nam);                                        % LWP 
        
        if j == 1
            TKE_CP(i) = trapz(Z,TKE);
            RAINNC_CP(i) = mean(rainnc.Data, 'all');
            PRECR_CP(i) = mean(precr.Data, 'all');
            PRECG_CP(i) = mean(precg.Data, 'all');
            LWP_CP(i) = trapz(Z,LWP);
        else
            TKE_NOCP(i) = trapz(Z,TKE);
            RAINNC_NOCP(i) = mean(rainnc.Data, 'all');
            PRECR_NOCP(i) = mean(precr.Data, 'all');
            PRECG_NOCP(i) = mean(precg.Data, 'all');
            LWP_NOCP(i) = trapz(Z,LWP);
        end        

    end
        time_hours(i) = 0 + (i-1)*nam.dt; 
%    catch
%        fprintf('File %s failed\n',file)
%        outlwp(i) = 0;
%        outlist(i) = 0; 
%    end

end


%%
figure();

subplot(2,3,1); 
plot(time_hours, TKE_CP, 'Linewidth', 1.5); hold on; 
plot(time_hours, TKE_NOCP, 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel('VTKE [m^2 s^{-2}]','LineWidth',1.5,'FontSize',15);
legend('CP', 'NOCP','Location', 'northwest');
title('Vert. integrated TKE', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal'); 


subplot(2,3,2); 
plot(time_hours, LWP_CP, 'Linewidth', 1.5); hold on; 
plot(time_hours, LWP_NOCP, 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel('LWP [kg m^{-2}]','LineWidth',1.5,'FontSize',15);
legend('CP', 'NOCP', 'Location', 'northwest'); 
ylim([50, 100])
title('Vert. integrated LWP', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');

subplot(2,3,3); 
plot(time_hours, RAINNC_CP, 'Linewidth', 1.5); hold on; 
plot(time_hours, RAINNC_NOCP, 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel(append('RAINNC', ' [', rainnc.Units, ']'),'LineWidth',1.5,'FontSize',15);
title('Accumulative precipitation', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
legend('CP', 'NOCP', 'Location', 'northwest') 

subplot(2,3,4); 
plot(time_hours, PRECR_CP, 'Linewidth', 1.5); hold on; 
plot(time_hours, PRECR_NOCP, 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel(append('PRECR', ' [', precr.Units, ']'),'LineWidth',1.5,'FontSize',15);
title('Rain precip rate', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
legend('CP', 'NOCP', 'Location', 'northwest') 

subplot(2,3,5); 
plot(time_hours, PRECG_CP, 'Linewidth', 1.5); hold on; 
plot(time_hours, PRECG_NOCP, 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel(append('PRECG', ' [', precr.Units, ']'),'LineWidth',1.5,'FontSize',15);
title('Graupel precip rate', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
legend('CP', 'NOCP', 'Location', 'northwest') 