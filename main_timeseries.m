% Script for computing timeseries properties of the convection data
% Updated by J. W. Skinner (2022-12-14)
%
% INPUT: fname - filename pointing to WRF-LES output
% OUTPUT: (qtvar,thilvar) - structures with variance terms

% Define constants
R = 287.04;                 % [J/kg/K]
cp = 1004.67;               % [J/kg/K]
g = 9.81;                   % [m/s^2]
Ll = 2.50e6;                % latent heat of evaporation (vapor:liquid) at 0C [J/kg]
Li = 2.83e6;                % latent heat of sublimation (vapor:solid) at 0C [J/kg]
T0 = 300;                   % base state temperature [K]
P0 = 1.e5;                  % base state pressure [Pa]
dt = 0.5;                   % Output frequency in hours

% Define the folders for each case
cold_pools = {"./data/small_domain/CP_OUT/", "./data/small_domain/NOCP_OUT/"};

% Get a list of all files in each folder
files_cp = dir(strcat(cold_pools{1}, 'wrfout*'));
files_nocp = dir(strcat(cold_pools{2}, 'wrfout*'));

% Preallocate arrays for storing the computed values; index 1 is for CP or
% NOCP
num_files = length(files_nocp);
TKE_out = zeros(2, num_files);
LWP_out = zeros(2, num_files);
RAINNC_out = zeros(2, num_files);
PRECR_out = zeros(2, num_files);
PRECG_out = zeros(2, num_files);
CLDFR_out = zeros(2, num_files);

time_hours = zeros(1, num_files);

% Import and calculate vertical structure variables for one of the files
[Z, p, H] = vert_struct(strcat(folders(1),files_all(1).name), nam);

% Loop over the files and cases
for i = 1:num_files

    for j = 1:2

        cp = cold_pools(j);

        fprintf('Filename: %s\n', files_all(i).name);

        file = files_all(i).name;
        fname=strcat(folders(j),file);

        [rainnc, tke, precr, precg, cldfr] = loadNetCDF(fname, 'RAINNC', 'TKE', 'PRECR', 'PRECG', 'CLDFRA');

        s=size(tke.Data); n=s(1); m=s(2); l=s(3); nm=n*m;
        TKE = mean(reshape(tke.Data,nm,l));                                % TKE
        LWP = mean_LWP(fname, nam);                                        % LWP
        CLDFR = mean(reshape(cldfr.Data,nm,l));


        TKE_out(j, i) = trapz(Z,TKE);
        RAINNC_out(j, i) = mean(rainnc.Data, 'all');
        PRECR_out(j, i) = mean(precr.Data, 'all');
        PRECG_out(j, i) = mean(precg.Data, 'all');
        LWP_out(j, i) = trapz(Z,LWP);
        CLDFR_out(j, i) = trapz(Z,CLDFR);


    end
    time_hours(i) = 0 + (i-1)*nam.dt;


end


%%
figure();

subplot(2,3,1);
plot(time_hours, TKE_out(1, :), 'Linewidth', 1.5); hold on;
plot(time_hours, TKE_out(2, :), 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel('VTKE [m^2 s^{-2}]','LineWidth',1.5,'FontSize',15);
legend('CP', 'NOCP','Location', 'northwest');
title('Vert. integrated TKE', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');


subplot(2,3,2);
plot(time_hours, LWP_out(1,:), 'Linewidth', 1.5); hold on;
plot(time_hours, LWP_out(2,:), 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel('LWP [kg m^{-2}]','LineWidth',1.5,'FontSize',15);
legend('CP', 'NOCP', 'Location', 'northwest');
ylim([50, 100])
title('Vert. integrated LWP', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');

subplot(2,3,3);
plot(time_hours, RAINNC_out(1, :), 'Linewidth', 1.5); hold on;
plot(time_hours, RAINNC_out(2, :), 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel(append('RAINNC', ' [', rainnc.Units, ']'),'LineWidth',1.5,'FontSize',15);
title('Accumulative precipitation', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
legend('CP', 'NOCP', 'Location', 'northwest')

subplot(2,3,4);
plot(time_hours, PRECR_out(1, :), 'Linewidth', 1.5); hold on;
plot(time_hours, PRECR_out(2, :), 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel(append('PRECR', ' [', precr.Units, ']'),'LineWidth',1.5,'FontSize',15);
title('Rain precip rate', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
legend('CP', 'NOCP', 'Location', 'northwest')

subplot(2,3,5);
plot(time_hours, PRECG_out(1,:), 'Linewidth', 1.5); hold on;
plot(time_hours, PRECG_out(2,:), 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel(append('PRECG', ' [', precg.Units, ']'),'LineWidth',1.5,'FontSize',15);
title('Graupel precip rate', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
legend('CP', 'NOCP', 'Location', 'northwest')

subplot(2,3,6);
plot(time_hours, CLDFR_out(1,:), 'Linewidth', 1.5); hold on;
plot(time_hours, CLDFR_out(2,:), 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel(append('CLDFR', ' [', cldfr.Units, ']'),'LineWidth',1.5,'FontSize',15);
title('Cloud Fraction', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
legend('CP', 'NOCP', 'Location', 'northwest')