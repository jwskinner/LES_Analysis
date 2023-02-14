% Script for computing timeseries properties of the convection data
% Updated by J. W. Skinner (2022-12-14)
%
% INPUT: fname - filename pointing to WRF-LES output
% OUTPUT: (qtvar,thilvar) - structures with variance terms

% Define constants
nam.R = 287.04;                 % [J/kg/K]
nam.cp = 1004.67;               % [J/kg/K]
nam.g = 9.81;                   % [m/s^2]
nam.Ll = 2.50e6;                % latent heat of evaporation (vapor:liquid) at 0C [J/kg]
nam.Li = 2.83e6;                % latent heat of sublimation (vapor:solid) at 0C [J/kg]
nam.T0 = 300;                   % base state temperature [K]
nam.P0 = 1.e5;                  % base state pressure [Pa]
nam.dt = 4.0;                   % Output frequency in hours

% Define the folders for each case
cold_pools = ["./data/large_domain/CP_OUT/", "./data/large_domain/NOCP_OUT/"];

% Get a list of all files in each folder
files_cp = dir(strcat(cold_pools{1}, 'wrfout*'));
files_nocp = dir(strcat(cold_pools{2}, 'wrfout*'));

files_all = files_cp; % choose the folder with shortest output

% Preallocate arrays for storing the computed values; index 1 is for CP or
% NOCP
num_files = 6 %length(files_nocp);
TKE_out = zeros(2, num_files);
LWP_out = zeros(2, num_files);
RAINNC_out = zeros(2, num_files);
PRECR_out = zeros(2, num_files);
PRECG_out = zeros(2, num_files);
CLDFR_out = zeros(2, num_files);
HFX_out = zeros(2, num_files);
QFX_out = zeros(2, num_files);
LH_out = zeros(2, num_files);

TPC_out = zeros(2, 256); % two point correlation function

time_hours = zeros(1, num_files);

% Import and calculate vertical structure variables for one of the files
[Z, p, H] = vert_struct(strcat(cold_pools(1),files_cp(1).name), nam);

% Loop over the files and cases
parfor i = 1:num_files

    for j = 1:2

        cp = cold_pools(j);
        try 
        fprintf('Filename: %s\n', files_all(i).name);

        file = files_all(i).name;
        fname=strcat(cp,file);

        [rainnc, tke] = loadNetCDF(fname, 'RAINNC', 'TKE');

        [hfx, qfx, lh] = loadNetCDF(fname, 'HFX', 'QFX', 'LH');

        s=size(tke.Data); n=s(1); m=s(2); l=s(3); nm=n*m;
        TKE = mean(reshape(tke.Data,nm,l));                                % TKE
        LWP = mean_LWP(fname, nam);                                        % LWP
%         CLDFR = mean(reshape(cldfr.Data,nm,l));


        TKE_out(j, i) = trapz(Z,TKE);
        RAINNC_out(j, i) = mean(rainnc.Data, 'all');
%         PRECR_out(j, i) = mean(precr.Data, 'all');
%         PRECG_out(j, i) = mean(precg.Data, 'all');
        LWP_out(j, i) = trapz(Z,LWP);
%         CLDFR_out(j, i) = trapz(Z,CLDFR);

        HFX_out(j, i) = mean(hfx.Data, 'all');
        QFX_out(j, i) = mean(qfx.Data, 'all');
        LH_out(j, i) = mean(lh.Data, 'all');

        %% 
        % Compute two point correlation functions
%         [u] = loadNetCDF(fname, 'U');
%         TPC_out(j, :) = two_point_cor(u, nam);
% 
%         figure()
%         plot(TPC_out(1, :), 'Linewidth', 1.5); hold on;
%         plot(TPC_out(2, :), 'Linewidth', 1.5); hold on;

        catch
        TKE_out(j, i) = -1;
        RAINNC_out(j, i) = -1;
        LWP_out(j, i) = -1;
        HFX_out(j, i) = -1;
        QFX_out(j, i) = -1;
        LH_out(j, i) = -1;
        fprintf('File %s failed\n',file)
        end


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

% subplot(2,4,4);
% plot(time_hours, PRECR_out(1, :), 'Linewidth', 1.5); hold on;
% plot(time_hours, PRECR_out(2, :), 'Linewidth', 1.5);
% xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
% ylabel(append('PRECR', ' [', precr.Units, ']'),'LineWidth',1.5,'FontSize',15);
% title('Rain precip rate', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
% legend('CP', 'NOCP', 'Location', 'northwest')
% 
% subplot(2,4,5);
% plot(time_hours, PRECG_out(1,:), 'Linewidth', 1.5); hold on;
% plot(time_hours, PRECG_out(2,:), 'Linewidth', 1.5);
% xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
% ylabel(append('PRECG', ' [', precg.Units, ']'),'LineWidth',1.5,'FontSize',15);
% title('Graupel precip rate', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
% legend('CP', 'NOCP', 'Location', 'northwest')

subplot(2,3,4);
plot(time_hours, HFX_out(1,:), 'Linewidth', 1.5); hold on;
plot(time_hours, HFX_out(2,:), 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel(append('HFX', ' [', hfx.Units, ']'),'LineWidth',1.5,'FontSize',15);
title(hfx.Desc, 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
legend('CP', 'NOCP', 'Location', 'northwest')

subplot(2,3,5);
plot(time_hours, LH_out(1,:), 'Linewidth', 1.5); hold on;
plot(time_hours, LH_out(2,:), 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel(append('LH', ' [', lh.Units, ']'),'LineWidth',1.5,'FontSize',15);
title(lh.Desc, 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
legend('CP', 'NOCP', 'Location', 'northwest')

subplot(2,3,6);
plot(time_hours, QFX_out(1,:), 'Linewidth', 1.5); hold on;
plot(time_hours, QFX_out(2,:), 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel(append('QFX', ' [', qfx.Units, ']'),'LineWidth',1.5,'FontSize',15);
title(qfx.Desc, 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
legend('CP', 'NOCP', 'Location', 'northwest')