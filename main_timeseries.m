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
nam.dt = 0.5;                   % Output frequency in hours

% Define the folders for each case
cold_pools = ["./data/small_domain/CP_OUT/", "./data/small_domain/NOCP_OUT/"];
%cold_pools = ["/scratch/05999/mkurowsk/GATE_CP_CONSTFLX/", "/scratch/05999/mkurowsk/GATE_NOCP_CONSTFLX/"]

% Get a list of all files in each folder
files_cp = dir(strcat(cold_pools{1}, 'wrfout*'));
files_nocp = dir(strcat(cold_pools{2}, 'wrfout*'));

files_all = files_cp; % choose the folder with shortest output
num_files = length(files_cp);

t_length = 30 %num_files;               % Length of the timeseries (usually num files but can be shorter if the heart desires)

% Preallocate arrays for storing the computed values; index 1 is for CP or
% NOCP
TKE_out = zeros(2, t_length);
LWP_out = zeros(2, t_length);
RAINNC_out = zeros(2, t_length);
PRECR_out = zeros(2, t_length);
PRECG_out = zeros(2, t_length);
CLDFR_out = zeros(2, t_length);
HFX_out = zeros(2, t_length);
QT_out = zeros(2, t_length);
LH_out = zeros(2, t_length);

TPC_out = zeros(2, 256); % two point correlation function

time_hours = zeros(1, t_length);

% Import and calculate vertical structure variables for one of the files
[Z, p, H] = vert_struct(strcat(cold_pools(1),files_cp(1).name), nam);

% Loop over the files and cases
for i = 1:t_length
    
    for j = 1:2

        cp = cold_pools(j);
        fprintf('Filename: %s\n', files_all(i).name);

        file = files_all(i).name;
        fname=strcat(cp,file);
        
        % Read in the things we want 
        try
        rainnc = ncread(fname,'RAINNC');
        tke = ncread(fname,'TKE');
        qv=ncread(fname,'QVAPOR');                                         % Water vapor mixing ratio [kg/kg]
        qc=ncread(fname,'QCLOUD');                                         % Cloud water mixing ratio [kg/kg]
        qr=ncread(fname,'QRAIN');                                          % Rain water mixing ratio [kg/kg]
        qi=ncread(fname,'QICE');                                           % Ice mixing ratio [kg/kg]
        qs=ncread(fname,'QSNOW');                                          % Snow mixing ratio
        qg=ncread(fname,'QGRAUP');                                         % Graupel mixing ratio

        s=size(tke); n=s(1); m=s(2); l=s(3); nm=n*m;
        TKE = mean(reshape(tke,nm,l));                                     % TKE
        LWP = mean_LWP(fname, nam);                                        % LWP


        TKE_out(j, i) = trapz(Z,TKE);
        RAINNC_out(j, i) = mean(rainnc, 'all');
        LWP_out(j, i) = mean(LWP, 'all');
        
        qt=qv+qc+qi;                                                       % total water mixing ratio [kg/kg] (no precipitating elements)

        QT=mean(reshape(qt,nm,l));                                         % Horizontally average moisture variance
        QT_out(j, i) = trapz(Z,QT);                                    % Vertically integrate the moisture vairance (compare to Schemann & Seifert, 2017)
% 
        catch 
% 
        TKE_out(j, i) = -1;
        RAINNC_out(j, i) = -1;
        LWP_out(j, i) = -1;                                                                                          
        QT_out(j, i) = -1; 
        
        end

    end

    time_hours(i) = 0 + (i-1)*nam.dt;


end



%%
figure();

subplot(1,3,1);
plot(time_hours, TKE_out(1, :), 'Linewidth', 1.5); hold on;
plot(time_hours, TKE_out(2, :), 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel('VTKE [m^2 s^{-2}]','LineWidth',1.5,'FontSize',15);
legend('CP', 'NOCP','Location', 'northwest');
title('Vert. integrated TKE', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');

subplot(1,3,2);
plot(time_hours, LWP_out(1,:), 'Linewidth', 1.5); hold on;
plot(time_hours, LWP_out(2,:), 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel('LWP [kg m^{-2}]','LineWidth',1.5,'FontSize',15);
legend('CP', 'NOCP', 'Location', 'northwest');
ylim([45, 60])
title('Vert. integrated LWP', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');

subplot(1,3,3);
plot(time_hours, RAINNC_out(1, :), 'Linewidth', 1.5); hold on;
plot(time_hours, RAINNC_out(2, :), 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
% ylabel(append('RAINNC', ' [', rainnc.Units, ']'),'LineWidth',1.5,'FontSize',15);
ylabel(append('RAINNC', ' [mm]'),'LineWidth',1.5,'FontSize',15);
title('Accumulative precipitation', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
legend('CP', 'NOCP', 'Location', 'northwest')

%%
figure()
plot(time_hours, QT_out(1, :)*1000, 'Linewidth', 1.5); hold on;
plot(time_hours, QT_out(2, :)*1000, 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel('[(g/kg)^2]','LineWidth',1.5,'FontSize',15);
legend('CP', 'NOCP','Location', 'northwest');
title('Vert. avg. Q_t', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');

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

% subplot(2,3,4);
% plot(time_hours, HFX_out(1,:), 'Linewidth', 1.5); hold on;
% plot(time_hours, HFX_out(2,:), 'Linewidth', 1.5);
% xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
% ylabel(append('HFX', ' [Wm^{-2}]'),'LineWidth',1.5,'FontSize',15);
% title('Upward Heat Flux at Surface', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
% legend('CP', 'NOCP', 'Location', 'northwest')
% 
% subplot(2,3,5);
% plot(time_hours, LH_out(1,:), 'Linewidth', 1.5); hold on;
% plot(time_hours, LH_out(2,:), 'Linewidth', 1.5);
% xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
% ylabel(append('LH', ' [Wm^{-2}]'),'LineWidth',1.5,'FontSize',15);
% title('Latent Heat Flux at Surface', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
% legend('CP', 'NOCP', 'Location', 'northwest')
% 
% subplot(2,3,6);
% plot(time_hours, QFX_out(1,:), 'Linewidth', 1.5); hold on;
% plot(time_hours, QFX_out(2,:), 'Linewidth', 1.5);
% xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
% ylabel(append('QFX', ' [kg m^{-2} s^{-1}]'),'LineWidth',1.5,'FontSize',15);
% title('Upward Moisture Flux at Surface', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
% legend('CP', 'NOCP', 'Location', 'northwest')