% LES Analysis for studying LES data in MATLAB
% This main script loops over the data and performs analysis frame by frame
% there are multiple options (e.g., Fields, Budgets, LWP, 1D profiles,
% Vertical Structure, etc.) that can be chosen with the plot_out parameter.

clear variables

addpath(['../cmocean-main/'])
addpath('./functs/')

% scratch = "/scratch/05999/mkurowsk/"; % For Tacc
scratch = "./data/"; % For jacks laptop

 folders = ["./small_domain/CP_OUT/", "./small_domain/NOCP_OUT/"]
%folders = ["./GATE_NOCP_int_domain/", "./GATE_CP_int_domain/"];
%folders = ["ocean_cp/", "ocean_nocp/"];
% folders = ["GATE_NOEVP1km_CONSTFLX_50km/", "GATE_NOEVP1.3km_CONSTFLX_100km/", "GATE_NOEVP1km_CONSTFLX_100km/"];
% folders = ["GATE_NOEVP1km_CONSTFLX_50km/"];

% Takes the first folder for loading in params.
folder = strcat(scratch, folders(1));  %'./data/large_domain/CP_OUT/' %'/scratch/05999/mkurowsk/ocean_nocp/' on TACC;
output = './plots/';

% Returns a list of all files in the folder
files_all = dir(strcat(folder, 'wrfout*'));
sample_file = strcat(folder,files_all(1).name);                            % File for loading in all the numerical parameters of the simulation

% Setup physical and numerical constants
nam.R=287.04;                                                              % [J/kg/K]
nam.cp=1004.67;                                                            % [J/kg/K]
nam.g=9.81;                                                                % [m/s2]
nam.Ll=2.50e6;                                                             % latent heat of evaporation (vapor:liquid) at 0C [J/kg]
nam.Li=2.83e6;                                                             % latent heat of sublimation (vapor:solid) at 0C [J/kg]
nam.T0=300;                                                                % ncread(fname,'T00'); % base state temperature [K]
nam.P0=1.e5;                                                               % ncread(fname,'P00'); % base state pressure [Pa]
nam.dx=double(ncreadatt(sample_file,'/','DX'));                            % [m]
nam.dy=double(ncreadatt(sample_file,'/','DY'));                            % [m]
nam.dt = 0.5;                                                              % Output frequency [hours]
nam.levs = size(ncread(sample_file,'U'), 3);                               % Number of vertical levels in the simulation
nam.nx = size(ncread(sample_file,'U'), 1);                                 % Number of x grid points in simulation
nam.ny = size(ncread(sample_file,'U'), 2);                                 % Number of y grid points in simulation
nam.txt = 'NOCP';

% Plot diagnostic
plot_out = 2;  % 0: Fields, 1: Budgets, 2:LWP, 3:1D profiles, 4: Vertical Structure, 5: KE Spectra

% For creating movies frame by frame
for i = 1:103 %length(files_all)

    file = files_all(i).name;
    fname=strcat(folder,file);

    % Sets up the paramters for the plotting
    fprintf('Filename: %s\n', files_all(i).name);

    time = (i-1)*nam.dt;

    %% -- Compute and plot files frame by frame into a movie --
    if plot_out == 0

        params.time = time;
        params.name = file;
        params.absum = 0;

        % Load the data
        [variable, variable_u, variable_v, variable_qc] = loadNetCDF(fname, 'W', 'U', 'V', 'QSNOW');
        params.save_folder = append(output, variable.Name, '/');

        % Paramters specific to this plot
        params.cmax = 100; %max(variable.Data(:));
        params.cmin = 0;% min(variable.Data(:));
        params.autoscale = 'True';
        params.save_movie = nam.txt;                                       % Output a movie and name it CP or NOCP
        plot_fields(variable_u, variable_v, variable, variable_qc, nam, params)
    end

    %% -- Plot Budget Terms --

    if plot_out == 1
        plot_budget_terms(file, folder, nam, output, i)
    end

    %% -- Compute and plot LWP --

    if plot_out == 2
        params.output = output;
        params.time = time;
        params.export = 'mov';                                             % export as 'frames' or 'movie'.
        params.cmap = cmocean('dense');
        plot_lwp(file, folder, nam, params)
    end

    %% -- Compute and Plot Vertical Profiles --

    if plot_out == 3

        % Loops over the selected folders, one line for each dataset
        for j = 1:size(folders, 2)

            % Import the data from currently selected folder
            files_all = dir(strcat(scratch,folders(j),'wrfout*'));
            fprintf('Filename: %s\n', files_all(i).name);
            file = files_all(i).name;
            fname=strcat(scratch,folders(j),file);

            % Call the 1D profiles script to compute the profiles
            [Z, U(j,:), V(j, :), W(j,:), TH(j,:), QT(j,:), QC(j,:), QR(j,:),...
                TKE(j,:), TKE_HOR(j,:), TKE_W(j,:), TEMP(j,:), CLD_FRC(j, :)...
                TOT_WAT(j,:)] = oned_profiles(fname, nam);
        end

        %% Setup the plots for the profiles
        params = {U, V, W, TKE, TKE_HOR, TH, TEMP, QT*1000, QC*1000, ...
            QR*1000, CLD_FRC, TOT_WAT, Z}; % converted Qt and Qc to [g/kg] from [kg/kg]

        xlimits = zeros(10, 2);
        xlimits(1,:) = [-1.5, 0.5];
        xlimits(2,:) = [-0.1, 0.1];
        xlimits(3,:) = [-0.001, 0.001];
        xlimits(4,:) = [-0.1, 0.1];
        xlimits(5,:) = [-0.1, 0.5];
        xlimits(6,:) = [290,360];
        xlimits(7,:) = [190, 300];
        xlimits(8,:) = [0, 18];
        xlimits(9,:) = [0, 0.05];
        xlimits(10,:) = [0, 0.025];
        xlimits (11,:) = [0, 0.01];
        xlimits(12,:) = [0, 10^7];

        xlabels = {"$\langle u \rangle$ (ms$^{-1}$)", ...
            "$\langle v \rangle$ (ms$^{-1}$)", ...
            "$\langle w \rangle$ (ms$^{-1}$)", ...
            "TKE (m$^2$s$^{-2}$)", ...
            "$1/2(u'^2+v'^2)$ (m$^2$s$^{-2}$)", ...
            "$\langle \theta \rangle$ (K)", ...
            "T (K)", ...
            "$\langle q_t \rangle$ (g/kg)", ...
            "$\langle q_c \rangle$ (g/kg)", ...
            "$\langle q_r \rangle$ (g/kg)", ...
            "Cloud Fraction", ...
            "$T_{\rm water}$[kg]"};

        titles = {"$\langle U \rangle$", ...
            "$\langle V \rangle$", ...
            "$\langle W \rangle$", ...
            "TKE", ...
            "TKE Horizontal", ...
            "$\langle \theta \rangle$", ...
            "$\langle T \rangle $", ...
            "$\langle q_t \rangle $", ...
            "$\langle q_c \rangle$", ...
            "$\langle q_r \rangle$", ...
            "Cloud Fraction", ...
            "Total Water"};

        legendLabels = {"OLD GATE NOCP", "GATE NOCP 100km"};

        %% Plot the vertical profiles with lines for each dataset
        plot_1d_profs(params, xlabels, titles, legendLabels, xlimits, time)
    end


    %% -- Plot out the vertical structure --
    %
    if plot_out == 4  % Output the vertical structure
        figure()
        for j = 1:size(folders, 2)

            % Get a list of all files in the current directory
            files_all = dir(strcat(folders(j), 'wrfout*'));
            file = files_all(i).name;
            fname=strcat(folder,file);

            % Preallocate arrays
            nam.levs = size(ncread(strcat(folders(j),file),'U'), 3);
            Z = zeros(2, nam.levs);
            P = zeros(2, nam.levs);
            H = zeros(2, nam.levs);

            % Import the data
            cp = cold_pools(j);
            fname=strcat(folders(j),file);
            [Z(j, :), P(j,:), H(j,:)] = vert_struct(fname, nam);

            %% Plot the vertical structure
            subplot(1, 2, j)
            plot(P(j,:)/10^5, Z(j,:)/10^3, '-o', 'LineWidth', 2); grid on;
            xlabel('P [bar]')
            ylabel('z [km]')
            set(gca, 'XScale', 'log')
            ylim([0, 4])
            title(cold_pools(j))

        end
    end

    %% -- Plot Spectra --

    if plot_out == 5

        % Computes 2D kinetic energy spectra at specified height, z
        z = 1;
        plot_KE_spectra(fname, z)
    end

    % Adds in exceptions handling for the old large domain simulations.
    %     catch
    %         fprintf('File %s failed\n',file)
    %     end
end




