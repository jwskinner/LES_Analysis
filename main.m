% LES Analysis for studying LES data in MATLAB
% This main script loops over the data and performs analysis frame by frame
% there are multiple options (e.g., Fields, Budgets, LWP, 1D profiles,
% Vertical Structure, etc.) that can be chosen with the plot_out parameter.

clear variables

%scratch = "/scratch/05999/mkurowsk/"; % For Tacc 
scratch = "./data/"; % For jacks laptop

folders = ["./large_domain/CP_OUT/", "./large_domain/NOCP_OUT/"]
%folders = ["./GATE_NOCP_int_domain/", "./GATE_CP_int_domain/"];
%folders = ["ocean_cp/", "ocean_nocp/"]; 
%folders = ["GATE_NOCP_CONSTFLX/", "GATE_NOCP_CONSTFLX_100km/", "GATE_NOCP_CONSTFLX_50km/"];
 

% Takes the first folder for loading in params. 
folder = strcat(scratch, folders(2));  %'./data/large_domain/CP_OUT/' %'/scratch/05999/mkurowsk/ocean_nocp/' on TACC;
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
plot_out = 3;  % 0: Fields, 1: Budgets, 2:LWP, 3:1D profiles, 4: Vertical Structure, 5: KE Spectra

% For creating movies frame by frame
for i = 16:16 %length(files_all)

    file = files_all(i).name;
    fname=strcat(folder,file);

    % Sets up the paramters for the plotting
    fprintf('Filename: %s\n', files_all(i).name);

    time = (i-1)*nam.dt; 

    %% -- Compute and plot files frame by frame into a movie --
    if plot_out == 0
  
        params.time = time;
        params.name = file;
        params.save_folder = append(output, variable.Name, '/');

        % Load the data
        [variable] = loadNetCDF(fname, 'HFX');

        % Paramters specific to this plot 
        params.cmax = 40;
        params.cmin = 0;
        params.save_movie = txt;                                           % Output a movie and name it CP or NOCP
        plot_fields(variable, nam, params)
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
                TKE(j,:), TKE_HOR(j,:), TKE_W(j,:), TEMP(j,:), TOT_WAT(j,:)] ... 
                = oned_profiles(fname, nam);
        end

        %% Setup the plots for the profiles
        params = {U, V, W, TKE, TKE_HOR, TH, TEMP, QT*1000, QC*1000, QR*1000, TOT_WAT, Z}; % converted Qt and Qc to [g/kg] from [kg/kg]
        
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
            "Total Water"};

        legendLabels = {"OLD GATE NOCP", "GATE NOCP 100km"};

        %% Plot the vertical profiles with lines for each dataset
        plot_1d_profs(params, xlabels, titles, legendLabels, time)
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




