% Jacks Package for analysing the LES data in MATLAB

clear variables

txt = 'NOCP';

folder = './data/large_domain/CP_OUT/'
%folder = '/scratch/05999/mkurowsk/ocean_nocp/';

% Which plot diagnostic to use: 
plot_out = 3; % 0: Fields, 1: Budgets, 2:LWP, 3:1D profiles, 4: Vertical Structure

output = './plots/';

% Get a list of all files in the current directory
files_all = dir(strcat(folder, 'wrfout*'));
sample_file = strcat(folder,files_all(1).name);                            % File for loading in all the sim parameters

nam.R=287.04;                                                              % [J/kg/K]
nam.cp=1004.67;                                                            % [J/kg/K]
nam.g=9.81;                                                                % [m/s2]
nam.Ll=2.50e6;                                                             % latent heat of evaporation (vapor:liquid) at 0C [J/kg]
nam.Li=2.83e6;                                                             % latent heat of sublimation (vapor:solid) at 0C [J/kg]
nam.T0=300;                                                                % ncread(fname,'T00'); % base state temperature [K]
nam.P0=1.e5;                                                               % ncread(fname,'P00'); % base state pressure [Pa]
nam.dx=double(ncreadatt(sample_file,'/','DX'));                            % [m]
nam.dy=double(ncreadatt(sample_file,'/','DY'));                            % [m]
nam.txt = txt;                                                             % A label for the data
nam.dt = 2.0;                                                              % Output frequency [hours]
nam.levs = size(ncread(sample_file,'U'), 3);                               % Number of vertical levels in the simulation
nam.nx = size(ncread(sample_file,'U'), 1);                                 % Number of x grid points in simulation
nam.ny = size(ncread(sample_file,'U'), 2);                                 % Number of y grid points in simulation

% For creating movies frame by frame
for i = 4:4 %length(files_all)

    %     try
    %% -- Compute and plot files frame by frame into a movie --

    if plot_out == 0

        fprintf('Filename: %s\n', files_all(i).name);
        file = files_all(i).name;
        fname=strcat(folder,file);

        % Load the data
        [variable] = loadNetCDF(fname, 'HFX');

        % Make plots with params
        params.cmax = 40;
        params.cmin = 0;
        params.time = (i-1)*nam.dt;
        params.name = file;
        params.save_folder = append(output, variable.Name, '/');
        params.save_movie = txt;                                               % Output a movie and name it CP or NOCP
        plot_fields(variable, nam, params)
    end

    %% -- Plot Budget Terms --

    if plot_out == 1
        plot_budget_terms(file, folder, nam, output, i)
    end

    %% -- Compute and plot LWP --

    if plot_out == 2
        plot_lwp(file, folder, nam, output, i)
    end

    %% -- Compute and plot 1D profiles --

    cold_pools = ["Large Domain, NOCP", "Large Domain, CP"];
    folders = ["./data/large_domain/NOCP_OUT/", "./data/large_domain/CP_OUT/"];

    if plot_out == 3
        for j = 1:size(folders, 2) % Loop over the profile data
            % Import the data from folder 
            files_all = dir(strcat(folders(j), 'wrfout*'));
            fprintf('Filename: %s\n', files_all(i).name);
            file = files_all(i).name;

            cp = cold_pools(j);
            fname=strcat(folders(j),file);
            [Z, U(j,:),TH(j,:), QT(j,:), QC(j,:), TKE(j,:), TKE_HOR(j,:), ...
                TKE_W(j,:), LWP(j,:)] = oned_profiles(fname, nam);
        end

        params = {U, TH, QT, QC, TKE, TKE_HOR, TKE_W, LWP, Z};
        xlabels = {"u (ms^{-1})", "\theta (K)", "q_t (gkg^{-1})", "q_c (gkg^{-1})", ...
            "TKE (m^2s^{-2})", "1/2(u'^2+v'^2) (m^2s^{-2})", "1/2(w'^2) (m^2s^{-2})", ...
            "kgm^{-2}"};
        titles = {"U", "TH", "QT", "QC", "TKE", "TKE HOR", "TKE W", "LWP"}; 
        legendLabels = {'CP', 'NO CP'};
        time = i*nam.dt;
        %%
        plot_1d_profs(params, xlabels, titles, legendLabels, time)
    end


    %% -- Plot out the vertical structure --
    %
    if plot_out == 4  % Output the vertical structures
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

            % Plot the vertical structure
            subplot(1, 2, j)
            plot(P(j,:)/10^5, Z(j,:)/10^3, '-o', 'LineWidth', 2); grid on;
            xlabel('P [bar]')
            ylabel('z [km]')
            set(gca, 'XScale', 'log')
            set(gca, 'XDir','reverse')
            ylim([0, 4])
            title(cold_pools(j))

        end
    end

    %     catch
    %         fprintf('File %s failed\n',file)
    %     end
end




