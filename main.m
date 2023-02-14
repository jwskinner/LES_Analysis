% Jacks Package for analysing the LES data in MATLAB

clear variables

txt = 'NOCP';

folder = './data/small_domain/CP_OUT/'
%folder = '/scratch/05999/mkurowsk/ocean_nocp/';

output = './plots/';

% Get a list of all files in the current directory
files_all = dir(strcat(folder, 'wrfout*'));

nam.R=287.04;                                                              % [J/kg/K]
nam.cp=1004.67;                                                            % [J/kg/K]
nam.g=9.81;                                                                % [m/s2]
nam.Ll=2.50e6;                                                             % latent heat of evaporation (vapor:liquid) at 0C [J/kg]
nam.Li=2.83e6;                                                             % latent heat of sublimation (vapor:solid) at 0C [J/kg]
nam.T0=300;                                                                % ncread(fname,'T00'); % base state temperature [K]
nam.P0=1.e5;                                                               % ncread(fname,'P00'); % base state pressure [Pa]
nam.dx=double(ncreadatt(strcat(folder,files_all(1).name),'/','DX'));       % [m]
nam.dy=double(ncreadatt(strcat(folder,files_all(1).name),'/','DY'));       % [m]
nam.txt = txt;                                                             % A label for the data
nam.dt = 0.5;                                                              % Output frequency [hours]

% For creating movies frame by frame
for i = 1:length(files_all)

    %     try

    fprintf('Filename: %s\n', files_all(i).name);
    file = files_all(i).name;
    fname=strcat(folder,file);

    % Make 2D Movie plots
        [variable] = loadNetCDF(fname, 'HFX');

%     Make plots with params
        params.cmax = 40;
        params.cmin = 0;
        params.time = (i-1)*nam.dt;
        params.name = file;
        params.save_folder = append(output, variable.Name, '/');
        params.save_movie = txt;                                               % Output a movie and name it CP or NOCP
        plot_fields(variable, nam, params)


    % Plot budget Terms
    %         plot_budget_terms(file, folder, nam, output, i)

    % Plot Liquid Water Path
    %         plot_lwp(file, folder, nam, output, i)

    %% Compute and plot 1D profiles
% 
%     cold_pools = ["CP", "NOCP"];
%     folders = ["./data/small_domain/CP_OUT/", "./data/small_domain/NOCP_OUT/"];
% 
%     for j = 1:2
%         cp = cold_pools(j);
%         fname=strcat(folders(j),file);
%         [Z, U(j,:),TH(j,:), QT(j,:), QC(j,:), TKE(j,:), TKE_HOR(j,:), ...
%             TKE_W(j,:), WAT_FLUX(j,:)] = oned_profiles(fname, nam);
%     end
% 
%     params = {U, TH, QT, QC, TKE, TKE_HOR, TKE_W, WAT_FLUX, Z};
%     xlabels = {"u (ms^{-1})", "\theta (K)", "q_t (gkg^{-1})", "q_c (gkg^{-1})", ...
%         "TKE (m^2s^{-2})", "1/2(u'^2+v'^2) (m^2s^{-2})", "1/2(w'^2) (m^2s^{-2})", ...
%         "\times 10^{-2} (w'q_t') (mgs^{-1}kg^{-1})"};
%     legendLabels = {'CP', 'NO CP'};
% 
%     plot_1d_profs(params, xlabels, legendLabels)


    %      catch
    %          fprintf('File %s failed\n',file)
    %      end


end
%%

