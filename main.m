% Jacks Package for analysing the LES data in MATLAB

clear variables 

txt = 'CP';
%folder = strcat('./', txt, '/');
folder = './new/TEST/NOCP_OUT/'

output = './plots/';

% Get a list of all files in the current directory                         
files_all = dir(strcat(folder, 'wrf*'));

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
nam.dt = 4.0;                                                              % Output frequency [hours]

% For creating movies frame by frame 
for i = 1:length(files_all)
    
%     try 

    fprintf('Filename: %s\n', files_all(i).name);
    file = files_all(i).name;
    fname=strcat(folder,file);

% Make 2D Movie plots  
    [variable] = loadNetCDF(fname, 'CLDFRA'); 

    test = max(variable.Data)
    % Make plots with params 
    params.cmax = 0.01 
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

% Plot Spectra 
%         two_point_cor(file, folder, nam, output, i)
    
%% Compute and plot 1D profiles 
%    [Z, U,TH, QT, QC, TKE, TKE_HOR, TKE_W, WAT_FLUX] = oned_profiles(file, nam, output, 'NOCP'); 
%    [Z_CP, U_CP,TH_CP, QT_CP, QC_CP, TKE_CP, TKE_HOR_CP, TKE_W_CP, WAT_FLUX_CP] = oned_profiles(file, nam, output, 'CP');
%    plot_oned_profiles(Z, U,TH, QT, QC, TKE, TKE_HOR, TKE_W, WAT_FLUX, Z_CP, U_CP,TH_CP, QT_CP, QC_CP, TKE_CP, TKE_HOR_CP, TKE_W_CP, WAT_FLUX_CP)
%     
%      catch
%          fprintf('File %s failed\n',file)
%      end


end
%%

