% Script for computing timeseries properties of the convection data
% Updated by J. W. Skinner (2023-30-03)

% Define the folders for each case
% folders = ["./data/small_domain/CP_OUT/", "./data/small_domain/NOCP_OUT/"];
%folders = ["/scratch/05999/mkurowsk/GATE_NOCP_CONSTFLX/", "/scratch/05999/mkurowsk/GATE_NOCP_CONSTFLX_100km/", "/scratch/05999/mkurowsk/GATE_NOCP_CONSTFLX_50km/"]
legend_text = ["NOEVP1", "NOEVP1.3", "GATE NOCP CONSTFLX"];
folders = ["/scratch/05999/mkurowsk/GATE_NOEVP1km_CONSTFLX_50km/","/scratch/05999/mkurowsk/GATE_NOEVP1.3km_CONSTFLX_50km/", "/scratch/05999/mkurowsk/GATE_NOCP_CONSTFLX_50km/"]; 

% Get a list of all files in each folder
files_cp = dir(strcat(folders{1}, 'wrfout*'));
files_nocp = dir(strcat(folders{2}, 'wrfout*'));

files_all = files_cp; % choose the folder with shortest output
num_files = length(files_cp);

% Define constants
nam = get_constants(strcat(folders(1),files_all(1).name));

t_length = 100 %num_files; % Length of the timeseries
f_length = length(folders);  % Number of datasets to loop over 

% Preallocate 3D array for storing the computed values
data = zeros(f_length, t_length, 15);

% Pre-assign each variable to a slice of the 3D array (this is more memory
% efficient than predefining all of these seperately)
TKE_out = data(:,:,1);
LWP_out = data(:,:,2);
RAINNC_out = data(:,:,3);
PRECR_out = data(:,:,4);
PRECG_out = data(:,:,5);
CLDFR_out = data(:,:,6);
HFX_out = data(:,:,7);
QT_out = data(:,:,8);
LH_out = data(:,:,9);
CLWP_out = data(:,:,10);
RWP_out = data(:,:,11);
CloudFrac_out = data(:,:,12);
maxu_out = data(:,:, 13);
maxv_out = data(:,:,14); 
maxw_out = data(:,:,15);

time_hours = zeros(1, t_length);

% Import and calculate vertical structure variables for one of the files
[Z, p, H] = vert_struct(strcat(folders(1),files_cp(1).name), nam);

% Loop over the files and cases
parfor i = 1:t_length
    
    for j = 1:f_length
        
        cp = folders(j);
        fprintf('Filename: %s\n', files_all(i).name);
        file = files_all(i).name;
        fname=strcat(cp,file);

        % Read in the things we want, make a combined ncread call for speed 
        rainnc = ncread(fname, 'RAINNC');
        tke = ncread(fname, 'TKE');
        qv = ncread(fname, 'QVAPOR');
        qc = ncread(fname, 'QCLOUD');
        qr = ncread(fname, 'QRAIN');
        qi = ncread(fname, 'QICE');
        qs = ncread(fname, 'QSNOW');
        qg = ncread(fname, 'QGRAUP');   
        th =ncread(fname,'T'  )+nam.T0;                                    % Potential temperature [K]
        u = ncread(fname, 'U');
        v = ncread(fname, 'V');
        w = ncread(fname, 'W');

        s=size(tke); n=s(1); m=s(2); l=s(3); nm=n*m;
        exn=(p/nam.P0).^(nam.R/nam.cp);                                    % exner function
        qt=qv+qc+qi;                                                       % total water mixing ratio [kg/kg] (no precipitating elements)
        
        % Compute the temperatures 
        t = zeros(size(th)); 
        tv = zeros(size(th));
        rho = zeros(size(th)); 
        for k=1:size(exn)                                                  
            t(:,:,k) = th(:,:,k).*exn(k);                                  % temperature [k]
            tv(:,:,k) = t(:,:,k).*(1+0.608*qv(:,:,k));                     % virtual temperature, bouyancy is tv - ql (eq. 1 Marcin)
            rho=p(k)./(nam.R*tv(:,:,k));                                   % density
        end                                    
           
        % Compute water paths 
        LWP = trapz(Z',qt.*rho,3);                                         % Liquid water path total                                     
        CLWP = trapz(Z',qc*1000.*rho,3);                                   % Cloud liquid water path converted to [g/m^-2]
        RWP = trapz(Z',qr*1000.*rho,3);                                    % Rain water path converted to [g/m^-2]

        % Compute cloud cover
        qcloud_threshold = 0.001;                                          % Define cloud mixing ratio threshold

        % Create a Cloud Mask 
        Cmask = qc >= qcloud_threshold;
        nonzero_ratio = nnz(Cmask) / numel(Cmask);                         % Fraction of grid cells that contain clouds
        
        meanCmask = mean(Cmask, 3);                                        % Vertically mean the Cloud Mask to get a 'top down' look
        cloudCover = nnz(meanCmask) / numel(meanCmask);                    % Cloud cover of the domain looking top down

        % cloud_cover_fraction = nnz(mean(qc, 3) >= qcloud_threshold) / numel(qc); % fraction over the vertical mean

        % Compute mean values adn output from the loop  
        TKE = mean(reshape(tke,nm,l));

        TKE_out(j, i) = trapz(Z,TKE);
        RAINNC_out(j, i) = mean(rainnc, 'all');
        LWP_out(j, i) = mean(LWP, 'all');
        CLWP_out(j, i) = mean(CLWP, 'all');
        RWP_out(j, i) = mean(RWP, 'all');
        CloudFrac_out(j, i) = cloudCover; 

        QT=mean(reshape(qt,nm,l));
        QT_out(j, i) = trapz(Z,QT);

        maxu_out(j, i) = max(u(:)); 
        maxv_out(j, i) = max(v(:));
        maxw_out(j, i) = min(w(:));
 
    end

    time_hours(i) = 0 + (i-1)*nam.dt;

end

%% Image plotting to visually check the cloud fraction thresholding. 
% figure()
% subplot(1, 2, 1)
% imagesc(mean(Cmask, 3))
% colorbar
% title("Mean Cloud Mask qc >= 0.001")
% subplot(1, 2, 2)
% imagesc(mean(qc, 3))
% title("Mean Qcloud field")
% colorbar

%% Output the LWP, CLWP and cloud fraction
figure(); 
subplot(3, 1, 1);
generate_subplot(LWP_out, time_hours,  legend_text, '', 'Liquid Water Path', folders);

% subplot(3, 1, 2);
% generate_subplot(CLWP_out, time_hours,  legend_text, '', 'Cloud Water Path', folders);

subplot(3, 1, 2);
generate_subplot(RWP_out, time_hours,  legend_text, '', 'Rain Water Path', folders);

subplot(3, 1, 3);
generate_subplot(CloudFrac_out, time_hours,  legend_text, '', 'Cloud Fraction', folders);

%%
figure();

subplot(1, 3, 1);
generate_subplot(TKE_out, time_hours,  legend_text, 'Vert. integrated TKE',...
    'VTKE [m^2 s^{-2}]', folders);

subplot(1, 3, 2);
generate_subplot(LWP_out, time_hours, legend_text, 'Vert. integrated LWP', ...
    'LWP [kg m^{-2}]', folders, [60, 70]);

subplot(1, 3, 3);
generate_subplot(RAINNC_out, time_hours, legend_text, ...
    'Accumulative precipitation', 'RAINNC [mm]', folders);

figure();

subplot(1, 3, 1);
generate_subplot(maxu_out, time_hours,  legend_text, 'Max(u)',...
    'Max(u) [m s^{-1}]', folders);

subplot(1, 3, 2);
generate_subplot(maxv_out, time_hours, legend_text, 'Max(v)', ...
    'Max(v) [m s^{-1}]', folders);

subplot(1, 3, 3);
generate_subplot(maxw_out, time_hours, legend_text, ...
    'Min(w)', 'Min(w) [m s^{-1}]', folders);

%%
function generate_subplot(data, time_hours, legend_text, title_text, ...
    y_label, folders, ylim_range)

    cmap = lines(length(folders)+1);

    for i = 1:length(folders)
        plot(time_hours, data(i, :), 'Color', cmap(i,:), 'Linewidth', 1.5); 
        hold on;
    end
    xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
    ylabel(y_label,'LineWidth',1.5,'FontSize',15);
    legend(legend_text, 'Location', 'northwest');
    title(title_text, 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
    if nargin > 6
        ylim(ylim_range);
    end
    
end


%% -- OLD CODE -- 
% figure();
% 
% subplot(1,3,1);
% for i=1:length(folders)
% plot(time_hours, TKE_out(i, :), 'Linewidth', 1.5); hold on;
% end
% xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
% ylabel('VTKE [m^2 s^{-2}]','LineWidth',1.5,'FontSize',15);
% legend(legend_text,'Location', 'northwest');
% title('Vert. integrated TKE', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
% 
% subplot(1,3,2);
% for i=1:length(folders)
% plot(time_hours, LWP_out(i,:), 'Linewidth', 1.5); hold on;
% end
% xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
% ylabel('LWP [kg m^{-2}]','LineWidth',1.5,'FontSize',15);
% legend(legend_text, 'Location', 'northwest');
% ylim([45, 60])
% title('Vert. integrated LWP', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
% 
% subplot(1,3,3);
% for i=1:length(folders)
% plot(time_hours, RAINNC_out(i, :), 'Linewidth', 1.5); hold on;
% end
% xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
% % ylabel(append('RAINNC', ' [', rainnc.Units, ']'),'LineWidth',1.5,'FontSize',15);
% ylabel(append('RAINNC', ' [mm]'),'LineWidth',1.5,'FontSize',15);
% title('Accumulative precipitation', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');
% legend(legend_text, 'Location', 'northwest')

%%
% figure()
% for i=1:length(folders)
% plot(time_hours, QT_out(i, :)*1000, 'Linewidth', 1.5); hold on;
% end
% xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
% ylabel('[(g/kg)^2]','LineWidth',1.5,'FontSize',15);
% legend(legend_text,'Location', 'northwest');
% title('Vert. avg. Q_t', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');

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