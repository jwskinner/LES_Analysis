
% This is a simpler script for plotting the CKE budgets from kardunov and
% randall paper

addpath('../cmocean-main/')

% folder = "/Users/jwskinner/Desktop/LES Analysis/data/large_domain/NOCP_OUT/"
folder = "/scratch/05999/mkurowsk/GATE_CP_CONSTFLX/"


files = [13, 49, 97];
%files = [96,118];
%files = [4, 7, 13, 16]
time = [6, 24, 48];
%time = [48, 58];

files_all = dir(strcat(folder, 'wrfout*'));

% Get the screen size
scrsz = get(0,'ScreenSize');

figure1 = figure('Position', [0 0 scrsz(3)*0.5 scrsz(3)*0.4]);

for i = 1:length(files)

    file_i = files_all(files(i)).name;
    fname=strcat(folder,file_i)

    % Import the constants and vertical structure variables
    [nam] = get_constants(fname);

    %% ========================================================================
    % read all relevant 3-D variables

    ph =ncread(fname,'PH' );                                                   % geopotential perturbation [m2/s2]
    phb=ncread(fname,'PHB');                                                   % base geopotential [m2/s2)
    p  =ncread(fname,'P'  );                                                   % pressure perturbation [Pa]
    pb =ncread(fname,'PB' );                                                   % base pressure [Pa]
    th =ncread(fname,'T'  )+nam.T0;                                            % Potential temperature [K]

    qv=ncread(fname,'QVAPOR');                                                 % Water vapor mixing ratio [kg/kg]
    qc=ncread(fname,'QCLOUD');                                                 % Cloud water mixing ratio [kg/kg]
    qr=ncread(fname,'QRAIN');                                                  % Rain water mixing ratio [kg/kg]
    qi=ncread(fname,'QICE');                                                   % Ice mixing ratio [kg/kg]

    %%

    %% ========================================================================
    % initial calculations

    s=size(th);
    n=s(1);
    m=s(2);
    l=s(3);
    nm=n*m;

    HS=mean(reshape(ph+phb,nm,l+1));
    H=0.5*(HS(:,1:end-1)+HS(:,2:end));
    Z=H./nam.g';                                                               % height at mass-levels [m]
    p=p+pb;                                                                    % pressure
    exn=(p/nam.P0).^(nam.R/nam.cp);                                            % exner function
    qt=qv+qc+qi;                                                               % total water mixing ratio [kg/kg] (no precipitating elements)

    ql = qr + qc; %liquid water in the system

    t=th.*exn;                                                                 % temperature
    tv=t.*(1+0.608*qv);                                                        % virtual temperature, bouyancy is tv - ql (eq. 1 Marcin)
    rho=p./(nam.R*tv);                                                         % density

    % Total water LWP
    %qtrho = qt.*rho;
    qtrho = ql.*rho; % Just the liquids

    LWP = trapz(Z',qtrho,3);                                                   % vertically integrated LWP

    x_km = linspace(0, (n*nam.dx)/1000, n);
    y_km = linspace(0, (n*nam.dy)/1000, n);

    %% Plot Liquid Water Path

    subplot(1, 3, i); % add this line to specify the subplot

    cmap = cmocean('rain');

    imagesc(x_km, y_km, LWP); % use imagesc instead of image to use the logarithmic scale
    set(gca, 'YDir', 'normal'); % flip the y-axis to match the image orientation
    ylabel('$y$ [km]','LineWidth',1,'FontSize',14, 'Interpreter','latex');
    xlabel('$x$ [km]','LineWidth',1,'FontSize',14, 'Interpreter','latex');
    title(num2str(time(i),'%.1f')+" hours")
    colormap(cmap)

    % Create colorbar
    c = colorbar;
    caxis([0 1]); % set the color limits
    c.Label.String = 'LWP'; % add a label to the colorbar

    xlim([min(x_km) max(x_km)])
    ylim([min(y_km) max(y_km)])
    set(gca, 'DataAspectRatio', [1 1 1]);

    grid off

end