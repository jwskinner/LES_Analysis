% Script for computing moisture variance timeseries and timeseries of the
% budget terms (to be added later on)
%
% [Caution: This is an early draft script!]
%
% Updated by J. W. Skinner (2023-28-02)

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
cold_pools = ["./data/small_domain/CP_OUT/", "./data/small_domain/NOCP_OUT/"];
%cold_pools = ["/scratch/05999/mkurowsk/GATE_CP_CONSTFLX/", "/scratch/05999/mkurowsk/GATE_NOCP_CONSTFLX/"]

% Get a list of all files in each folder
files_cp = dir(strcat(cold_pools{1}, 'wrfout*'));
files_nocp = dir(strcat(cold_pools{2}, 'wrfout*'));

files_all = files_cp; % choose the folder with shortest output

% Preallocate arrays for storing the computed values; index 1 is for CP or
% NOCP
num_files = length(files_nocp);

turb_out = zeros(2, num_files);
micro_out = zeros(2, num_files);
prod_out = zeros(2, num_files);
diss_out = zeros(2, num_files);
vert_av_qt = zeros(2, num_files);

time_hours = zeros(1, num_files);

% Import and calculate vertical structure variables for one of the files
[Z, p, H] = vert_struct(strcat(cold_pools(1),files_cp(1).name), nam);

% Loop over the files and cases
for i = 1:num_files

    for j = 1:2

        cp = cold_pools(j);
        %         try
        fprintf('Filename: %s\n', files_all(i).name);

        file = files_all(i).name;
        fname=strcat(cp,file);

        %% ========================================================================
        % initial calculations

        [u, v, w] = loadNetCDF(fname, 'U', 'V', 'W');
        u = u.Data; v = v.Data; w = w.Data;

        u=0.5*(u(1:end-1,:,:,:)+u(2:end,:,:,:)); %C->A grid
        v=0.5*(v(:,1:end-1,:,:)+v(:,2:end,:,:));
        w=0.5*(w(:,:,1:end-1,:)+w(:,:,2:end,:));

        ph =ncread(fname,'PH' );                                           % geopotential perturbation [m2/s2]
        phb=ncread(fname,'PHB');                                           % base geopotential [m2/s2]
        tke=ncread(fname,'TKE');                                           % TKE [m2/s2]
        p  =ncread(fname,'P'  );                                           % pressure perturbation [Pa]
        pb =ncread(fname,'PB' );                                           % base pressure [Pa]
        th =ncread(fname,'T'  )+nam.T0;                                    % Potential temperature [K]

        qv=ncread(fname,'QVAPOR');                                         % Water vapor mixing ratio [kg/kg]
        qc=ncread(fname,'QCLOUD');                                         % Cloud water mixing ratio [kg/kg]
        qr=ncread(fname,'QRAIN');                                          % Rain water mixing ratio [kg/kg]
        qi=ncread(fname,'QICE');                                           % Ice mixing ratio [kg/kg]
        qs=ncread(fname,'QSNOW');                                          % Snow mixing ratio
        qg=ncread(fname,'QGRAUP');                                         % Graupel mixing ratio

        kh=ncread(fname,'XKHH');                                           % eddy-diffusivity for scalars [m2/s]

        qvtend=ncread(fname,'QVTEND_MICRO');                               % qv microphysical tendency [kg/kg 1/s]
        qctend=ncread(fname,'QCTEND_MICRO');                               % qc microphysical tendency [kg/kg 1/s]
        qitend=ncread(fname,'QITEND_MICRO');                               % qi microphysical tendency [kg/kg 1/s]
        thtend=ncread(fname,'THTEND_MICRO');                               % th microphysical tendency [K/s]

        s=size(w);
        n=s(1);
        m=s(2);
        l=s(3);
        nm=n*m;

        HS=mean(reshape(ph+phb,nm,l+1));
        H=0.5*(HS(:,1:end-1)+HS(:,2:end));

        ZS=HS./nam.g;                                                      % height at w-levels [m]
        Z=H./nam.g';                                                       % height at mass-levels [m]
        p=p+pb;                                                            % pressure
        exn=(p/nam.P0).^(nam.R/nam.cp);                                    % exner function
        qt=qv+qc+qi;                                                       % total water mixing ratio [kg/kg] (no precipitating elements)
        qt = qt * 1000; % convert to [g/kg]

        TH=mean(reshape(th,nm,l));
        t=th.*exn;                                                         % temperature
        tv=t.*(1+0.608*qv);                                                % virtual temperature, bouyancy is tv - ql (eq. 1 Marcin)
        thil=th-(nam.Ll*qc+nam.Li*qi)./(nam.cp*exn);                       % liquid water potential temperature

        rho=p./(nam.R*tv);
        RHO=mean(reshape(thil,nm,l));

        %% ========================================================================
        % QT variance budget terms -- Eq (1) Thomas et al. 2021

        W=mean(reshape(w,nm,l));
        wpr(1,1,:)=W;
        wpert=bsxfun(@minus,w,wpr);

        QT=mean(reshape(qt,nm,l));                                         % Horizontally average moisture variance
        qtpr(1,1,:)=QT;
        qtpert=bsxfun(@minus,qt,qtpr);

        % Caculate the gradient source terms in budget equation
        wqt  = wpert.*qtpert;
        wqt2 = wpert.*(qtpert.^2);

        WQT =mean(reshape(wqt ,nm,l));
        WQT2=mean(reshape(wqt2,nm,l));

        QTVAR =mean(reshape(qtpert.^2 ,nm,l));

        Qqt=qvtend+qctend+qitend;                                          % Microphysical tendency term for sources

        [dqtdx,dqtdy,dqtdz]=gradient((qtpert),nam.dx,nam.dy,50);

        % all budget terms converted to units of [1/s*(g/kg)^2]
        qtvar.prod=-2*WQT'.*gradient(QT',Z)*1e6;                           % Gradient Production Term
        qtvar.trns=-1./RHO.*gradient(RHO.*WQT2,Z)*1e6;                     % Turbulent transport term

        qtvar.srcs=2*squeeze(mean(mean((qtpert).*Qqt)))*1e6;               % Microphysical sources and sinks of the moisture variance
        qtvar.diss=-2*squeeze(mean(mean(kh.*(dqtdx.^2+dqtdy.^2+dqtdz.^2))))*1e6;   % Eq. (2) of Schemann 2016

        % Vertically averaged budget terms

        vert_av_qt(j, i) = mean(QT);
        turb_out(j, i) = mean(qtvar.trns);
        micro_out(j, i) = mean(qtvar.srcs);
        prod_out(j, i) = mean(qtvar.prod);
        diss_out(j, i) = mean(qtvar.diss);


    end
    time_hours(i) = 0 + (i-1)*nam.dt;

end


%%
figure();

subplot(1,2,1);
plot(time_hours, vert_av_qt(1, :), 'Linewidth', 1.5); hold on;
plot(time_hours, vert_av_qt(2, :), 'Linewidth', 1.5);
xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel('q_t variance (vert. avg) ([g/kg])^2','LineWidth',1.5,'FontSize',15);
legend('CP', 'NOCP','Location', 'northwest');
title('Vert. average q_t variance', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');

subplot(1,2,2);

colors = {'blue', [0 0.5 0], 'red', [1 0.5 0]};

plot(time_hours, turb_out(1,:), 'Linewidth', 1.5, 'Color', colors{1}); hold on;
plot(time_hours, turb_out(2,:), 'Linewidth', 1.5, 'Color', colors{1}, 'linestyle','--');
plot(time_hours, micro_out(1,:), 'Linewidth', 1.5, 'Color', colors{2}); hold on;
plot(time_hours, micro_out(2,:), 'Linewidth', 1.5, 'Color', colors{2}, 'linestyle','--');
plot(time_hours, prod_out(1,:), 'Linewidth', 1.5, 'Color', colors{3}); hold on;
plot(time_hours, prod_out(2,:), 'Linewidth', 1.5, 'Color', colors{3}, 'linestyle','--');
plot(time_hours, diss_out(1,:), 'Linewidth', 1.5, 'Color', colors{4}); hold on;
plot(time_hours, diss_out(2,:), 'Linewidth', 1.5, 'Color', colors{4}, 'linestyle','--');

xlabel('Time [hours]','LineWidth',1.5,'FontSize',15);
ylabel('Magnitude [1/s*(g/kg)^2]','LineWidth',1.5,'FontSize',15);
legend('TURB. CP', 'TURB. NOCP', 'MICRO. CP', 'MICRO. NOCP', 'PROD. CP', 'PROD. NOCP', 'DISS. CP', 'DISS. NOCP');
title('Vert. Avg. Budget Terms', 'LineWidth',1,'FontSize',13, 'FontWeight','Normal');




