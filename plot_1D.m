
% This is a simpler and more broad script for plotting 1d profiles of any
% variabes from a single file rather than doing it in the loop over all
% files in main.m

fname = "./data/large_domain/CP_OUT/wrfout_d01_2020-01-01_22:00:00";
time = 16.5; % time in hours form the filename

% Import the constants and vertical structure variables
[nam] = get_constants(fname);

%% ========================================================================
% read all relevant 3-D variables

u=ncread(fname,'U'); u=0.5*(u(1:end-1,:,:,:)+u(2:end,:,:,:)); %C->A grid
v=ncread(fname,'V'); v=0.5*(v(:,1:end-1,:,:)+v(:,2:end,:,:));
w=ncread(fname,'W'); w=0.5*(w(:,:,1:end-1,:)+w(:,:,2:end,:));

ph =ncread(fname,'PH' );                                                   % geopotential perturbation [m2/s2]
phb=ncread(fname,'PHB');                                                   % base geopotential [m2/s2)
tke=ncread(fname,'TKE');                                                   % TKE [m2/s2]
p  =ncread(fname,'P'  );                                                   % pressure perturbation [Pa]
pb =ncread(fname,'PB' );                                                   % base pressure [Pa]
th =ncread(fname,'T'  )+nam.T0;                                            % Potential temperature [K]

qv=ncread(fname,'QVAPOR');                                                 % Water vapor mixing ratio [kg/kg]
qc=ncread(fname,'QCLOUD');                                                 % Cloud water mixing ratio [kg/kg]
qr=ncread(fname,'QRAIN');                                                  % Rain water mixing ratio [kg/kg]
qi=ncread(fname,'QICE');                                                   % Ice mixing ratio [kg/kg]
qs=ncread(fname,'QSNOW');                                                  % Snow mixing ratio
qg=ncread(fname,'QGRAUP');                                                 % Graupel mixing ratio

kh=ncread(fname,'XKHH');                                                   % eddy-diffusivity for scalars [m2/s]

qvtend=ncread(fname,'QVTEND_MICRO');                                       % qv microphysical tendency [kg/kg 1/s]
qctend=ncread(fname,'QCTEND_MICRO');                                       % qc microphysical tendency [kg/kg 1/s]
qitend=ncread(fname,'QITEND_MICRO');                                       % qi microphysical tendency [kg/kg 1/s]
thtend=ncread(fname,'THTEND_MICRO');                                       % th microphysical tendency [K/s]

%%

%% ========================================================================
% initial calculations

s=size(w);
n=s(1);
m=s(2);
l=s(3);
nm=n*m;

HS=mean(reshape(ph+phb,nm,l+1));
H=0.5*(HS(:,1:end-1)+HS(:,2:end));

ZS=HS./nam.g;                                                              % height at w-levels [m]
Z=H./nam.g';                                                               % height at mass-levels [m]
p=p+pb;                                                                    % pressure
exn=(p/nam.P0).^(nam.R/nam.cp);                                            % exner function
qt=qv+qc+qi;                                                               % total water mixing ratio [kg/kg] (no precipitating elements)

QT = mean(reshape(qt,nm,l));
QC = mean(reshape(qc,nm,l));
QR = mean(reshape(qr,nm,l));

TKE = mean(reshape(tke,nm,l));

qtpr(1,1,:)=QT;
qtpert=bsxfun(@minus,qt,qtpr);

TH=mean(reshape(th,nm,l));

t=th.*exn;                                                                 % temperature [K]
TEMP = mean(reshape(t,nm,l));

tv=t.*(1+0.608*qv);                                                        % virtual temperature, bouyancy is tv - ql (eq. 1 Marcin)
thil=th-(nam.Ll*qc+nam.Li*qi)./(nam.cp*exn);                               % liquid water potential temperature

U=mean(reshape(u,nm,l));
upr(1,1,:)=U;
upert=bsxfun(@minus,u,upr);

V=mean(reshape(v,nm,l));
vpr(1,1,:)=V;
vpert=bsxfun(@minus,v,vpr);

W=mean(reshape(w,nm,l));
wpr(1,1,:)=W;
wpert=bsxfun(@minus,w,wpr);
WPERT=mean(reshape(wpert,nm,l));
WPERT2 = mean(reshape(wpert.*wpert,nm,l));
WPERT3 = mean(reshape(wpert.*wpert.*wpert,nm,l));

TKE_hor = 0.5*(upert.^2 + vpert.^2);
TKE_HOR=mean(reshape(TKE_hor,nm,l));
TKE_W = 0.5*mean(reshape(wpert,nm,l)).^2;

% Compute cloud cover
qcloud_threshold = 0.001;                                                  % Define cloud mixing ratio threshold

% Create a Cloud Mask
for i= 1: size(qc, 3) % loop over the heights
    Cmask = qc(:,:,i) >= qcloud_threshold;
    CLD_FRC(i) = nnz(Cmask) / numel(Cmask);                                 % Fraction of grid cells that contain clouds
end

rho=p./(nam.R*tv);                                                         % density
RHO=mean(reshape(thil,nm,l));
qtrho = qt.*rho;

v_air = (nam.dx * nam.dy * (max(Z)/nam.levs));                             % Volume of each grid box [m^3]
M_AIR = RHO .* v_air;
TOT_WAT = (M_AIR .* QT);                                                   % Total water [Kg]


%% Setup the plots for the profiles
% params = {U, V, W, TKE, TKE_HOR, TH, TEMP, QT*1000, QC*1000, ...
%     QR*1000, CLD_FRC, TOT_WAT, Z}; % converted Qt and Qc to [g/kg] from [kg/kg]
%
% xlimits = zeros(10, 2);
% xlimits(1,:) = [-1.5, 0.5];
% xlimits(2,:) = [-0.1, 0.1];
% xlimits(3,:) = [-0.01, 0.01];
% xlimits(4,:) = [-0.1, 0.1];
% xlimits(5,:) = [-0.1, 0.5];
% xlimits(6,:) = [290,360];
% xlimits(7,:) = [190, 300];
% xlimits(8,:) = [0, 18];
% xlimits(9,:) = [0, 0.05];
% xlimits(10,:) = [0, 0.025];
% xlimits (11,:) = [0, 0.01];
% xlimits(12,:) = [0, 10^7];
%
% xlabels = {"$\langle u \rangle$ (ms$^{-1}$)", ...
%     "$\langle v \rangle$ (ms$^{-1}$)", ...
%     "$\langle w \rangle$ (ms$^{-1}$)", ...
%     "TKE (m$^2$s$^{-2}$)", ...
%     "$1/2(u'^2+v'^2)$ (m$^2$s$^{-2}$)", ...
%     "$\langle \theta \rangle$ (K)", ...
%     "T (K)", ...
%     "$\langle q_t \rangle$ (g/kg)", ...
%     "$\langle q_c \rangle$ (g/kg)", ...
%     "$\langle q_r \rangle$ (g/kg)", ...
%     "Cloud Fraction", ...
%     "$T_{\rm water}$[kg]"};
%
% titles = {"$\langle U \rangle$", ...
%     "$\langle V \rangle$", ...
%     "$\langle W \rangle$", ...
%     "TKE", ...
%     "TKE Horizontal", ...
%     "$\langle \theta \rangle$", ...
%     "$\langle T \rangle $", ...
%     "$\langle q_t \rangle $", ...
%     "$\langle q_c \rangle$", ...
%     "$\langle q_r \rangle$", ...
%     "Cloud Fraction", ...
%     "Total Water"};
%
legendLabels = {"OLD GATE NOCP", "GATE NOCP 100km"};
%
% plot_1d_profs(params, xlabels, titles, legendLabels, xlimits, time)


%% Setup the plots for the profiles
params = {W, WPERT, WPERT2, WPERT3}; % converted Qt and Qc to [g/kg] from [kg/kg]

xlimits = zeros(10, 2);
xlimits(1,:) = [-0.005, 0.005];
xlimits(2,:) = [-0.005, 0.005];
xlimits(3,:) = [0, 5];
xlimits(4,:) = [0, 5];

ylimits = [0,12];

xlabels = {"$\langle w \rangle$", ...
    "$\langle w' \rangle$", ...
    "$\langle w'w' \rangle$", ...
    "$\langle w'w'w' \rangle$", ...
    };

titles = xlabels; 

for i = 1:(size(params, 2)) % removes the Z variable from the sizing
    subplot(1,4,i);
    for j = 1:size(params{i}, 1) % Loops over lines to plot
        P(i) = plot(params{i}(j,:), Z/10^3, 'LineWidth', 2); ...
            hold on; grid on;
    end
    %         legend(legendLabels, 'Location', 'best')
    xlabel(xlabels{i},'Interpreter','latex')
    ylabel('$z$ [km]','Interpreter','latex')
    title(titles{i},'FontWeight','Normal','Interpreter','latex')
    xlim([xlimits(i,1), xlimits(i,2)])
    ylim(ylimits)
    legend("CP 50km (old)", "CP 200km (old)", "CP 100 km (GATE)")
end
sgtitle(strcat(num2str(time), ' hours'))
exportgraphics(gcf,strcat('./plots/profiles/', 'profiles', '.gif'),...
    'Resolution',150, 'Append',true)


