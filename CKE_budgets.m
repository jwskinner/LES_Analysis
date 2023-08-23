
% This is a simpler script for plotting the CKE budgets from kardunov and
% randall paper 

 fname = "/Users/jwskinner/Desktop/LES Analysis/data/small_domain/NOCP_OUT/wrfout_d01_0001-01-01_15:30:00";
%fname = "/Volumes/Fortress L3/WRF_Data/int_domain/GATE_NOCP/wrfout_d01_2020-01-01_16:00:00";

time = 16.5; % time in hours form the filename

legendtext = ["CP", "NOCP"];

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
thpr(1,1,:)=TH;
thpert=bsxfun(@minus,th,thpr);

t=th.*exn;                                                                 % temperature [K]
TEMP = mean(reshape(t,nm,l));

tv=t.*(1+0.608*qv);                                                        % virtual temperature, bouyancy is tv - ql (eq. 1 Marcin)
thil=th-(nam.Ll*qc+nam.Li*qi)./(nam.cp*exn);                               % liquid water potential temperature

THIL=mean(reshape(thil,nm,l));
thilpr(1,1,:)=THIL;
thilpert=bsxfun(@minus,thil,thilpr);

U=mean(reshape(u,nm,l));
upr(1,1,:)=U;
upert=bsxfun(@minus,u,upr);

V=mean(reshape(v,nm,l));
vpr(1,1,:)=V;
vpert=bsxfun(@minus,v,vpr);

% vertical velocity terms
W=mean(reshape(w,nm,l));
wpr(1,1,:)=W;
wpert=bsxfun(@minus,w,wpr);
WPERT=mean(reshape(wpert,nm,l));
WPERT2 = mean(reshape(wpert.*wpert,nm,l));
WPERT3 = mean(reshape(wpert.*wpert.*wpert,nm,l));

% TKE terms 
TKE_hor = 0.5*(upert.^2 + vpert.^2);
TKE_HOR=mean(reshape(TKE_hor,nm,l));
TKE_W = 0.5*mean(reshape(wpert,nm,l)).^2;


% Pressure 
P=mean(reshape(p,nm,l));
ppr(1,1,:)=P;
ppert=bsxfun(@minus,p,ppr);

% Density 
rho=p./(nam.R*tv);              
RHO=mean(reshape(rho,nm,l));

% Virtual potential temperature 
TV=mean(reshape(tv,nm,l));
tvpr(1,1,:)=TV;
tvpert=bsxfun(@minus,tv,tvpr);

% CTKE terms 
cke = 0.5*rho.*(upert.^2 + vpert.^2+wpert.^2);
CKE=mean(reshape(cke,nm,l));
ckepr(1,1,:)=CKE;
cke_pert=bsxfun(@minus,cke,ckepr);

%% ========================================================================
% CKE Budget terms from Equation (6) of Khairoutdinov and Randall 2001

% Transport term 
WE_pert = mean(reshape(wpert .* cke_pert,nm,l));
dWE_dz = gradient(WE_pert, Z);
CTKE_T_Term = -dWE_dz;

% Pressure Correlation term 
WP = mean(reshape(wpert.*ppert,nm,l));
dWP_dz = gradient(WP, Z);
CTKE_P_Term = - dWP_dz; 

% Mean vertical wind shear term
UW_pert = mean(reshape(upert .* wpert,nm,l));
du_dz = gradient(U, Z);
VW_pert = mean(reshape(vpert .* wpert,nm,l));
dv_dz = gradient(V, Z);
CTKE_S_Term = -RHO .* (UW_pert.*du_dz + VW_pert.*dv_dz);


% Bouyancy production term 
Beta = nam.g ./ TH; % JACK double check this term 

B_flx = mean(reshape(wpert.*tvpert,nm,l));  % buoyancy flux
CTKE_B_Term = -RHO.*Beta.*B_flx;

% Viscous dissipation term 
[dthildx,dthildy,dthildz]=gradient((thilpert),nam.dx,nam.dy,50);
D_Term=-2*squeeze(mean(mean(kh.*(dthildx.^2+dthildy.^2+dthildz.^2))));


% Residuals 

CTKE_Res = CTKE_T_Term + CTKE_P_Term + CTKE_S_Term + CTKE_B_Term + D_Term; 

%% Fluxes 
M_flx = mean(reshape(qtpert .* wpert,nm,l)); % Moisture flux  

%% ========================================================================
% TKE Budget terms from Stull 1988

TKE_e = 0.5*(upert.^2 + vpert.^2 +wpert.^2);

TKE_B_Term = (nam.g./TH).*(mean(reshape(wpert.*thpert,nm,l))); % Bouyant production Term
TKE_S_Term = -(mean(reshape(upert.*wpert,nm,l))).*gradient(U, Z); % Shear Term
TKE_T_Term = -gradient(mean(reshape(wpert.*TKE_e,nm,l)), Z); % Turbulent transport of TKE
TKE_P_Term = -(1./RHO).*gradient(mean(reshape(wpert.*ppert,nm,l)), Z); % Pressure correlation term

%% Plot out the Fluxes that we need: 

% CTKE Budget 
figure('Renderer', 'painters', 'Position', [10 10 900 600]) % makes paper format figure

subplot(1,2,1)
plot(B_flx,Z/10^3,'LineWidth',2.0); hold on; grid on
ylabel('$Z$ [km]','FontSize',15, 'Interpreter','latex')
xlabel('Value','FontSize',15, 'Interpreter','latex')
legend(legendtext);
title("Bouyancy Flux, $\langle w'\theta_v' \rangle$", 'Interpreter','latex')
xlim([-0.1 0.1])
ylim([0 12])

subplot(1,2,2)
plot(M_flx,Z/10^3,'LineWidth',2.0); hold on; grid on
ylabel('$Z$ [km]','FontSize',15, 'Interpreter','latex')
xlabel('Value','FontSize',15, 'Interpreter','latex')
legend(legendtext);
title("Total Moisture Flux, $\langle w' q_t' \rangle$", 'Interpreter','latex')
xlim([-0.0005 0.0005])
ylim([0 12])

%% ========================================================================
%% non dimensionalize the equations by (rho* w*^3) /z* 
% Here the height will be scaled by z* -- defined as the height at which 
% the buoyancy flux near the  cloud  top  is  most  negative.

% plot out buoyancy flux to pinpoint location for z*
z_idx = 93; %84                                                                % index of height z* where bouyancy flux is most negative. 
z_star = Z(z_idx);

DLen = 5*z_star; 
D_Term = mean(reshape(TKE_e,nm,l)).^(-3/2)./ DLen;

TKE_Res = TKE_B_Term + TKE_S_Term + TKE_T_Term + TKE_P_Term + D_Term; 

figure()
plot(B_flx);

%% Compute rho_* term 

rho_star = (1/Z(z_idx)).*trapz(RHO(1:z_idx), Z(1:z_idx)); 

%% compute w^3_* term: 
z_trop = find(Z>=0 & Z<=12000); %isolate the tropospheric region 

% Compute mean tropospheric temperature
temp_trop = squeeze(mean(t(:,:,z_trop), [1,3])); % Average temperature over tropospheric altitudes
mean_temp_trop = mean(temp_trop); % Mean temperature of troposphere - THETA term in eq (1) of Khairoutdinov and Randall 2001

w_3_star = 2.5*(nam.g/(mean_temp_trop*rho_star)).*trapz(RHO(1:z_idx).*B_flx(1:z_idx), Z(1:z_idx)); 

scaling = (rho_star * w_3_star)/Z(z_idx);

%% Make all the plots 

% CTKE Budget 
figure('Renderer', 'painters', 'Position', [10 10 900 600]) % makes paper format figure

subplot(1,6,1)
plot(CTKE_T_Term/scaling,(Z/z_star),'LineWidth',2.0); hold on; grid on
ylabel('$z/z_\ast$','FontSize',15, 'Interpreter','latex')
xlabel('Dimensionless by $\rho_\ast w_\ast^3 / z_\ast$','FontSize',15, 'Interpreter','latex')
legend(legendtext);
title("Transport", 'Interpreter','latex')
xlim([-1 1])
ylim([0 1.5])

subplot(1,6,2)
plot(CTKE_P_Term/scaling,(Z/z_star),'LineWidth',2.0); hold on; grid on
ylabel('$z/z_\ast$','FontSize',15, 'Interpreter','latex')
xlabel('Dimensionless by $\rho_\ast w_\ast^3 / z_\ast$','FontSize',15, 'Interpreter','latex')
title("Pressure", 'Interpreter','latex')
xlim([-2.5 2.5])
ylim([0 1.5])
legend(legendtext);

subplot(1,6,3)
plot(CTKE_S_Term/scaling,(Z/z_star),'LineWidth',2.0); hold on; grid on
ylabel('$z/z_\ast$','FontSize',15, 'Interpreter','latex')
xlabel('Dimensionless by $\rho_\ast w_\ast^3 / z_\ast$','FontSize',15, 'Interpreter','latex')
title("Shear", 'Interpreter','latex')
ylim([0 1.5])
xlim([-0.02 0.02])
legend(legendtext);

subplot(1,6,4)
plot(CTKE_B_Term/scaling,(Z/z_star),'LineWidth',2.0); hold on; grid on
ylabel('$z/z_\ast$','FontSize',15, 'Interpreter','latex')
xlabel('Dimensionless by $\rho_\ast w_\ast^3 / z_\ast$','FontSize',15, 'Interpreter','latex')
title("Bouyancy", 'Interpreter','latex')
ylim([0 1.5])
xlim([-1.0 1.0])
legend(legendtext);

subplot(1,6,5)
plot(D_Term/scaling,(Z/z_star),'LineWidth',2.0); hold on; grid on
ylabel('$z/z_\ast$','FontSize',15, 'Interpreter','latex')
xlabel('Dimensionless by $\rho_\ast w_\ast^3 / z_\ast$','FontSize',15, 'Interpreter','latex')
title("Dissipation", 'Interpreter','latex')
ylim([0 1.5])
xlim([-1.0 1.0])
legend(legendtext);

subplot(1,6,6)
plot(mean(CTKE_Res)/scaling,(Z/z_star),'LineWidth',2.0); hold on; grid on
ylabel('$z/z_\ast$','FontSize',15, 'Interpreter','latex')
xlabel('Dimensionless by $\rho_\ast w_\ast^3 / z_\ast$','FontSize',15, 'Interpreter','latex')
title("Residual", 'Interpreter','latex')
ylim([0 1.5])
xlim([-1.0 1.0])
legend(legendtext);

sgtitle("$\bar{\mathcal{E}} \equiv 0.5 \bar{\rho}\left(\overline{u^{\prime 2}}+\overline{v^{\prime 2}}+\overline{w^{\prime 2}}\right)$", 'Interpreter', 'latex')
hold off; 


%% Next figure for the TKE Budgets 

% TKE Budget 
figure('Renderer', 'painters', 'Position', [10 10 900 600]) % makes paper format figure

subplot(1,6,1)
plot(TKE_T_Term,Z/10^3,'LineWidth',2.0); hold on; grid on
ylabel('$z$','FontSize',15, 'Interpreter','latex')
legend(legendtext);
title("Transport", 'Interpreter','latex')
xlim([-10^-3 10^-3])
ylim([0 12])

subplot(1,6,2)
plot(TKE_P_Term,Z/10^3,'LineWidth',2.0); hold on; grid on
ylabel('$z$','FontSize',15, 'Interpreter','latex')
title("Pressure", 'Interpreter','latex')
ylim([0 12])
xlim([-10^-3 10^-3])
legend(legendtext);

subplot(1,6,3)
plot(TKE_S_Term,Z/10^3,'LineWidth',2.0); hold on; grid on
ylabel('$z$','FontSize',15, 'Interpreter','latex')
title("Shear", 'Interpreter','latex')
ylim([0 12])
xlim([-10^-5 10^-5])
legend(legendtext);

subplot(1,6,4)
plot(TKE_B_Term,Z/10^3,'LineWidth',2.0); hold on; grid on
ylabel('$z$','FontSize',15, 'Interpreter','latex')
title("Bouyancy", 'Interpreter','latex')
ylim([0 12])
xlim([-10^-3 10^-3])
legend(legendtext);

subplot(1,6,5)
plot(D_Term,Z/10^3,'LineWidth',2.0); hold on; grid on
ylabel('$z$','FontSize',15, 'Interpreter','latex')
title("Dissipation", 'Interpreter','latex')
ylim([0 12])
xlim([-10^-3 10^-3])
legend(legendtext);

subplot(1,6,6)
plot(mean(TKE_Res),Z/10^3,'LineWidth',2.0); hold on; grid on
ylabel('$z$','FontSize',15, 'Interpreter','latex')
title("Residual", 'Interpreter','latex')
ylim([0 12])
xlim([-10^-3 10^-3])
legend(legendtext);

sgtitle("$\bar{\mathcal{E}} \equiv 0.5 \left(\overline{u^{\prime 2}}+\overline{v^{\prime 2}}+\overline{w^{\prime 2}}\right)$", 'Interpreter', 'latex')

