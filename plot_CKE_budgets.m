% J.W.Skinner 07/27/23 
% This is a test script used to develop the functions for plotting the CKE budgets from kardunov and
% randall paper 

addpath('../cmocean-main/')
addpath('./functs/')

data_loc = '/data1/jwskinner/';
 folder = 'GATE_NOEVP1.3km_CONSTFLX_50km/';
% folder = 'GATE_NOCP_CONSTFLX/';

% Takes the first folder for loading in params.
f_in = strcat(data_loc, folder); 

% Returns a list of all files in the folder
files_all = dir(strcat(f_in, 'wrfout*'));

i = 96; 
file = files_all(i).name
fname=strcat(f_in,file);

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
Beta = nam.g ./ TV; % JACK double check this term 

B_flx = mean(reshape(wpert.*tvpert,nm,l));  % buoyancy flux
CTKE_B_Term = -RHO.*Beta.*B_flx;

% Viscous dissipation term [[OLD WAY]]
[dthildx,dthildy,dthildz]=gradient((thilpert),nam.dx,nam.dy,50);
% D_Term=-2*squeeze(mean(mean(kh.*(dthildx.^2+dthildy.^2+dthildz.^2))));

%% Compute the TKE dissiaption from the WRF model code [[NEW WAY]]
dx =nam.dx; dy = nam.dy; 

[ni, nj, nk] = size(tke); % get sizes for the array loops
    
    delxy = sqrt(dx * dy);
    
    tketmp = max(tke);
    ce1 = (0.54 / 0.10) * 0.19;
    ce2 = max(0.0, 0.93 - ce1);
    
    % Compute the Length Scale and coeficients
    l_scale = zeros(ni, nj, nk);
    coefc = zeros(ni, nj, nk); 
    c1 = zeros(nk); c2 = zeros(nk);  
    
    dz8w = gradient(Z); %dz between full levels (m)
    rdzw = 1./dz8w; % Spacing between full levels in the model
   
    deltas = (dx * dy ./ rdzw) .^ 0.33333333; % calculate a length scale proportional to grid spacing
    for k = 1:nk
        c1(k) = 1.0; c2(k) = 1.0; % Setting these to 1 following "subroutine pbl_diffusion_em(k, p_hl, dudt, dvdt)"
        l_scale(:,:,k) = deltas(k); 
        coefc(:,:, k) = ce1 + ce2 * l_scale(:,:,k) / deltas(k);
    end
    
     % compute the mu term from da_transform_xtoxa.inc [We need this!]
     %mu(i, j)= 1.0; %-(grid%xa%psfc(i,j)+grid%xb%psac(i,j)*sdmd)/s1md
                
     % Rough estimation of mu for now
     kappa = 0.4;    % von Karman constant
     epsilon = 0.1;  % Turbulent dissipation rate
     mu = kappa .* (tke.^2) ./ epsilon;

     % Compute the tendency
     tend = zeros(ni, nj, nk);
     mean_tend = zeros(nk); % vertical mean of tendency  
    
    for i2 = 1:ni
        for j2 = 1:nj
            for k2 = 1:nk
                tketmp = max(tke(i2,j2,k2), 1.0e-6);
                
                tend(i2,j2,k2) = tend(i2,j2,k2) - (c1(k2)*mu(i2,j2)+c2(k2)) * ...
                    coefc(i2,j2,k2) * tketmp^1.5 / l_scale(i2,j2,k2);
            end
        end
    end
    
    [n, m, l] = size(tend); 
    mean_tend = squeeze(mean(reshape(tend, n*m, l))); 
    
    D_Term = mean_tend;


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

%% Plot out the Budgets that we need: 

figure(); 
plot(mean_tend, Z/1000, 'LineWidth', 1.5);hold on;
plot(D_Term, Z/1000, 'LineWidth', 1.5);hold on;
hold off;
title('TKE Dissipation');
xlabel('Tendency');
ylabel('Height (z) [km]');
xlim([-1e-3, 1e-3]);

% CTKE Budget 
figure('Renderer', 'painters', 'Position', [10 10 900 600]) % makes paper format figure

subplot(1,2,1)
plot(B_flx,Z/10^3,'k','LineWidth',2.0); hold on; 
ylabel('$z$ [km]','FontSize',15, 'Interpreter','latex')
xlabel('Value','FontSize',15, 'Interpreter','latex')
% legend(legendtext);
title("Bouyancy Flux, $\langle w'\theta_v' \rangle$", 'Interpreter','latex')
xlim([-0.05 0.05])
ylim([0 18])

subplot(1,2,2)
plot(M_flx,Z/10^3,'k','LineWidth',2.0); hold on; 
ylabel('$z$ [km]','FontSize',15, 'Interpreter','latex')
xlabel('Value','FontSize',15, 'Interpreter','latex')
% legend(legendtext);
title("Total Moisture Flux, $\langle w' q_t' \rangle$", 'Interpreter','latex')
xlim([-0.00025 0.00025])
ylim([0 18])

%% ========================================================================
%% non dimensionalize the equations by (rho* w*^3) /z* 
% Here the height will be scaled by z* -- defined as the height at which 
% the buoyancy flux near the  cloud  top  is  most  negative.

% plot out buoyancy flux to pinpoint location for z*
z_idx = 17; %84                                                            % index of height z* where bouyancy flux is most negative. 
z_star = Z(z_idx);

DLen = 5*z_star; 
% D_Term = mean(reshape(TKE_e,nm,l)).^(-3/2)./ DLen; % Jack why??

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
xlim([-2 2])
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
xlim([-2.0 2.0])
legend(legendtext);

subplot(1,6,5)
plot(D_Term/scaling,(Z/z_star),'LineWidth',2.0); hold on; grid on
ylabel('$z/z_\ast$','FontSize',15, 'Interpreter','latex')
xlabel('Dimensionless by $\rho_\ast w_\ast^3 / z_\ast$','FontSize',15, 'Interpreter','latex')
title("Dissipation", 'Interpreter','latex')
ylim([0 1.5])
xlim([-2.0 2.0])
legend(legendtext);

subplot(1,6,6)
plot(CTKE_Res/scaling,(Z/z_star),'LineWidth',2.0); hold on; grid on
ylabel('$z/z_\ast$','FontSize',15, 'Interpreter','latex')
xlabel('Dimensionless by $\rho_\ast w_\ast^3 / z_\ast$','FontSize',15, 'Interpreter','latex')
title("Residual", 'Interpreter','latex')
ylim([0 1.5])
xlim([-2.0 2.0])
legend(legendtext);

sgtitle("$\bar{\mathcal{E}} \equiv 0.5 \bar{\rho}\left(\overline{u^{\prime 2}}+\overline{v^{\prime 2}}+\overline{w^{\prime 2}}\right)$", 'Interpreter', 'latex')
hold off; 


%% Next figure for the TKE Budgets 

% TKE Budget 
figure('Renderer', 'painters', 'Position', [10 10 900 600]) % makes paper format figure
plot(TKE_T_Term,Z/10^3,'LineWidth',2.0); hold on; grid on
plot(TKE_P_Term,Z/10^3,'LineWidth',2.0); hold on; 
plot(TKE_S_Term,Z/10^3,'LineWidth',2.0); hold on;
plot(TKE_B_Term,Z/10^3,'LineWidth',2.0); hold on;
plot(D_Term,Z/10^3,'LineWidth',2.0); hold on;
plot(TKE_Res,Z/10^3,'LineWidth',2.0); hold on;
ylabel('$z$','FontSize',15, 'Interpreter','latex')
legend("TRANS.", "PROD.", "SHEAR", "BOUY.", "DISS.", "RES.");
sgtitle("$\bar{\mathcal{E}} \equiv 0.5 \left(\overline{u^{\prime 2}}+\overline{v^{\prime 2}}+\overline{w^{\prime 2}}\right)$", 'Interpreter', 'latex')
xlim([-10^-3 10^-3])
ylim([0 18])

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
plot(TKE_Res,Z/10^3,'LineWidth',2.0); hold on; grid on
ylabel('$z$','FontSize',15, 'Interpreter','latex')
title("Residual", 'Interpreter','latex')
ylim([0 12])
xlim([-10^-3 10^-3])
legend(legendtext);

sgtitle("$\bar{\mathcal{E}} \equiv 0.5 \left(\overline{u^{\prime 2}}+\overline{v^{\prime 2}}+\overline{w^{\prime 2}}\right)$", 'Interpreter', 'latex')

