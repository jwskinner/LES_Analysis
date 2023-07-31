function [cke_var] = comp_cke_budgets(fname, nam, i)

time = nam.dt * i;

%% ========================================================================
% read all relevant 3-D variables

[u, v, w] = loadNetCDF(fname, 'U', 'V', 'W'); 

u = u.Data; v = v.Data; w = w.Data; 

u=0.5*(u(1:end-1,:,:,:)+u(2:end,:,:,:)); %C->A grid
v=0.5*(v(:,1:end-1,:,:)+v(:,2:end,:,:)); 
w=0.5*(w(:,:,1:end-1,:)+w(:,:,2:end,:));

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

TH=mean(reshape(th,nm,l));
t=th.*exn;                                                                 % temperature
tv=t.*(1+0.608*qv);                                                        % virtual temperature, bouyancy is tv - ql (eq. 1 Marcin)
thil=th-(nam.Ll*qc+nam.Li*qi)./(nam.cp*exn);                               % liquid water potential temperature

rho=p./(nam.R*tv);                                                         % density
RHO=mean(reshape(rho,nm,l));

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
%TKE_test = TKE - mean(reshape(0.5*(u.^2+v.^2),nm,l));

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

% KE terms 
ke = 0.5*rho.*(upert.^2 + vpert.^2);
KE=mean(reshape(ke,nm,l));
kepr(1,1,:)=KE;
ke_pert=bsxfun(@minus,ke,kepr);

%% ========================================================================
% CKE Budget terms from Equation (5) of Khairoutdinov and Randall 2001

% Transport term 
WE_pert = mean(reshape(wpert .* cke_pert,nm,l));
dWE_dz = gradient(WE_pert, Z);
CTKE_T_Term = -dWE_dz; % Jack checked ok.

% Pressure Correlation term 
WP = mean(reshape(wpert.*ppert,nm,l));
dWP_dz = gradient(WP, Z);
CTKE_P_Term = - dWP_dz; % Jack the sign of plot isn't consistent with K&R fig 6 -- chase this down.

% Mean vertical wind shear term
UW_pert = mean(reshape(upert .* wpert,nm,l));
du_dz = gradient(U, Z);
VW_pert = mean(reshape(vpert .* wpert,nm,l));
dv_dz = gradient(V, Z);
CTKE_S_Term = -RHO .* (UW_pert.*du_dz + VW_pert.*dv_dz); % Jack checked ok.


% Bouyancy production term 
Beta = nam.g ./ TV; % From Stull page 152 

B_flx = mean(reshape(wpert.*tvpert,nm,l));  % buoyancy flux
CTKE_B_Term = RHO.*Beta.*B_flx;

% Viscous dissipation term [old method]
% [dthildx,dthildy,dthildz]=gradient((thilpert),nam.dx,nam.dy,50);
% D_Term=-2*squeeze(mean(mean(kh.*(dthildx.^2+dthildy.^2+dthildz.^2))));

%% Compute the TKE dissiaption from the WRF model code [new method]
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
    
     % compute the mu term from da_transform_xtoxa.inc
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

cke_var.Bterm = CTKE_B_Term; 
cke_var.Sterm = CTKE_S_Term;
cke_var.Tterm = CTKE_T_Term;
cke_var.Pterm = CTKE_P_Term;
cke_var.Dterm = D_Term;
cke_var.Res = CTKE_Res; 


end