function [Z, U, V, W, TH, QT, QC, QR, TKE, TKE_HOR, TKE_W, TEMP, CLD_FRC, TOT_WAT] = oned_profiles(fname, nam)

%% ========================================================================
% read all relevant 3-D variables

u=ncread(fname,'U'); u=0.5*(u(1:end-1,:,:,:)+u(2:end,:,:,:)); %C->A grid
v=ncread(fname,'V'); v=0.5*(v(:,1:end-1,:,:)+v(:,2:end,:,:));
w=ncread(fname,'W'); w=0.5*(w(:,:,1:end-1,:)+w(:,:,2:end,:));

% Can put this is in to import the variables because it retains the units
% and description for each one which is useufl later.
% [ph, phb, tke, p, pb, th] = loadNetCDF(fname, 'PH', 'PHB', 'TKE', 'P', ...
%     'PB', 'T');
%
% [qv, qc, qr, qi, qs, qg, kh, qvtend, qctend, qitend, thtend] = ...
%     loadNetCDF(fname, 'QVAPOR','QCLOUD','QRAIN', 'QICE', 'QSNOW',...
%     'QGRAUP', 'XKHH','QVTEND_MICRO','QCTEND_MICRO','QITEND_MICRO', ...
%     'THTEND_MICRO');


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

TKE_hor = 0.5*(upert.^2 + vpert.^2);
TKE_HOR=mean(reshape(TKE_hor,nm,l));
TKE_W = 0.5*mean(reshape(wpert,nm,l)).^2;

%WAT_FLUX = (mean(reshape(wpert,nm,l)).*mean(reshape(qtpert,nm,l)));       % Total water flu

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

% Calculation of total water from air mass and total water mixing
% ratio.
% q_t = mass_w / mass_air  (density water vapour / air density)
% q_t = mass_w / (rho*V)
% mass_w = q_t *(rho*V)

v_air = (nam.dx * nam.dy * (max(Z)/nam.levs));                             % Volume of each grid box [m^3]
M_AIR = RHO .* v_air;
TOT_WAT = (M_AIR .* QT);                                                   % Total water [Kg]


