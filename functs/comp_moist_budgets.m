function [qtvar, thilvar] = comp_moist_budgets(fname, nam, i)

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

rho=p./(nam.R*tv);              % density
RHO=mean(reshape(rho,nm,l));


%% ========================================================================
% QT variance budget terms -- Eq (1) Thomas et al. 2021


W=mean(reshape(w,nm,l));           
wpr(1,1,:)=W;
wpert=bsxfun(@minus,w,wpr);

QT=mean(reshape(qt,nm,l));                                                 % Horizontally average moisture variance 
qtpr(1,1,:)=QT;
qtpert=bsxfun(@minus,qt,qtpr);

% Caculate the gradient source terms in budget equation
wqt  = wpert.*qtpert;
wqt2 = wpert.*(qtpert.^2);

WQT =mean(reshape(wqt ,nm,l));
WQT2=mean(reshape(wqt2,nm,l));

QTVAR =mean(reshape(qtpert.^2 ,nm,l));

Qqt=qvtend+qctend+qitend;                                                  % Microphysical tendency term for sources

[dqtdx,dqtdy,dqtdz]=gradient((qtpert),nam.dx,nam.dy,50);

% What is this block for?
%zmat=shiftdim(repmat(diff(ZS),1,n,m),1);
%dqtdz=dqtdz.*50./zmat;          % Jack: Issue here with incompatible sizes

% all budget terms converted to units of [1/s*(g/kg)^2]
qtvar.prod=-2*WQT'.*gradient(QT',Z)*1e6;                                   % Gradient Production Term
%qtvar.trns=-1./RHO.*gradient(RHO.*WQT2',Z)*1e6;                           % old wasn't working because of the matrix transpose: '
qtvar.trns=(-1./RHO.*gradient(RHO.*WQT2,Z)*1e6)';                             % Turbulent transport term

qtvar.srcs=2*squeeze(mean(mean((qtpert).*Qqt)))*1e6;                       % Microphysical sources and sinks of the moisture variance
qtvar.diss=-2*squeeze(mean(mean(kh.*(dqtdx.^2+dqtdy.^2+dqtdz.^2))))*1e6;   % Eq. (2) of Schemann 2016


%% ========================================================================
% THIL variance budget terms -- Eq (5) Schemann 2017

THIL=mean(reshape(thil,nm,l));
thilpr(1,1,:)=THIL;
thilpert=bsxfun(@minus,thil,thilpr);

wthil=wpert.*thilpert;
wthil2=(wpert.^2).*(thilpert.^2);

WTHIL =mean(reshape(wthil ,nm,l));
WTHIL2 =mean(reshape(wthil2,nm,l));

WTHIL3 = mean(reshape((wpert).*(thilpert.^2),nm,l));                       % Transport term in Schemann 2017

Qthil=thtend-nam.Ll/nam.cp*qctend-nam.Li/nam.cp*qitend;

[dthildx,dthildy,dthildz]=gradient((thilpert),nam.dx,nam.dy,50);
%dthildz=dthildz.*50./zmat;

%thilvar.prod=-2*WTHIL'.*gradient(WTHIL',Z);   
thilvar.prod=-2*WTHIL'.*gradient(THIL',Z);                                 % Liquid water potential temperature, equation (5) Schemann 2017
%thilvar.trns=-1./RHO.*gradient(RHO.*WTHIL2,Z);                            
thilvar.trns=(-1./RHO.*gradient(RHO.*WTHIL3,Z))'; 
thilvar.srcs=2*squeeze(mean(mean((thilpert).*Qthil)));
thilvar.diss=-2*squeeze(mean(mean(kh.*(dthildx.^2+dthildy.^2+dthildz.^2))));


end