function [outputArg1,outputArg2] = plot_budget_terms(file, folder, nam, output)

fname=strcat(folder,file);

%% ========================================================================
% read all relevant 3-D variables

ph =ncread(fname,'PH' );                                                   % geopotential perturbation [m2/s2]
phb=ncread(fname,'PHB');                                                   % base geopotential [m2/s2)
p  =ncread(fname,'P'  );                                                   % pressure perturbation [Pa]
pb =ncread(fname,'PB' );                                                   % base pressure [Pa]
th =ncread(fname,'T'  )+nam.T0;                                            % Potential temperature [K]

qv=ncread(fname,'QVAPOR');                                                 % Water vapor mixing ratio [kg/kg]
qc=ncread(fname,'QCLOUD');                                                 % Cloud water mixing ratio [kg/kg]
%qr=ncread(fname,'QRAIN');                                                  % Rain water mixing ratio [kg/kg]
qi=ncread(fname,'QICE');                                                   % Ice mixing ratio [kg/kg]

%% 

%% ========================================================================
% initial calculations

s=size(th);
n=s(1); 
m=s(2); 
l=s(3); 
nm=n*m;

HS=mean(reshape(ph+phb,nm,l+1));                                           % Jack added +1 here to fix reshape [check] 
%H=0.5*(HS(1:end-1,:)+HS(2:end,:));                                        % Why is this an empty matrix (?)
H=0.5*(HS(:,1:end-1)+HS(:,2:end));                                         % Jack flipped the indexing to get the half values here 

Z=H./nam.g';                                                        % height at mass-levels [m]
p=p+pb;                                                                    % pressure
exn=(p/nam.P0).^(nam.R/nam.cp);                                            % exner function
qt=qv+qc+qi;                                                               % total water mixing ratio [kg/kg] (no precipitating elements)

t=th.*exn;                                                                 % temperature
tv=t.*(1+0.608*qv);                                                        % virtual temperature, bouyancy is tv - ql (eq. 1 Marcin)
rho=p./(nam.R*tv);                                                         % density

QTRHO = mean(reshape(qt.*rho,nm,l));

LWP = trapz(Z',QTRHO);                                                     % vertically integrated LWP

return LWP

end