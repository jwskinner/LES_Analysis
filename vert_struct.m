function [Z, P, H, p, exn] = vert_struct(fname, nam)

    ph =ncread(fname,'PH' );                                               % geopotential perturbation [m2/s2]
    phb=ncread(fname,'PHB');                                               % base geopotential [m2/s2)
    p  =ncread(fname,'P'  );                                               % pressure perturbation [Pa]
    pb =ncread(fname,'PB' );                                               % base pressure [Pa]
    th =ncread(fname,'T'  )+nam.T0;                                        % Potential temperature [K]

    s=size(phb); 
    n=s(1); m=s(2); l=s(3); nm=n*m;

    HS=mean(reshape(ph+phb,nm,l));                                             
    H=0.5*(HS(:,1:end-1)+HS(:,2:end));                                         

    ZS=HS./nam.g;                                                          % height at w-levels [m]
    Z=H./nam.g';                                                           % height at mass-levels [m]
    p=p+pb;                                                                % pressure
  
    P = mean(reshape(p,nm,l-1)); 
    exn=(p/nam.P0).^(nam.R/nam.cp);                                        % exner function

