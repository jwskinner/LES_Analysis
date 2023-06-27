%% J.W.Skinner 2023 
% This script reads in data from specified location and computes the
% dissiapton using the subroutines in the WRF code. 

clear variables 

data_loc = '/data1/jwskinner/';
folder = 'GATE_NOEVP1.3km_CONSTFLX_50km/' %'GATE_NOCP_CONSTFLX/';

% index of file to start and finish the averaging over; set to same for
% instantanious histograms 
i_start = 90; 
i_end = 90; % 48 hours

% Takes the first folder for loading in params.
f_in = strcat(data_loc, folder); 

% Returns a list of all files in the folder
files_all = dir(strcat(f_in, 'wrfout*'));

% Load relevant physical constants into nam structure 
[nam] = get_constants(strcat(f_in, files_all(1).name));

% Load the vertical structure from a sample file
[Z, P, H, p, exn] = vert_struct(strcat(f_in, files_all(1).name), nam);

% Pre-allocate arrays for the time averaging
u_t = zeros((i_end - i_start), nam.nx-1, nam.ny, nam.levs); 
v_t = zeros((i_end - i_start), nam.nx-1, nam.ny, nam.levs);
w_t = zeros((i_end - i_start), nam.nx-1, nam.ny, nam.levs);
th_t = zeros((i_end - i_start), nam.nx-1, nam.ny, nam.levs);
tke_t = zeros((i_end - i_start), nam.nx-1, nam.ny, nam.levs);

tend = zeros((i_end - i_start), nam.nx-1, nam.ny, nam.levs);

j = 1;
for i = i_start: i_end

    file = files_all(i).name
    fname=strcat(f_in,file);

    [u, v, w, th, tke, qc, qr, qi, qv] =  ...
        loadNetCDF(fname, 'U', 'V', 'W', 'T', 'TKE', 'QCLOUD', 'QRAIN', 'QICE', 'QVAPOR');

    % Extract data component and convert to grid
    u = u.Data; v = v.Data; w = w.Data; th = th.Data + nam.T0; ...
        tke = tke.Data; qc = qc.Data; qr = qr.Data; qi = qi.Data; qv = qv.Data; 
    
    %C->A grid
    u=0.5*(u(1:end-1,:,:,:)+u(2:end,:,:,:)); 
    v=0.5*(v(:,1:end-1,:,:)+v(:,2:end,:,:)); 
    w=0.5*(w(:,:,1:end-1,:)+w(:,:,2:end,:));
    
    s=size(w); n=s(1); m=s(2); l=s(3); nm=n*m;
    
    TH=mean(reshape(th,nm,l));
    t=th.*exn;                                                             % temperature
    tv=t.*(1+0.608*qv);                                                    % virtual temperature, bouyancy is tv - ql (eq. 1 Marcin)
    thil=th-(nam.Ll*qc+nam.Li*qi)./(nam.cp*exn);                           % liquid water potential temperature
    rho=p./(nam.R*tv);                                                     % density

    % Load the Lat-Lon arrays from WRF 
    xlat =ncread(fname,'XLAT' );                                           % geopotential perturbation [m2/s2]
    xlon=ncread(fname,'XLONG');
    
    % Compute the corriolis parameter (zero in the LES code?)
    f = 2 * nam.omega *sin(xlat); 
    
%     % Calculate the diffusion terms
%     du_dt = -(u.*gradient(u)) + f.*v - 1/rho.*gradient(p) + K.*del2(u);
%     dv_dt = -(v.*gradient(v)) + u.*f - 1/rho.*gradient(q) + K.*del2(v);
%     dT_dt = -(u.*gradient(T)) - v.*w./rho + K.*del2(T);
    dx =nam.dx; dy = nam.dy; 
    
    % Check the size of the input array
    [rows, cols, levs] = size(u);

    % Create coordinate arrays for gradient calculation
    [x, y] = meshgrid(1:cols, 1:rows);
    
    % Preallocate arrays for gradient components
    du_dx = zeros(size(u));
    du_dy = zeros(size(u));
    dv_dx = zeros(size(v));
    dv_dy = zeros(size(v));

    %% Compute velocity gradients
    for k = 1:levs
        [du_dx(:,:,k),du_dy(:,:,k)]  = gradient(squeeze(u(:,:,k)), dx, dy);
%         du_dy(:,:,k) = gradient(squeeze(u(:,:,k)), dx, dy);
        [dv_dx(:,:,k), dv_dy(:,:,k)] = gradient(squeeze(v(:,:,k)), dx, dy);
%         dv_dy(:,:,k) = gradient(squeeze(v(:,:,k)), dx, dy);
    end

    % Parameters for Smagorinsky scheme
    Cs = 0.2; % Smagorinsky coefficient
    dx_squared = dx.^2;
    
    % Compute diffusion coefficient
    diff_coeff = Cs * dx_squared;
    
    % Compute diffusion terms
    du_dt = diff_coeff .* del2(du_dx) + diff_coeff .* del2(du_dy);
    dv_dt = diff_coeff .* del2(dv_dx) + diff_coeff .* del2(dv_dy);
    
    % Compute the mean of the diffusion terms along the horizontal dimensions
    du_dx_mean = squeeze(mean(du_dx, [1, 2]));
    du_dy_mean = squeeze(mean(du_dy, [1, 2]));
    dv_dx_mean = squeeze(mean(dv_dx, [1, 2]));
    dv_dy_mean = squeeze(mean(dv_dy, [1, 2]));
    
    du_dt_mean = squeeze(mean(du_dt, [1, 2]));
    dv_dt_mean = squeeze(mean(dv_dt, [1, 2]));
    
    %% Compute the TKE from github version of WRF code
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
    
    j = j + 1; 
end

%% plot vertical profiles of tend
figure(); 
plot(mean_tend, Z/1000, 'LineWidth', 1.5);hold on;
hold off;
title('TKE Dissipation');
xlabel('Tendency');
ylabel('Height (z) [km]');
xlim([-1e-3, 1e-3]);



%% Assuming du_dt and dv_dt are computed diffusion terms

% Create a grid of x and y coordinates
[x, y] = meshgrid(1:size(du_dt, 2), 1:size(du_dt, 1));

lev = 1;

% Plot du_dt
figure();
subplot(1, 2, 1)
contourf(x, y, mean(du_dt(:,:,:), 3),'LineStyle', 'none');
colorbar;
title('Vertical Mean Diffusion Term du/dt');
xlabel('x');
ylabel('y');

% Plot dv_dt
subplot(1, 2, 2)
contourf(x, y, mean(dv_dt(:,:,lev), 3),'LineStyle', 'none');
colorbar;
title('Vertical Mean Diffusion Term dv/dt');
xlabel('x');
ylabel('y');

%% plot vertical profiles of hoirzontally averaged diffusion terms 
figure(); 
subplot(1, 2, 1)
plot(du_dx_mean, Z/1000, 'LineWidth', 1.5);hold on;
plot(du_dy_mean, Z/1000, 'LineWidth', 1.5);hold on;
plot(dv_dx_mean, Z/1000, 'LineWidth', 1.5);hold on;
plot(dv_dy_mean, Z/1000, 'LineWidth', 1.5);hold on;
hold off;

legend('du/dx', 'du/dy', 'dv/dx', 'dv/dy');
title('Horizontally-Averaged Diffusion Terms');
xlabel('Diffusion Term');
ylabel('Height (z) [km]');

subplot(1, 2, 2)
plot(du_dt_mean, Z/1000, 'LineWidth', 1.5);hold on;
plot(dv_dt_mean, Z/1000, 'LineWidth', 1.5);hold on;
hold off;

legend('du/dt', 'dv/dt');
title('Horizontally-Averaged Diffusion Terms');
xlabel('Diffusion Term');
ylabel('Height (z) [km]');


% function [l_diss] = tke_dissip(tendency, mu, c1, c2, tke, dx, dy, rdzw, msftx, msfty, hpbl, dlk, l_diss)
% 
%     ce1 = (0.54 / 0.10) * 0.19;
%     ce2 = max(0.0, 0.93 - ce1);
% 
%     [ni, nj, nk] = size(tendency);
% 
%     l_scale = zeros(ni, nk, nj);
% 
%     for j = 1:nj
%         for k = 1:nk
%             for i = 1:ni
%                 deltas = (dx/msftx(i,j) * dy/msfty(i,j) / rdzw(i,k,j)) ^ 0.33333333;
%                 tketmp = max(tke(i,k,j), 1.0e-6);
% 
%                 coefc = ce1 + ce2 * l_scale(i,k,j) / deltas;
% 
%                 delxy = sqrt(dx/msftx(i,j) * dy/msfty(i,j));
%                 pu1 = pu(delxy,hpbl(i,j));
%                 coefc_m = 0.3;
%                 l_diss(i,k,j) = (1.0-pu1)*l_scale(i,k,j)/coefc + pu1*dlk(i,k,j)/coefc_m;
% 
%                 tendency(i,k,j) = tendency(i,k,j) - (c1(k)*mu(i,j)+c2(k)) * coefc * tketmp^1.5 / l_scale(i,k,j);
%             end
%         end
%     end
% end


% function [l_scale] = calc_l_scale(tke, BN2, dx, dy, rdzw, msftx, msfty, i_start, i_end, j_start, j_end)

% This subroutine calculates the length scale, based on stability,
% for TKE parameterization of subgrid-scale turbulence. It takes in:

% - tke: turbulent kinetic energy
% - BN2: turbulent buoyancy flux squared
% - l_scale: length scale
% - i_start, i_end: starting and ending indices for i dimension
% - ktf: highest index in k dimension
% - j_start, j_end: starting and ending indices for j dimension
% - dx, dy: grid spacings in x and y directions, respectively
% - rdzw: reciprocal of cell height
% - msftx, msfty: mask for solid-fluid cells
% - ids, ide: starting and ending indices for i dimension of the entire grid
% - jds, jde: starting and ending indices for j dimension of the entire grid
% - kds, kde: starting and ending indices for k dimension of the entire grid
% - ims, ime: starting and ending indices for i dimension of the subgrid
% - jms, jme: starting and ending indices for j dimension of the subgrid
% - kms, kme: starting and ending indices for k dimension of the subgrid
% - its, ite: starting and ending indices for i dimension of the inner subgrid
% - jts, jte: starting and ending indices for j dimension of the inner subgrid
% - kts, kte: starting and ending indices for k dimension of the inner subgrid
% - msftx, msfty: mask for solid-fluid cells in the subgrid

%     [ni, nk, nj] = size(tke);
%     l_scale = zeros(ni, nk, nj);
%     
%     for j = j_start:j_end
%         for k = 1:nk
%             for i = i_start:i_end
%                 deltas = (dx/msftx(i,j) * dy/msfty(i,j) / rdzw(i,k,j)) ^ 0.33333333;
%                 l_scale(i,k,j) = deltas;
%                 
%                 if BN2(i,k,j) > 1.0e-6
%                     tmp = sqrt(max(tke(i,k,j), 1.0e-6));
%                     l_scale(i,k,j) = 0.76 * tmp / sqrt(BN2(i,k,j));
%                     l_scale(i,k,j) = min(l_scale(i,k,j), deltas);
%                     l_scale(i,k,j) = max(l_scale(i,k,j), 0.001 * deltas);
%                 end
%             end
%         end
%     end
% end
