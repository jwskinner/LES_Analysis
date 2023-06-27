%% J.W.Skinner 2023
% This script reads in data from specified location and time averages
% the 1D profiles and plots the profiles and fluxes for the paper

clear variables

data_loc = '/data1/jwskinner/';
folder_cp = 'GATE_CP_CONSTFLX/';
folder_nocp = 'GATE_NOEVP1.3km_CONSTFLX_100km/';

% index of file to start and finish the averaging over
i_start = 89;
i_end = 97; % 48 hours

% Takes the first folder for loading in params.
f_in_cp = strcat(data_loc, folder_cp);
f_in_nocp = strcat(data_loc, folder_nocp);

% Returns a list of all files in the folder
files_all = dir(strcat(f_in_cp, 'wrfout*'));

% Load relevant physical constants into nam structure
[nam] = get_constants(strcat(f_in_cp, files_all(1).name));

% Preallocate the time arrays for the CP and NOCP cases
data_cp_t = zeros((i_start - i_end), 11, nam.levs);
data_nocp_t = zeros((i_start - i_end), 11, nam.levs);

% Load the vertical structure
[Z, P, H, p, exn] = vert_struct(strcat(f_in_cp, files_all(1).name), nam);

j = 1;
for i = i_start: i_end
    
    file = files_all(i).name
    fname_cp=strcat(f_in_cp,file);
    fname_nocp=strcat(f_in_nocp,file);
    
    data_cp_t(i,:,:) = get_1d_profs(fname_cp,nam, (i_end - i_start), exn);
    data_nocp_t(i,:,:) = get_1d_profs(fname_nocp, nam,(i_end - i_start), exn);
    
    j = j + 1;
end

%% Average over the time dimension of the arrays
data_cp_tmean = mean(data_cp_t, 1);
data_nocp_tmean = mean(data_nocp_t, 1);

% Plot the budgets
figure('Renderer', 'painters', 'Position', [10 10 900 600]) % makes paper format figure

lw = 1.5;
colors = {'#0072BD', '#A2142F', '#EDB120', '#77AC30', '#000000'};
style = {'-', '-', '-', '-', '--'};

subplot(1, 5, 1)
plot(squeeze(data_cp_tmean(1, 4, :)),Z/10^3,'k','LineWidth',lw); hold on; 
plot(squeeze(data_nocp_tmean(1, 4, :)),Z/10^3,'b','LineWidth',lw); hold on; 
ylabel('$z$ [km]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
xlabel('[ms$^{-1}$]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
legend('CP', 'NOCP', 'Box', 'off');
title('$\langle {\rm U} \rangle$', 'FontSize',12, 'Interpreter', 'Latex'); hold off;
xlim([-0.15 0])
ylim([0 18])

subplot(1, 5, 2)
plot(squeeze(data_cp_tmean(1, 5, :)),Z/10^3,'k', 'LineWidth',lw); hold on; 
plot(squeeze(data_nocp_tmean(1, 5, :)),Z/10^3,'b','LineWidth',lw); hold on;
ylabel('$z$ [km]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
xlabel('[ms$^{-1}$]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
legend('CP', 'NOCP', 'Box', 'off');
title('$\langle {\rm V} \rangle$', 'FontSize',12, 'Interpreter', 'Latex'); hold off;
xlim([-1e-2 1e-2])
ylim([0 18])

subplot(1, 5, 3)
plot(squeeze(data_cp_tmean(1, 6, :)),Z/10^3,'k','LineWidth',lw); hold on; 
plot(squeeze(data_nocp_tmean(1, 6, :)),Z/10^3,'b','LineWidth',lw); hold on;
ylabel('$z$ [km]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
xlabel('[ms$^{-1}$]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
legend('CP', 'NOCP', 'Box', 'off');
title('$\langle {\rm W} \rangle$', 'FontSize',12, 'Interpreter', 'Latex'); hold off;
xlim([-1e-4 1e-4])
ylim([0 18])

subplot(1, 5, 4)
plot(squeeze(data_cp_tmean(1, 7, :)),Z/10^3,'k','LineWidth',lw); hold on; 
plot(squeeze(data_nocp_tmean(1, 7, :)),Z/10^3,'b','LineWidth',lw); hold on;
ylabel('$z$ [km]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
xlabel('[ms$^{-1}$]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
legend('CP', 'NOCP', 'Box', 'off');
title('$\langle {\rm u}^\prime {\rm w }^\prime \rangle$', 'FontSize',12, 'Interpreter', 'Latex'); hold off;
% xlim([-1e-4 1e-4])
ylim([0 18])

subplot(1, 5, 5)
plot(squeeze(data_cp_tmean(1, 8, :)),Z/10^3,'k','LineWidth',lw); hold on; 
plot(squeeze(data_nocp_tmean(1, 8, :)),Z/10^3,'b','LineWidth',lw); hold on;
ylabel('$z$ [km]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
xlabel('[ms$^{-1}$]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
legend('CP', 'NOCP', 'Box', 'off');

title('$\langle {\rm v}^\prime {\rm w }^\prime \rangle$', 'FontSize',12, 'Interpreter', 'Latex'); hold off;
% xlim([-1e-4 1e-4])
ylim([0 18])

%% Plot the higher moments of vertical velocity and the fluxes 
figure('Renderer', 'painters', 'Position', [10 10 900 600]) % makes paper format figure

lw = 1.5;
colors = {'#0072BD', '#A2142F', '#EDB120', '#77AC30', '#000000'};
style = {'-', '-', '-', '-', '--'};

subplot(1, 3, 1)
plot(squeeze(data_cp_tmean(1, 11, :)),Z/10^3,'k','LineWidth',lw); hold on; grid on; 
plot(squeeze(data_nocp_tmean(1, 11, :)),Z/10^3,'b','LineWidth',lw); hold on; 
ylabel('$z$ [km]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
xlabel('[ms$^{-1}$]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
legend('CP', 'NOCP', 'Box', 'off');
title('$\langle {\rm w}^\prime \theta^\prime \rangle$', 'FontSize',12, 'Interpreter', 'Latex'); hold off;
% xlim([-0.15 0])
ylim([0 18])

subplot(1, 3, 2)
plot(squeeze(data_cp_tmean(1, 9, :)),Z/10^3,'k','LineWidth',lw); hold on;grid on; 
plot(squeeze(data_nocp_tmean(1, 9, :)),Z/10^3,'b','LineWidth',lw); hold on; 
ylabel('$z$ [km]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
xlabel('[m$^2$s$^{-2}$]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
legend('CP', 'NOCP', 'Box', 'off');
title('$\langle {\rm w}^{\prime \, 2} \rangle$', 'FontSize',12, 'Interpreter', 'Latex'); hold off;
xlim([-0.01 0.04])
ylim([0 18])

subplot(1, 3, 3)
plot(squeeze(data_cp_tmean(1, 10, :)),Z/10^3,'k','LineWidth',lw); hold on;grid on; 
plot(squeeze(data_nocp_tmean(1, 10, :)),Z/10^3,'b','LineWidth',lw); hold on; 
ylabel('$z$ [km]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
xlabel('[m$^3$s$^{-3}$]','LineWidth',1.5,'FontSize',15, 'Interpreter', 'Latex')
legend('CP', 'NOCP', 'Box', 'off');
title('$\langle {\rm w}^{\prime \, 3} \rangle$', 'FontSize',12, 'Interpreter', 'Latex'); hold off;
 xlim([-0.1 0.1])
ylim([0 18])


% Customaisable function for pulling in all the data for the 1D profiles in
% time
function data_out_t = get_1d_profs(fname,nam, tlen, exn)
    
    nlev = nam.levs;
    data_out_t = zeros(11, nlev);
    
    % Sets up the paramters for the plotting
    [u, v, w, th, qv] = loadNetCDF(fname, 'U', 'V', 'W', 'T', 'QVAPOR');
    u = u.Data; v = v.Data; w = w.Data; th = th.Data+nam.T0; qv = qv.Data;
    
    u=0.5*(u(1:end-1,:,:,:)+u(2:end,:,:,:)); %C->A grid
    v=0.5*(v(:,1:end-1,:,:)+v(:,2:end,:,:));
    w=0.5*(w(:,:,1:end-1,:)+w(:,:,2:end,:));
    
    s=size(w);n=s(1);m=s(2);l=s(3);nm=n*m;
    
    % Compute the purturbations for the fluxes
    W=mean(reshape(w,nm,l));
    wpr(1,1,:)=W;
    wpert=bsxfun(@minus,w,wpr);
    
    U=mean(reshape(u,nm,l));
    upr(1,1,:)=U;
    upert=bsxfun(@minus,u,upr);
    
    V=mean(reshape(v,nm,l));
    vpr(1,1,:)=V;
    vpert=bsxfun(@minus,v,vpr);
    
    TH=mean(reshape(th,nm,l));
    thpr(1,1,:)=TH;
    thpert=bsxfun(@minus,th,thpr);
    
    t=th.*exn; % temperature 
    tv=t.*(1+0.608*qv);                                                    % virtual temperature, bouyancy is tv - ql (eq. 1 Marcin)
    % thil=th-(nam.Ll*qc+nam.Li*qi)./(nam.cp*exn);                         % liquid water potential temperature
    % Virtual potential temperature 
    TV=mean(reshape(tv,nm,l));
    tvpr(1,1,:)=TV;
    tvpert=bsxfun(@minus,tv,tvpr);
    
    %output the data in the time loop 
    data_out_t(1,:) = mean(reshape(wpert, nm, l));
    data_out_t(2,:) = mean(reshape(vpert, nm, l));
    data_out_t(3,:) = mean(reshape(wpert, nm, l));
    data_out_t(4,:) = U;
    data_out_t(5,:) = V;
    data_out_t(6,:) = W;
    data_out_t(7,:) = mean(reshape(upert.*wpert, nm, l));
    data_out_t(8,:) = mean(reshape(vpert.*wpert, nm, l));
    data_out_t(9,:) = mean(reshape(wpert.*wpert, nm, l));
    data_out_t(10,:) = mean(reshape(wpert.*wpert.*wpert, nm, l));
    data_out_t(11,:) = mean(reshape(wpert.*tvpert, nm, l)); % bouyancy flux
end