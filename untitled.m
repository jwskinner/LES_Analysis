% Load wind field data for the u component at one height z
fname = "./data/large_domain/NOCP_OUT/wrfout_d01_2020-01-01_22:00:00";
[u_nc, tke_nc] = loadNetCDF(fname, 'U', 'TKE')

u = u_nc.Data; % take u field at some height z

% Define the size of the velocity field
[M, N, L] = size(u);
TKE=mean(reshape(tke_nc.Data,N*N,L)); % mean TKE profile

% Find the height at which TKE reduces to 5% of its near-surface value for
% the SBL height
dh = 1; % vertical grid spacing
h = find(TKE < 0.05*TKE(1), 1, 'first')*dh;

% Define the spatial separation vector (rx, ry, rz)
rx = 200;
ry = 200;

% Preallocate arrays
dx = 10;                   % x spacing
dy = 10;                   % y spacing
Rij = zeros(N,M);           % initialize Rij array
num = zeros(N,M);           % initialize Rij array
den = zeros(N,M);           % initialize Rij array
x = ((0:N-1)*dx);           % x coordinate array
y = ((0:M-1)*dy);           % y coordinate array
z = round(0.2*h);                    % height to perform the analysis at

% Compute the two-point correlation function R11
[R] = two_point_correlation(u, z, rx, ry);
[R1X] = oned_correlation(u, z, 1000, 0, rx);
[R1Y] = oned_correlation(u, z, 1000, 1, ry);

% Compute the autocorrelation function with matlab
u_pert = u(:,:,z) - mean(u(:)); 
corr_x = xcorr(u_pert(1000,:));
corr_y = xcorr(u_pert(:,1000));

%% Plot the result
figure()
subplot(3,2,1)
imagesc(R);
colormap jet;
colorbar;
caxis([-1.0 1.0])
title('Correlation Function');
xlabel('x');
ylabel('y');

subplot(3,2,2)
imagesc(u(:,:,z));
colormap jet;
colorbar;
xlabel('x');
ylabel('y');
title("u field");

subplot(3,2,3)
plot(R1X);
title("Cross Correlation, x");

subplot(3,2,4)
plot(R1Y);
title("Cross Correlation, y");

subplot(3,2,5)
plot(corr_x);
title("Matlab Correlation, x");
xlabel('Lag');
ylabel('Correlation');


subplot(3,2,6)
plot(corr_y);
title("Matlab Correlation, y");
xlabel('Lag');
ylabel('Correlation');

%% All these useful functions 
function R = two_point_correlation(u, z, rx, ry)
% Compute two-point correlation function R
mean_u = mean(u(:));  % Compute mean over entire velocity field
[M, N, P] = size(u);
R = zeros(N, M);

for i = 1:M
    for j = 1:N
        xi = i; yi = j;
        % check the roving location is within index bounds
        if xi+rx > 0 && xi+rx <= M && yi+ry > 0 && yi+ry <= N

            % Compute the local velocity fluctuations at the reference
            % location ui'(x, y, z) and roving location uj'(x+rx, y+rx, z+rz)
            u_ref = u(i,j,z) - mean_u;
            u_rov = u(xi+rx,yi+ry,z) - mean_u;

            % Compute the turbulence intensity, <ui'^2(x, y, z)> at the reference location
            hu2_ref = mean(u_ref.^2);

            % Compute the turbulence intensity term hu2 at the roving
            % location <uj'^2(x+rx, y+ry, z)>
            hu2_rov = mean(u_rov .^ 2);

            % Compute the numerator of the correlation function Rij, which is the
            % spatial average of the product of the velocity fluctuations,
            % <ui^2'(x,y,z)*uj^2'(x+rx, y+ry, z)>
            num = mean(u_ref .* u_rov);

            % Compute the denominator of the correlation function Rij, which is the
            % product of the turbulence intensity terms at the reference
            % locations (<ui'^2(z)><uj'^2(z)>)^0.5
            den = sqrt(hu2_ref * hu2_ref);

            % Compute the correlation function Rij
            R(j,i) = num / den;
        end
    end
end
end

% Compute two-point correlation function R; u = data; z = height; index =
% position; xy = direction (either x or y direction); r = box size; 
%
% This is a function in MATLAB for computing the one-dimensional two-point 
% correlation function, given input data, height, index, direction, and 
% box size. The correlation function is computed using local velocity 
% fluctuations and turbulence intensity terms at the reference and roving 
% locations, and returned as output. The code first computes the mean over
% the entire velocity field and initializes the output matrix R with zeros.
% The correlation function is then computed separately for the 
% x and y directions.

function R = oned_correlation(u, z, index, xy, r)

mean_u = mean(u(:));  % Compute mean over entire velocity field
[M, N, P] = size(u);
R = zeros(N);

if xy ==0 
    rx = r; 
    for i = 1:M
        xi = i;
        % check the roving location is within index bounds
        if xi+rx > 0 && xi+rx <= M
            % Compute the local velocity fluctuations at the reference
            % location ui'(x, y, z) and roving location uj'(x+rx, y+rx, z+rz)
            u_ref = u(i, index,z) - mean_u;
            u_rov = u(xi+rx,index,z) - mean_u;

            % Compute the turbulence intensity, <ui'^2(x, y, z)> at the reference location
            hu2_ref = mean(u_ref.^2);

            % Compute the turbulence intensity term hu2 at the roving
            % location <uj'^2(x+rx, y+ry, z)>
            hu2_rov = mean(u_rov .^ 2);

            % Compute the numerator of the correlation function Rij, which is the
            % spatial average of the product of the velocity fluctuations,
            % <ui^2'(x,y,z)*uj^2'(x+rx, y+ry, z)>
            num = mean(u_ref .* u_rov);

            % Compute the denominator of the correlation function Rij, which is the
            % product of the turbulence intensity terms at the reference
            % locations (<ui'^2(z)><uj'^2(z)>)^0.5
            den = sqrt(hu2_ref * hu2_ref);

            % Compute the correlation function Rij
            R(i) = num / den;
        end
    end
else
    ry = r; 
    for i = 1:N
        yi = i;
        % check the roving location is within index bounds
        if yi+ry > 0 && yi+ry <= N
            % Compute the local velocity fluctuations at the reference
            % location ui'(x, y, z) and roving location uj'(x+rx, y+rx, z+rz)
            u_ref = u(index, i,z) - mean_u;
            u_rov = u(index,yi+ry,z) - mean_u;

            % Compute the turbulence intensity, <ui'^2(x, y, z)> at the reference location
            hu2_ref = mean(u_ref.^2);

            % Compute the turbulence intensity term hu2 at the roving
            % location <uj'^2(x+rx, y+ry, z)>
            hu2_rov = mean(u_rov .^ 2);

            % Compute the numerator of the correlation function Rij, which is the
            % spatial average of the product of the velocity fluctuations,
            % <ui^2'(x,y,z)*uj^2'(x+rx, y+ry, z)>
            num = mean(u_ref .* u_rov);

            % Compute the denominator of the correlation function Rij, which is the
            % product of the turbulence intensity terms at the reference
            % locations (<ui'^2(z)><uj'^2(z)>)^0.5
            den = sqrt(hu2_ref * hu2_ref);

            % Compute the correlation function Rij
            R(i) = num / den;
        end
    end
end
end

%% Jack old code clean this up! :) 


% % Pre-allocate the output variable
% Rij = zeros(M, N);

% % Compute the cross-correlation between the fluctuating velocity field at two points
% Rij = mean(mean(hu .* circshift(hu, [rx, ry]))) ./ (sqrt(hu2) .* sqrt(hu2));
%
% % Correct for the spatial mean
% Rij = Rij - mean(mean(Rij));


% % Define the spatial separation between points
% dx = 20;
% dy = 20;
%
% % Initialize the two-point correlation matrix
% R = zeros(m, n);
%
% % Calculate the two-point correlation function
% for i = 1:m
%     for j = 1:n
%         % Nested loop to calculate correlation between each pair of points
%         for ii = 1:m
%             for jj = 1:n
%                 % Calculate the relative spatial position of the two points
%                 rx = (i - ii) * dx;
%                 ry = (j - jj) * dy;
%
%                 % Check if the relative position is within the bounds of the matrix
%                 if (i + rx > 0) && (i + rx <= m) && (j + ry > 0) && (j + ry <= n)
%                     % Add the correlation value between the two points to the sum
%                     R(i, j) = R(i, j) + perturbation(ii, jj) * perturbation(i + rx, j + ry);
%                 end
%             end
%         end
%         % Normalize the correlation sum for each pair of points to get the
%         % spatial average
%         R(i, j) = R(i, j) / (m * n);
%     end
% end
%

%% Loop over x and y values
%
% for i = 1:M
%     for j = 1:N
%         xi = i; yi = j;
%         % check the roving location is within index bounds
%         if xi+rx > 0 && xi+rx <= M && yi+ry > 0 && yi+ry <= N
%
%             % Define the general perturbation from the velocity field
%             % u'= u(x, y, z) - mean(u(x, y, z))
%             %mean_u = mean(reshape(u,N*M,L));  % Mean u value over all z.
%             mean_u = mean(mean(mean(u))); % Compute mean over entire velocity field
%             u_pert = u - mean_u;
%
%             % Compute the local velocity fluctuations at the reference
%             % location u'(x, y, z) and roving location u'(x+rx, y+ry, z+rz)
%             u_ref = u_pert(xi, yi, z);
%             u_rov = u(xi+rx, yi+ry, z) - mean_u;
%
%             % Compute the turbulence intensity, <u'^2(z)> at the reference location
%             hu2_ref = mean(u_ref.^2);
%
%             % Compute the turbulence intensity term hu2 at the roving location
%             hu2_rov = mean(u_rov .^ 2);
%
%             % Compute the numerator of the correlation function Rij, which is the
%             % spatial average of the product of the velocity fluctuations,
%             % <u'(x,y,z)u'(rx, ry, rz)>
%             num(i,j) = mean(u_ref .* u_rov);
%
%             % Compute the denominator of the correlation function Rij, which is the
%             % product of the turbulence intensity terms at the reference and
%             % roving locations, (<ui'^2(z)><uj'^2(z)>)^0.5
%             den(i,j) = sqrt(hu2_ref * hu2_rov);
%
%             % Compute the correlation function Rij
%             Rij(i,j) = num(i,j) / den(i,j);
%         end
%       end
% end
%%
% figure;
% subplot(2,3,1)
% imagesc(Rij);
% colormap jet;
% colorbar;
% title('Spatial Correlation Function');
% xlabel('x');
% ylabel('y');
%
% subplot(2,3,2)
% imagesc(u(:,:,z));
% colormap jet;
% colorbar;
% xlabel('x');
% ylabel('y');
% title("u field");
%
% subplot(2,3,3)
% plot(x, u(1,:,z));
% xlabel('x');
% ylabel('y');
% title("u_{i=1}");
%
% subplot(2,3,4)
% plot(y, ui1_pert);
% xlabel('x');
% ylabel('y');
% title("u'_{i=1}");
%
% subplot(2,3,5)
% plot(x, uj1_pert);
% xlabel('x');
% ylabel('y');
% title("u'_{j=1}");
%
% subplot(2,3,5)
% contourf(R1);
% xlabel('x');
% ylabel('y');
% colormap jet;
% colorbar;
% title("R11");
