% MChinita (2/14/2023)
% This code is a template for the auto-correlation function. It assumes
% the following variables have been loaded:
% z (vertical grid), u (3D zonal wind component), v, and w 

% nz, nx, and ny are the number of grid points in the 3D domain, e.g., 
% from Chinita et al 2022, run B = 150 x 256 x 256 grid points, so nz =
% 150, nx = 256, ny = 256 grid points.

% Instantaneous u-,v-, and w-components 
u2=zeros(nx,ny,nz);
v2=zeros(nx,ny,nz);
w2=zeros(nx,ny,nz);

% The instantaneous fields of UConn's LES model are in a Arakawa C (staggered) grid, 
% e.g., u = nz x nx+1 x ny, so here it converts to nz x nx x ny
for i=1:nx
   for j=1:ny
       for k=1:nz
           u2(i,j,k)=0.5*(u(k,i,j)+u(k,i+1,j));
           v2(i,j,k)=0.5*(v(k,i,j)+v(k,i,j+1));
           w2(i,j,k)=0.5*(w(k,i,j)+w(k+1,i,j));
       end
   end 
end


%% Compute turbulent part of the flow (i.e., fluctuations)
u2_mean= squeeze(mean(mean(u2)));
v2_mean= squeeze(mean(mean(v2)));
w2_mean= squeeze(mean(mean(w2)));

for i=1:nx
   for j=1:ny
       for k=1:nz
           uf(i,j,k)=u2(i,j,k)-u2_mean(k);
           vf(i,j,k)=v2(i,j,k)-v2_mean(k);
           wf(i,j,k)=w2(i,j,k)-w2_mean(k);
       end
   end 
end

% Calculate variance
uf2 = uf.*uf;
variance_uf2 = squeeze(mean(mean(uf2)));

vf2 = vf.*vf;
variance_vf2 = squeeze(mean(mean(vf2)));

wf2 = wf.*wf;
variance_wf2 = squeeze(mean(mean(wf2)));

% Two-point autocorrelation function using rx by ry grid points 
rx = 20; % [-20:20] grid points
ry = 20;

% Let's create a frame of zeros around the turbulent quantities
uf_framed = zeros(nx+rx*2,ny+rx*2,hz); 
uf_framed(rx+1:end-rx,ry+1:end-ry,:) = uf(:,:,1:hz); % center
% The domain is doubly periodic so let's replace the zeros with the values
% of the opposite size
uf_framed(1:rx,ry+1:end-ry,:) = uf(end-rx+1:end,:,1:hz); % left side uses the right side of uf
uf_framed(end-rx+1,ry+1:end-ry,:) = uf(1:rx,:,1:hz); % right side uses the left side of uf
uf_framed(rx+1:end-rx,1:ry,:) = uf(:,end-ry+1:end,1:hz); % top uses the bottom of uf
uf_framed(rx+1:end-rx,end-ry+1,:) = uf(:,1:ry,1:hz); % top uses the bottom of uf

% Vertical level at the top of the PBL (use method of preference)
hz = 28;

vecx = (-rx:rx); 
vecy = (-rx:rx);

for izz=1:hz
    iz = izz;
    cc = 0;
    for i = 1:size(vecx,2)
        for j = 1:size(vecy,2)
            for ii=rx+1:nx+rx
                for jj = ry+1:ny+ry
                    cc = cc+1;
                    helpR(cc) = uf_framed(ii,jj,iz).*uf_framed(ii+vecx(i),jj+vecy(j),iz);
                end
            end
            Ruu(i,j) = sum(helpR)./cc;
            Ruu_norm(i,j) = Ruu(i,j)./variance_uf2(iz);
            cc = 0;
        end
    end
        
    Ruu_norm_z(izz,:,:) = Ruu_norm;
end

% Integral length scale;  dx is the grid resolution
for izz = 1:hz
    L1(izz) = sum(Ruu_norm_z(izz,:,rx+1).*dx);
end

% Figure
hsbl = z(28); % normalize by PBL height

[X,Y] = meshgrid((vecy.*dz)./hsbl,(vecx.*dz)./hsbl);
hfig=figure; 
[C,h]=contour(Y,X,squeeze(Ruu_norm_z(6,:,:)),'LevelList',(0.1:0.2:0.7),'LineColor','k','ShowText','on'); hold on
clabel(C,h,'Interpreter','Latex','FontSize',12)
% If there is overturning motions, plot those using dashed lines
[C,h]=contour(Y,X,squeeze(Ruu_norm_z(6,:,:)),'--','LevelList',(-0.1),'LineColor','k','ShowText','on')
clabel(C,h,'Interpreter','Latex','FontSize',12)
xlim([-0.3 0.3]); ylim([-0.3 0.3]);
ylabel('$r_y$','Interpreter','Latex','Fontsize',18); 
xlabel('$r_x$','Interpreter','Latex','Fontsize',18); 
text(4.5,5,'$R_{11}$', 'Interpreter', 'latex','Fontsize',12);
title('z/h = 0.20','Interpreter', 'latex','Fontsize',12)
set(gca,'FontSize',12,'TickLabelInterpreter','Latex')
pbaspect([1 1 1])

% Repeat for vv and ww 