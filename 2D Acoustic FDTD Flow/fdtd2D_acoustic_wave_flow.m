clc;clear all; close all;

%% user modifiable parameters
nx = 201; ny = 200; %set number of nodes in x and y directions

%WIND
vx_max = 300.0; vy_max = 300.0; %flow velocity vector field
c_water = 1501.0; c_air = 343.0; %speed of sound m/s
cmax = max([c_water c_air])+sqrt(vx_max^2+vy_max^2);
rho0_water = 1000.0; rho0_air = 0.1*rho0_water; %density kg/m^3

% excitation parameters
freq = 100.0; omega = 2*pi*freq; lambda = cmax/freq;

%distance between nodes
dx = 0.5; dy = dx; max_iter = 500;

%% Initialise the p and u fields and set coefficients
p = zeros(nx,ny); wx = zeros(nx-1,ny); wy = zeros(nx,ny-1);
vx = vx_max*ones(nx-1,ny); vy = vy_max*ones(nx,ny-1); %flow velocity

b = zeros(nx,ny); %buoyancy
kappa = zeros(nx,ny); %adiabatic bulk modulus

p_n = zeros(nx,ny); wx_n = zeros(nx-1,ny); wy_n = zeros(nx,ny-1);
p_mod = zeros(nx,ny);wx_mod = zeros(nx-1,ny);wy_mod= zeros(nx,ny-1);

dt = dx/(sqrt(ndims(p))*cmax); %time step
dt = 0.2*dt;%time step correction (experimental)

%set the material parameters
center_x = ceil(nx/2); center_y = ceil(ny/2);

%water
kappa = ones(nx,ny)*rho0_water*c_water^2;
b = ones(nx,ny)/rho0_water;

%% FDTD loop
for n = 1:max_iter;
    % Update the wx values:
    i=[2:nx-1];j=[2:ny-1];
    wx(i,j) = wx(i,j) - (dt*(b(i+1,j) + b(i,j))/2.0).*...
        ((p(i,j) - p(i-1,j))/dx) ...
        - dt*(wx(i,j).*(vx(i,j) - vx(i-1,j))/dx +...
        wy(i,j).*(vx(i,j) - vx(i,j-1))/dy) ...
        - dt*(vx(i,j).*(wx(i,j) - wx(i-1,j))/dx +...
        vy(i,j).*(wx(i,j) - wx(i,j-1))/dy);

    % Update the wy values
    i=[2:nx-1];j=[2:ny-1];
    wy(i,j) = wy(i,j) - (dt*(b(i,j+1) + b(i,j))/2.0).*...
        (-(p(i,j) - p(i,j-1))/dy) ...
        - dt*(wx(i,j).*(vy(i,j) - vy(i-1,j))/dx +...
        wy(i,j).*(vy(i,j) - vy(i,j-1))/dy) ...
        - dt*(vx(i,j).*(wy(i,j) - wy(i-1,j))/dx +...
        vy(i,j).*(wy(i,j) - wy(i,j-1))/dy);

	i=[1:nx-2];j=[1:ny-2];
	vx1 = (vx(i+1,j)+vx(i,j))/2.0;vy1 = (vy(i,j+1)+vy(i,j))/2.0;
	p(i,j)= p(i,j)-dt*kappa(i,j).*((wx(i+1,j) - wx(i,j))/dx -...
        (wy(i,j+1) - wy(i,j))/dy) ...
        - dt*(vx1.*(p(i+1,j) - p(i,j))/dx +...
        vy1.*(p(i,j+1) - p(i,j))/dy);

    % excittaion source
    p(ceil(nx/2),ceil(ny/2)) = p(ceil(nx/2),ceil(ny/2)) +...
    sin(2.0*pi*n*dt*freq); %insert hard source

    %Mur1 ABC
    c = c_water;%! water at the infinity
    i=[1:nx];j=[1:ny];
    p(1,j)=p_n(2,j)+((c*dt-dx)/(c*dt+dx)).*(p(2,j)-p_n(1,j));
    p(nx-1,j)=p_n(nx-2,j)+((c*dt-dx)/(c*dt+dx)).*...
    (p(nx-2,j)-p_n(nx-1,j));
    p(i,1)=p_n(i,2)+((c*dt-dy)/(c*dt+dy)).*(p(i,2)-p_n(i,1));
    p(i,ny-1)=p_n(i,ny-2)+((c*dt-dy)/(c*dt+dy)).*...
    (p(i,ny-2)-p_n(i,ny-1));

    %save 'old' fields for Mur's ABC
    wx_n=wx;wy_n=wy;p_n = p;

  if ( mod(n,50) == 0 && n>1)
    figure(1);subplot(1,1,1);
    imagesc(p(:,:)/0.1); colorbar; caxis([-1 1]);
	hold on; set(gca,'FontSize', 22);
    axis image; colormap(redblue); title(n); drawnow;
  end;

  p_mod = p_mod + p*exp(1i*n*dt*2*pi*freq);
  wx_mod = wx_mod + wx*exp(1i*n*dt*2*pi*freq);
  wy_mod = wy_mod + wy*exp(1i*n*dt*2*pi*freq);

  p_det(n) = p(ceil(3*nx/4),ceil(3*ny/4));

end;% end FDTD body
