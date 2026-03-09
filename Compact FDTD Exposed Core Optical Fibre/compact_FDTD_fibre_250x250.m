clc; clear all;close all;

% physical constants
C = 299792458.0; MU0 = 1.256637e-6; EPS0 =  8.8541878e-12;

%%% definition
freq = C/1550e-9;%Hz
omega = 2*pi*freq; lambda = C/freq;

%centre of the computational domain
nx=250+40;ny=250+40;
dx = 100e-9;dy = dx;%meters
center_x = round(nx/2);
center_y = round(ny/2);
eps_inf = 1; neff = 1.43;
wavevector=(2*pi/lambda)*neff;

ex=zeros(nx,ny); ey=zeros(nx,ny);
ez=zeros(nx,ny); hx=zeros(nx-1,ny);
hy=zeros(nx,ny-1); hz=zeros(nx,ny);
epsl=zeros(nx,ny);%dielectric constant
sgm=zeros(nx,ny);%electrical conductivity

%definition of the waveguide structure
xind = [1:nx]; yind = [1:ny];
sgm(xind,yind) = 0.0;
epsl(xind,yind) = 1.44^2;

ellipsex = 1.2;
 for xind=1:nx
     for yind=1:ny
         ydist = (center_y-yind-99*2);
         xdist = (center_x-xind);
         dist = sqrt((xdist/ellipsex/1.1)^2+ydist^2);
         if(dist<98*2)
             epsl(xind,yind)=1.0^2;
         end;
     end;
 end;

 for xind=1:nx
     for yind=1:ny
         ydist = (center_y-yind+78*2);
         xdist = (center_x-xind+117*2);
         dist = sqrt((xdist/ellipsex)^2+ydist^2);
         if(dist<98*2)
             epsl(xind,yind)=1.0;
         end;
     end;
 end;

 for xind=1:nx
     for yind=1:ny
         ydist = (center_y-yind+78*2);
         xdist = (center_x-xind-117*2);
         dist = sqrt((xdist/ellipsex)^2+ydist^2);
         if(dist<98*2)
             epsl(xind,yind)=1.0;
         end;
     end;
 end;

excit = epsl;

for plo_k=1:length(wavevector);
    k=wavevector(plo_k); %wave vector

    % Courant stability condition
    dt = 0.95/((C/sqrt(eps_inf))*...
    sqrt(1.0/(dx*dx)+1.0/(dy*dy)+(k/2)^2));
    max_iter = 4000;

    ex=zeros(nx,ny); ey=zeros(nx,ny);
    ez=zeros(nx,ny); hx=zeros(nx-1,ny);
    hy=zeros(nx,ny-1); hz=zeros(nx,ny);

    %old fields for Mur1 ABC
    ex_n=ex;ey_n=ey;ez_n=ez;

    %modal fields
    ex_m=ex;ey_m=ey;ez_m=ez;
    hx_m=hx;hy_m=hy;hz_m=hz;

    %initial condition
    ex(1:nx,:)=excit(1:nx,:)-1;

    %%% FDTD body
    for iteration = 0:max_iter
        SimTime = dt*iteration; iteration

        % Update hx values
        xi=[1:nx-1];yi=[1:ny-1];
        hx(xi,yi) = hx(xi,yi) + (dt/MU0)*(1i*k*ey(xi+1,yi) -...
        (ez(xi+1,yi+1) - ez(xi+1,yi))/dy);

        % Update hy values
        xi=[1:nx-1];yi=[1:ny-1];
        hy(xi,yi) = hy(xi,yi) + (dt/MU0)*((ez(xi+1,yi+1) -...
        ez(xi,yi+1))/dx - 1i*k*ex(xi,yi+1));

        % Update the hz values
        xi=[1:nx-1];yi=[1:ny-1];
        hz(xi,yi)=hz(xi,yi)+(dt/MU0)*((ex(xi,yi+1)-ex(xi,yi))/dy-...
        (ey(xi+1,yi) - ey(xi,yi))/dx);

        % Update ex values:
        xi=[1:nx];yi=[2:ny];
        c1(xi,yi) = 1.0 - (dt.*sgm(xi,yi))./(2.0.*epsl(xi,yi)*EPS0);
        c2(xi,yi) = 1.0 + (dt.*sgm(xi,yi))./(2.0.*epsl(xi,yi)*EPS0);
        ex(xi,yi) = (c1(xi,yi)./c2(xi,yi)).*ex(xi,yi) +...
        (dt./epsl(xi,yi)/EPS0)./c2(xi,yi).* ...
            ((hz(xi,yi) - hz(xi,yi-1))/dy - 1i*k*hy(xi,yi-1));

        % Update ey values
        xi=[2:nx];yi=[1:ny];
        c1(xi,yi) = 1.0 - (dt.*sgm(xi,yi))./(2.0.*epsl(xi,yi)*EPS0);
        c2(xi,yi) = 1.0 + (dt.*sgm(xi,yi))./(2.0.*epsl(xi,yi)*EPS0);
        ey(xi,yi) = (c1(xi,yi)./c2(xi,yi)).*ey(xi,yi) +...
        (dt./epsl(xi,yi)/EPS0)./c2(xi,yi).*...
            (1i*k*hx(xi-1,yi) - (hz(xi,yi) - hz(xi-1,yi))/dx);

        % Update ez values
        xi=[2:nx];yi=[2:ny];
        c1(xi,yi) = 1.0 - (dt.*sgm(xi,yi))./(2.0.*epsl(xi,yi)*EPS0);
        c2(xi,yi) = 1.0 + (dt.*sgm(xi,yi))./(2.0.*epsl(xi,yi)*EPS0);
        ez(xi,yi) = (c1(xi,yi)./c2(xi,yi)).*ez(xi,yi) +...
        (dt./epsl(xi,yi)/EPS0)./c2(xi,yi).*...
            ((hy(xi,yi-1) - hy(xi-1,yi-1))/dx  -...
            (hx(xi-1,yi) - hx(xi-1,yi-1))/dy);

        %Mur1 ABC

        %LEFT WALL:
        xi=[1:nx];yi=[1:ny];
        ez(1,yi)=ez_n(2,yi)+((C*dt/sqrt(eps_inf)-dx)/...
        (C*dt/sqrt(eps_inf)+dx)).*(ez(2,yi)-ez_n(1,yi));
        ey(1,yi)=ey_n(2,yi)+((C*dt/sqrt(eps_inf)-dx)/...
        (C*dt/sqrt(eps_inf)+dx)).*(ey(2,yi)-ey_n(1,yi));

        %RIGHT WALL:
        ez(nx,yi)=ez_n(nx-1,yi)+((C*dt/sqrt(eps_inf)-dx)/...
        (C*dt/sqrt(eps_inf)+dx)).*(ez(nx-1,yi)-ez_n(nx,yi));
        ey(nx,yi)=ey_n(nx-1,yi)+((C*dt/sqrt(eps_inf)-dx)/...
        (C*dt/sqrt(eps_inf)+dx)).*(ey(nx-1,yi)-ey_n(nx,yi));

        %TOP WALL
        ez(xi,1)=ez_n(xi,2)+((C*dt/sqrt(eps_inf)-dy)/...
        (C*dt/sqrt(eps_inf)+dy)).*(ez(xi,2)-ez_n(xi,1));
        ex(xi,1)=ex_n(xi,2)+((C*dt/sqrt(eps_inf)-dy)/...
        (C*dt/sqrt(eps_inf)+dy)).*(ex(xi,2)-ex_n(xi,1));

        %BOTTOM WALL
        ez(xi,ny)=ez_n(xi,ny-1)+((C*dt/sqrt(eps_inf)-dy)/...
        (C*dt/sqrt(eps_inf)+dy)).*(ez(xi,ny-1)-ez_n(xi,ny));
        ex(xi,ny)=ex_n(xi,ny-1)+((C*dt/sqrt(eps_inf)-dy)/...
        (C*dt/sqrt(eps_inf)+dy)).*(ex(xi,ny-1)-ex_n(xi,ny));

        %EDGES
        xi=[1:nx];ex(xi,1)=0.5*(ex(xi,1)+ex(xi,2));
        ex(xi,ny)=0.5*(ex(xi,ny-1)+ex(xi,ny));
        yi=[1:ny];ey(1,yi)=0.5*(ey(1,yi)+ey(2,yi));
        ey(nx,yi)=0.5*(ey(nx-1,yi)+ey(nx,yi));

        %save 'old' fields for Mur
        ex_n=ex;ey_n=ey;ez_n=ez;

        %modal fields
        ld_m = lambda;
        ex_m = ex_m + ex*exp(1i*SimTime*2*pi*C/ld_m);

    end;% end FDTD body

    figure(1);hold on;imagesc(abs(ex_m'));
    shading flat;colormap(jet);colorbar; hold on;
    ydist = (center_y-99*2); xdist = (center_x);
    ellipse(xdist,ydist,98*2*ellipsex*1.1,98*2,'w');
    ydist = (center_y+78*2); xdist = (center_x+117*2);
    ellipse(xdist,ydist,98*2*ellipsex,98*2,'w');
    ydist = (center_y+78*2);xdist = (center_x-117*2);
    ellipse(xdist,ydist,98*2*ellipsex,98*2,'w');
    axis([0 nx 0 ny]); set(gca, 'fontsize', 18); drawnow;

end;%plo_k
