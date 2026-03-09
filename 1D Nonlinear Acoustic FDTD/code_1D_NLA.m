clc; close all; clear all;

KE = 250;%size of the computtaional domain

%memory allocation for the arrays
g2 = ones(KE,1); g3 = ones(KE,1);
f2 = ones(KE,1); f3 = ones(KE,1);
ga = ones(KE,1); gb = ones(KE,1);

p = zeros(KE,1); epsilon = ones(KE,1);
rhoa = zeros(KE,1); beta = zeros(KE,1);
tmp_p = zeros(KE,1); pn_m1 = zeros(KE,1);
u = zeros(KE,1);
%material parameters
rho0 = ones(KE,1)*1000;
c = ones(KE,1)*1540; C = zeros(KE,1);
shear_visc = zeros(KE,1);
bulk_visc = zeros(KE,1);
cp = zeros(KE,1); cv = zeros(KE,1);
kappa = zeros(KE,1); delta_test = zeros(KE,1);

%PML-related parameters
npml = 10; n_pml = npml;

for k=0:n_pml-1
    xxn =    (npml-k)/npml; xn = 0.33*xxn^3;
    g2(k+1) = 1/(1+xn); g2(KE-k) = 1/(1+xn);
    g3(k+1) = (1-xn)/(1+xn);
    g3(KE-k) = (1-xn)/(1+xn);

    xxn = (npml-k-0.5)/npml; xn = 0.33*xxn^3;
    f2(k+1) = 1/(1+xn); f2(KE-k-1) = 1/(1+xn);
    f3(k+1) = (1-xn)/(1+xn);
    f3(KE-k-1) = (1-xn)/(1+xn);
end;

%nonlinear media parameters
kc = KE/2; beta(npml+11:KE-npml-11) = 5;

%some useful constants
Na = 6.022e23; eps0 = 8.854e-12;
M = 18e-3; alpha = 1.62e-40;%water
tmp_c = Na*alpha/3/eps0/M;

%space and time increments
cmax = max(c); dz =  18e-6; dt = 0.499*dz/cmax;

k=[1:KE]; ga(k) = dt.*rho0(k).*c(k).^2/dz;
gb(k) = dt./rho0(k)/dz;

shear_visc(k) = 1.002e-3;%Pa sec
bulk_visc(k) = 3.09e-3;%Pa sec
cp(k) = 4.182e3;%J/(kg K)
cv(k) = cp(k)/1.33; kappa(k) = 0.597;%W/(m K)
C(k) = kappa(k).*(1./cv(k)-1./cp(k))./rho0(k)./c(k).^4;

%main FDTD loop
T = 0;for n=1:4000
    T = T + 1; pn_m1 = tmp_p;  tmp_p = p;

    k=[2:npml+10];
    p(k)=g3(k).*p(k)+g2(k).*ga(k).*(u(k-1)-u(k));

    k=[npml+11:KE-npml-11];

    p(k)=-((c(k).^4*dz.*pn_m1(k)-...
    2*c(k).^4*dz.*p(k)).*rho0(k).*C(k)+...
        c(k).^4*dt^2.*rho0(k).^2.*(u(k-1)-u(k))+...
        c(k).^2*dt*dz.*p(k).*rho0(k)-...
        beta(k)*dt*dz.*p(k).*pn_m1(k)+...
        beta(k)*dt*dz.*p(k).^2)./...
        (c(k).^4*dz.*rho0(k).*C(k)-c(k).^2*dt*dz.*rho0(k)-...
        beta(k)*dt*dz.*pn_m1(k)+beta(k)*dt*dz.*p(k));

    k=[KE-npml-10:KE];
    p(k)=g3(k).*p(k)+g2(k).*ga(k).*(u(k-1)-u(k));

    pulse = 5e6*sin(2*pi*2e6*T*dt); p(n_pml+5) = pulse;

    k=[1:npml+10];u(k)=f3(k).*u(k)+f2(k).*gb(k).*(p(k)-p(k+1));

    k=[npml+11:KE-npml-11];
    u(k)=u(k)+gb(k).*(p(k)-p(k+1))-...
        ((bulk_visc(k)+4*shear_visc(k)/3)./...
        rho0(k).^2./c(k).^2/dz).*...
        (p(k+1)-p(k)-tmp_p(k+1)+tmp_p(k));

    k=[KE-npml-10:KE-1];
    u(k)=f3(k).*u(k)+f2(k).*gb(k).*(p(k)-p(k+1));

    p_det(n) = p(KE-npml-5);

    k=[1:KE-1];
    rhoa(k)=p(k)./c(k).^2-(beta(k)-1)./rho0(k)./c(k).^4.*p(k).^2-...
    C(k).*(p(k)-tmp_p(k))/dt;
    epsilon(k) = (1+2*tmp_c.*(rho0(k)+rhoa(k)))./...
    (1-tmp_c.*(rho0(k)+rhoa(k)));

    %visualisation in the time domain
    if ( mod(n,50) == 0)
        figure(1);subplot(2,1,1);plot(p,'r-');
        axis([npml KE-npml -5e6 5e6]);title(n);
        xlabel('x-coordinate'); ylabel('Pressure (Pa)');
        figure(1);subplot(2,1,2);plot(epsilon,'r-');
        xlabel('x-coordinate'); ylabel('Dielectric permittivity');
        axis([npml KE-npml min(epsilon(1:KE-1)) ...
                                      max(epsilon(1:KE-1))]);
        title(n);drawnow;
    end;

end;%fdtd

%FFT analysis
N = length(p_det); w1 = ones(N,1);
fft_p = abs(fft(w1'.*p_det));
T=1:N/2;f=T/(dt*N);ld=c(1)./f;
imin=1;imax=1;ld_min=c(1)/0.1e1;ld_max=c(1)/10e6;
[tmp imin] = min(abs(ld-ld_min));
[tmp imax] = min(abs(ld-ld_max));

figure(10); set(gca,'fontsize',32);
plot(c(1)./ld(imin:imax)/2e6,...
10*log10(abs(fft_p(imin:imax))/max(abs(fft_p(imin:imax)))),...
'r-','linewidth',1);hold on;
ylabel('Normalised amplitude, dB');
xlabel('f/f_{0}'); axis square;


