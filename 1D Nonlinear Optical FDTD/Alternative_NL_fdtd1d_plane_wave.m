% 1D Nonlinear-optical FDTD with alternative PML
clc, close all, clear; pkg load signal %Only for Octave

%*************************************************
%     Fundamental constants
%*************************************************
cc=2.99792458e8; %speed of light in free space
muz=4.0*pi*1.0e-7; %permeability of free space
epsz=1.0/(cc*cc*muz);  %permittivity of free space
Z0 = sqrt(muz/epsz);

freq_p=192e+12; %center frequency of source excitation
lambda_p=cc/freq_p; omega_p=2.0*pi*freq_p;
Ep = 5.1e7;

freq_s=195e+12; %center frequency of source excitation
lambda_s=cc/freq_s; omega_s=2.0*pi*freq_s;
Es = Ep/100;

%*************************************************
%     Grid parameters
%*************************************************
ie=1200;  %number of grid cells in x-direction
ienl1 = 100; ienl2 = 1100;

ib=ie+1; is=15;   %location of z-directed hard source
dx=1e-8;        %space increment of square lattice
dt=0.5*dx/cc;   %time step
nmax=round(4e-12/dt);  %total number of time steps
iebc=20;   %thickness of left and right PML region

rmax=0.00001; orderbc=2; ibbc=iebc+1;
iefbc=ie+2*iebc; ibfbc=iefbc+1;

is1=7; %connecting boundaries
is2=ie-is1-1;

%*************************************************
%     Material parameters
%*************************************************
media=1; xx3=[1e-18 1e-18]; eps=[1.0 1.0];
sig=[0.0 0.0];%sig=[0.0 1.0e+7];
mur=[1.0 1.0]; sim=[0.0 0.0];

%*************************************************
%     Wave excitation
%*************************************************
ey_inc=zeros(ib,1); hz_inc=zeros(ie,1);
ey_low_m1=0; ey_low_m2=0;
ey_high_m1=0; ey_high_m2=0;

%*************************************************
%     Field arrays
%*************************************************
ey=zeros(ib,1); ey_old=zeros(ib,1); hz=zeros(ie,1);

eybcl=zeros(iebc,1); %fields in left PML region
hzxbcl=zeros(iebc,1); hzybcl=zeros(iebc,1);

eybcr=zeros(ibbc,1);  %fields in right PML region
hzxbcr=zeros(iebc,1);  hzybcr=zeros(iebc,1);

%*************************************************
%     Update coefficients
%*************************************************
for i=1:media
  eaf  =dt*sig(i)/(2.0*epsz*eps(i));
  ca(i)=(1.0-eaf)/(1.0+eaf);
  cb(i)=dt/epsz/eps(i)/dx/(1.0+eaf);
  haf  =dt*sim(i)/(2.0*muz*mur(i));
  da(i)=(1.0-haf)/(1.0+haf);
  db(i)=dt/muz/mur(i)/dx/(1.0+haf);
end

%*************************************************
%    Main grid initialisation
%*************************************************
caey(1:ib)=ca(1); cbey(1:ib)=cb(1);
dahz(1:ie)=da(1); dbhz(1:ie)=db(1);

%  Add non-linear region
x3=zeros(ib,1); x3(ienl1:ienl2) = xx3(1);

%*************************************************
%     PML regions
%*************************************************
delbc=iebc*dx;
sigmam=-log(rmax/100.0)*epsz*cc*(orderbc+1)/(2*delbc);
bcfactor=eps(1)*sigmam/(dx*(delbc^orderbc)*(orderbc+1));

%  LEFT
caeybcl(1)=1.0; cbeybcl(1)=0.0;
for i=2:iebc
  x1=(iebc-i+1.5)*dx; x2=(iebc-i+0.5)*dx;
  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
  ca1=exp(-sigmax*dt/(epsz*eps(1)));
  cb1=(1-ca1)/(sigmax*dx);
  caeybcl(i)=ca1; cbeybcl(i)=cb1;  caeybcf(i)=ca1;
  cbeybcf(i)=cb1; caeybcb(i)=ca1; cbeybcb(i)=cb1;
end

sigmax=bcfactor*(0.5*dx)^(orderbc+1);
ca1=exp(-sigmax*dt/(epsz*eps(1)));
cb1=(1-ca1)/(sigmax*dx);
caey(1)=ca1; cbey(1)=cb1;
caeybcf(iebc+1)=ca1; cbeybcf(iebc+1)=cb1;
caeybcb(iebc+1)=ca1; cbeybcb(iebc+1)=cb1;

for i=1:iebc
  x1=(iebc-i+1)*dx; x2=(iebc-i)*dx;
  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
  sigmaxs=sigmax*(muz/(epsz*eps(1)));
  da1=exp(-sigmaxs*dt/muz);
  db1=(1-da1)/(sigmaxs*dx);
  dahzxbcl(i)=da1; dbhzxbcl(i)=db1;
  dahzxbcf(i)=da1; dbhzxbcf(i)=db1;
  dahzxbcb(i)=da1; dbhzxbcb(i)=db1;
  caexbcl(i)=ca(1); cbexbcl(i)=cb(1);
  dahzybcl(i)=da(1); dbhzybcl(i)=db(1);
end

%  RIGHT

caeybcr(ibbc)=1.0; cbeybcr(ibbc)=0.0;
for i=2:iebc
  x1=(i-0.5)*dx; x2=(i-1.5)*dx;
  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
  ca1=exp(-sigmax*dt/(epsz*eps(1)));
  cb1=(1-ca1)/(sigmax*dx);
  caeybcr(i)=ca1; cbeybcr(i)=cb1;
  caeybcf(i+iebc+ie)=ca1; cbeybcf(i+iebc+ie)=cb1;
  caeybcb(i+iebc+ie)=ca1; cbeybcb(i+iebc+ie)=cb1;
end

sigmax=bcfactor*(0.5*dx)^(orderbc+1);
ca1=exp(-sigmax*dt/(epsz*eps(1)));
cb1=(1-ca1)/(sigmax*dx);
caey(ib)=ca1; cbey(ib)=cb1;
caeybcf(iebc+ib)=ca1; cbeybcf(iebc+ib)=cb1;
caeybcb(iebc+ib)=ca1; cbeybcb(iebc+ib)=cb1;

for i=1:iebc
  x1=i*dx; x2=(i-1)*dx;
  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
  sigmaxs=sigmax*(muz/(epsz*eps(1)));
  da1=exp(-sigmaxs*dt/muz);
  db1=(1-da1)/(sigmaxs*dx);
  dahzxbcr(i) = da1; dbhzxbcr(i) = db1;
  dahzxbcf(i+ie+iebc)=da1; dbhzxbcf(i+ie+iebc)=db1;
  dahzxbcb(i+ie+iebc)=da1; dbhzxbcb(i+ie+iebc)=db1;
  caexbcr(i)=ca(1); cbexbcr(i)=cb(1);
  dahzybcr(i)=da(1); dbhzybcr(i)=db(1);
end

%*************************************************
% MAIN FDTD LOOP
%*************************************************
for n=1:nmax

%1D buffer
ey_inc(2:ib-1)=ey_inc(2:ib-1)-(dt/(epsz*dx)).*...
                                   (hz_inc(2:ie)-hz_inc(1:ie-1));
ey_inc(1)=ey_low_m2; ey_low_m2=ey_low_m1;
ey_low_m1=ey_inc(2); ey_inc(ib)=ey_high_m2;
ey_high_m2=ey_high_m1; ey_high_m1=ey_inc(ib-1);

ey_inc(3)=ey_inc(3)+Ep*sin(omega_p*dt*n) +...
                                            Es*sin(omega_s*dt*n);

ey_old = ey;%    save old ey fields

ey(2:ie)=caey(2:ie)'.*ey(2:ie)+...
           cbey(2:ie)'.*(hz(1:ie-1)-hz(2:ie));

ey(2:ie) = ey(2:ie)./(eps(1) +...
                              x3(2:ie).*abs(ey_old(2:ie)).^2);

%Ey field correction
ey(is1)=ey(is1)+(dt/(epsz*dx)).*hz_inc(is1-1);
ey(is2)=ey(is2)-(dt/(epsz*dx)).*hz_inc(is2);

% Ey in PML regions
% LEFT
eybcl(2:iebc)=caeybcl(2:iebc)'.*eybcl(2:iebc)-...
  cbeybcl(2:iebc)'.*(hzxbcl(2:iebc)+hzybcl(2:iebc)-...
                      hzxbcl(1:iebc-1)-hzybcl(1:iebc-1));
ey(1)=caey(1)'.*ey(1)-...
  cbey(1)'.*(hz(1)-hzxbcl(iebc)-hzybcl(iebc));

% RIGHT
eybcr(2:iebc)=caeybcr(2:iebc)'.*eybcr(2:iebc)-...
  cbeybcr(2:iebc)'.*(hzxbcr(2:iebc)+hzybcr(2:iebc)-...
                      hzxbcr(1:iebc-1)-hzybcr(1:iebc-1));
ey(ib)=caey(ib)'.*ey(ib)-...
  cbey(ib)'.*(hzxbcr(1)+hzybcr(1)- hz(ie));

%  Hz in main grid
hz_inc(1:ie)=hz_inc(1:ie)-(dt/(muz*dx)).*...
                                  (ey_inc(2:ib)-ey_inc(1:ib-1));

hz(1:ie)=dahz(1:ie)'.*hz(1:ie)+...
                              dbhz(1:ie)'.*(ey(1:ie)-ey(2:ib));

%hz field correction
hz(is1-1)=hz(is1-1)+(dt/(muz*dx)).*ey_inc(is1);
hz(is2)=hz(is2)-(dt/(muz*dx)).*ey_inc(is2);

%  Additional updates in PML regions
%  LEFT
hzxbcl(1:iebc-1)=dahzxbcl(1:iebc-1)'.*hzxbcl(1:iebc-1)-...
  dbhzxbcl(1:iebc-1)'.*(eybcl(2:iebc)-eybcl(1:iebc-1));
hzxbcl(iebc)=dahzxbcl(iebc)'.*hzxbcl(iebc)-...
  dbhzxbcl(iebc)'.*(ey(1)-eybcl(iebc));

% RIGHT
hzxbcr(2:iebc)=dahzxbcr(2:iebc)'.*hzxbcr(2:iebc)-...
  dbhzxbcr(2:iebc)'.*(eybcr(3:ibbc)-eybcr(2:iebc));
hzxbcr(1)=dahzxbcr(1)'.*hzxbcr(1)-...
  dbhzxbcr(1)'.*(eybcr(2)-ey(ib));

% Detector (to save results)
ey_time(n) = ey(ienl2+2);
if mod(n,5000)==0;
    timestep=int2str([n nmax])
end;

end%End of the main FDTD loop

% Fourier transformation
w1 = window(@gausswin,nmax,5);
ey_dft=fft(w1'.*ey_time(1:nmax));

T_dft=1:max(nmax)/2;freq_dft=T_dft/(dt*nmax);
ld=cc./freq_dft;om_dft=2*pi.*freq_dft;
imin=1;imax=1;ld_min=cc/175e12;ld_max=cc/210e12;
[tmp imin] = min(abs(ld_dft-ld_min));
[tmp imax] = min(abs(ld_dft-ld_max));
Port2_dft = abs(ey_dft(imin:imax));
mmax_dft = max(Port2_dft);

figure;fntsz=24;set(gca,'fontsize',fntsz);hold on;
plot(cc./ld_dft(imin:imax)/1e12,...
10*log10(Port2_dft/mmax_dft),'r-','linewidth',1);
grid on; axis([175 210 -80 0]);box on;
xlabel('Frequency (THz)'); ylabel('Spectrum (dB)');
