clc; close all; clear all;
N=16384*32; dx=20e-9; c=299792458; dt=dx/(2*c); a=40*dx;
f=(1:N/2)/(dt*N); idx=find(f>0.5*c/a,1):find(f>2.5*c/a,1);
files={'spectrum_NL1.dat','spectrum_LIN1.dat','spectrum_LIN1_NO.dat'};
styles={'r-','b--','k:'};

hold on;
for i=1:3
  d=load(files{i}); out=fft(d(:,1)).*fft(d(:,2));
  plot(1e9*c./f(idx), abs(out(idx)/max(out(idx))), styles{i}, 'linewidth', 1);
end

xlim([400 500]); ylim([0 1]); ylabel('Transmission'); xlabel('Wavelength (nm)');
legend('Nonlinear','Linear','Linear (no defect)','location','northwest');
set(gca,'fontsize',16); print('my_spectrum.png', '-dpng', '-r300');
