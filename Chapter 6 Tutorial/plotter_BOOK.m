clc;
clear all;

diel_profile = load('Eps.dat');

N = 5000;
dx = 20.0e-9;
c  = 299792458.0;
dt = dx/(2*c);

x = 0:1000;
x = x*dx;

for im = 0:100:N

    filename = sprintf('step%d.dat', im);
    loadfile = load(filename);

    Ez = loadfile(:,1);

    clf

    [ax, h1, h2] = plotyy(x/1e-6, diel_profile, x/1e-6, Ez);

    set(h1,'LineWidth',2)
    set(h2,'LineWidth',2,'Color','r')

    set(ax(1),'fontsize',16);
    set(ax(2),'fontsize',16);

    xlabel('x-coordinate (um)');

    ylabel(ax(1),'Dielectric permittivity');
    ylabel(ax(2),'Electric field (V/m)');

    xlim(ax(1),[0 x(end)/1e-6]);
    xlim(ax(2),[0 x(end)/1e-6]);
 %   ylim(ax(2),[-1.5e8 1.5e8]);

    title(sprintf('Time step %d',im));

    drawnow

end

% 2. Save to PNG
% '-dpng' specifies the PNG device
% '-r300' sets the resolution to 300 DPI (optional)
print('my_plot.png', '-dpng', '-r300');
