clc;close all;clear all;

%% Constants
pi_val = pi;
melec = 9.1093837e-31;
hbar = 1.054571817e-34;
elcharge = 1.60217663e-19;
kf = 6.241509e18;

DomainSize = 4001;
kstart = floor(DomainSize/2);
kcentre = floor(kstart/2);
NSTEPS = 100000;
ddx = 0.1e-11;
ra = 0.125;
lambda_val = 1.6e-10;
sigma = 1.6e-10;
U0 = 600.0; % eV

%% Arrays
PsiReal = zeros(1, DomainSize);
PsiImag = zeros(1, DomainSize);
vp = zeros(1, DomainSize);

%% Variables
dt = 0.25 * (melec / hbar) * ddx^2;
L = DomainSize * ddx;

%% Generate the potential well
for k = 1:DomainSize
    xx = (k-1) * ddx;
    Vpot = (4.0 * U0 / L^2) * (xx^2) - (4.0 * U0 / L) * xx;
    vp(k) = Vpot * elcharge;
end

%% Add a barrier to the potential well
for k = (kstart-10):(kstart+10)
    vp(k) = vp(k) + 400.0 * elcharge;
end

%% Initialise the wavefunction
ptot = 0;
for k = 2:(DomainSize-2)
    PsiReal(k) = cos(2*pi_val*ddx*(k-kcentre)/lambda_val) * ...
                 exp(-0.5 * ((ddx*(k-kcentre))/sigma)^2);
    PsiImag(k) = sin(2*pi_val*ddx*(k-kcentre)/lambda_val) * ...
                 exp(-0.5 * ((ddx*(k-kcentre))/sigma)^2);
    ptot = ptot + PsiReal(k)^2 + PsiImag(k)^2;
end

%% Normalise the wavefunction
norm_factor = sqrt(ptot);
PsiReal(2:DomainSize-2) = PsiReal(2:DomainSize-2) / norm_factor;
PsiImag(2:DomainSize-2) = PsiImag(2:DomainSize-2) / norm_factor;

%% Time evolution
for n = 1:(NSTEPS-1)
    PsiReal(2:end-1) = PsiReal(2:end-1) - ra*(PsiImag(3:end) - 2*PsiImag(2:end-1) + PsiImag(1:end-2)) + ...
                       (dt/hbar) * vp(2:end-1) .* PsiImag(2:end-1);

    PsiImag(2:end-1) = PsiImag(2:end-1) + ra*(PsiReal(3:end) - 2*PsiReal(2:end-1) + PsiReal(1:end-2)) - ...
                       (dt/hbar) * vp(2:end-1) .* PsiReal(2:end-1);
end

%% Energy calculation
KeReal = 0; KeImag = 0; Pe = 0;
for k = 2:(DomainSize-2)
    LapReal = PsiReal(k+1) - 2*PsiReal(k) + PsiReal(k-1);
    LapImag = PsiImag(k+1) - 2*PsiImag(k) + PsiImag(k-1);

    KeReal = KeReal + PsiReal(k)*LapReal + PsiImag(k)*LapImag;
    KeImag = KeImag + PsiReal(k)*LapImag - PsiImag(k)*LapReal;
    Pe = Pe + vp(k) * (PsiReal(k)^2 + PsiImag(k)^2);
end

kine = 0.5 * (hbar / melec) * (hbar / ddx^2) * sqrt(KeReal^2 + KeImag^2);
fprintf('Kinetic Energy (ke): %.16e, Potential Energy (pe): %.16e\n', kine*kf, Pe*kf);

%% Compute values for plotting
x_values = (0:(DomainSize-1)) * ddx; % Convert indices to positions
vp_scaled = vp / elcharge;           % Potential in eV
wavefunction_amplitude = sqrt(PsiReal.^2 + PsiImag.^2); % Combined amplitude

%% Plot using plotyy
[ax, h1, h2] = plotyy(x_values, wavefunction_amplitude, x_values, vp_scaled);
set(h1, 'LineWidth', 1.5, 'Color', 'b');
set(h2, 'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'g');
xlabel('Position (m)');
ylabel(ax(1), 'Wavefunction Amplitude');
ylabel(ax(2), 'Potential (eV)');
grid on;
title('Wavefunction and Potential');
legend([h1, h2], 'Wavefunction Amplitude','Potential (eV)');
