
% 2D-analytical model of a single layer with a rectangular slot
% Using analytical Fourier transform
% Johan Lundgren 2019

% Fundamental constants
c = 299792458;
mu0 = 4*pi*10^(-7);
eps0 = 1/(c^2*mu0);
Z0 = sqrt(mu0/eps0);
Y0 = 1/Z0;

px = 10e-3; % Periodicity of the structure in x-dir. [m]
py = 10e-3; % Periodicity of the structure in y-dir. [m]
wx = 1.2e-3; % Width of slit in x-dir [m]
wy = 6e-3; % Width of slit in y-dir [m]
% d1 = 1.5e-3; % Thickness of layer 1 [m]
% d2 = 6e-3; % Thickness of layer 2 [m]
% d3 = 1.5e-3; % Thickness of layer 3 [m]
% d = [d1 d2 d3];
er1 = 1; % relative permittivity in media 1
er2 = 1; % relative permittivity in media 2
er3 = 1; % relative permittivity in media 3
er = [er1 er2 er3];

f = linspace(16e9,32e9,1002); % Frequency band [Hz]
omega = 2*pi*f;
k0 = omega/c;

% Number of modes used
N = 50;
n = linspace(1,N,N);
n = [flip(-n) 0.0000001 n];
M = 50;
m = linspace(1,M,M);
m = [flip(-m) 0.0000001 m];

[N_m,M_m]=meshgrid(n,m);

kn = 2*pi/px*N_m;
km = 2*pi/py*M_m;

% Plotting field-profile
steps_x = 2^10;
lim = px/2;
dx = 2*lim/steps_x;
x = -lim:dx:lim-dx;

steps_y = 2^10;
lim = py/2;
dy = 2*lim/steps_y;
y = -lim:dy:lim-dy;

[X,Y]=meshgrid(x,y);

% TM
Field = cos(pi*X/wx).*(1 - (2*X/wx).^2).^(-1/2).*rectpuls(Y/wy).*rectpuls(X/wx);

figure(1)
imagesc(x/1e-3,y/1e-3,Field)
%xlim([-px/2e-3 px/2e-3])
%ylim([-py/2e-3 py/2e-3])
colormap jet;
colorbar
xlabel('x [mm]')
ylabel('y [mm]')
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
pbaspect([1 1 1])

% Fourier transform of field-profile
f_t = @(kx,ky) pi*wy/4*(besselj(0,wy/2*abs(ky + pi/wy))...
    + besselj(0,wy/2*abs(ky - pi/wy)))*wx.*sinc(wx*kx/(2*pi));

% Transformer ratios
N0 = f_t(0,0)/sqrt(2);
Nh_TM = f_t(kn,km).*(kn./sqrt(kn.^2 + km.^2));
Nh_TE = f_t(kn,km).*(km./sqrt(kn.^2 + km.^2));

% Equivalent admittance for each frequency
Yeq = zeros(1,length(f));
for i = 1:length(f)
    betah0 = -1i*sqrt(abs(kn).^2 + abs(km).^2 - er(1)*k0(i)^2);
    betah1 = -1i*sqrt(abs(kn).^2 + abs(km).^2 - er(2)*k0(i)^2);

    Yh0_TM = omega(i)*eps0*er(1)./betah0;
    Yh1_TM = omega(i)*eps0*er(2)./betah1;
    
    Yh0_TE = betah0./(omega(i)*mu0);
    Yh1_TE = betah1./(omega(i)*mu0);

    YhL_TM = Yh0_TM;
    YhR_TM = Yh1_TM;

    YhL_TE = Yh0_TE;
    YhR_TE = Yh1_TE;

    Yeq_TM(i) = sum(sum(abs(Nh_TM).^2.*(YhL_TM + YhR_TM))) - abs(N0)^2*Y0 - abs(N0)^2*Y0;
    Yeq_TE(i) = sum(sum(abs(Nh_TE).^2.*(YhL_TE + YhR_TE))) - abs(N0)^2*Y0 - abs(N0)^2*Y0;
end
Yeq = Yeq_TM + Yeq_TE;

% Calculation of scattering paramters
S11 = (abs(N0)^2*2*Y0 - abs(N0)^2*2*Y0 - Yeq)./(abs(N0)^2*2*Y0 + abs(N0)^2*2*Y0 + Yeq);
S21 = 2*abs(N0)^2*2*Y0./(abs(N0)^2*2*Y0 + abs(N0)^2*2*Y0 + Yeq);


figure(2)
plot(f/1e9,abs(S11),'r','Linewidth',1.2)
xlim([16 32])
grid ON
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
xlabel('Frequency (GHz)')
ylabel('Reflection coefficient')

figure(3)
plot(f/1e9,abs(S21),'r','Linewidth',1.2)
xlim([16 32])
grid ON
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
xlabel('Frequency (GHz)')
ylabel('Transmission coefficient')




