% Implementation of code for stacked equal 2D-layers with field profile
% imported from CST
% Johan Lundgren 2019


c = 299792458;
mu0 = 4*pi*10^(-7);
eps0 = 1/(c^2*mu0);
Z0 = sqrt(mu0/eps0);
Y0 = 1/Z0;

%%%%%%% Frequencies to caluculate %%%%%%%
f = linspace(16e9,32e9,1002); % Frequency band [Hz]
omega = 2*pi*f;
k0 = omega/c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%% General parameters for the layers %%%%%%%
px = 10e-3; % Periodicity of the structure in x-dir. [m]
py = 10e-3; % Periodicity of the structure in y-dir. [m]
path_TM = 'field_profiles/e-field (f=22.5)_TM.txt'; % Path to file w. inc. x-pol.
path_TE = 'field_profiles/e-field (f=22.5)_TE.txt'; % Path to file w. inc. y-pol.
%path_TM = 'TM';
%path_TE = 'TE';
n_layer = 2; % Number of layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Layer 1 to 2 %%%%%%%
d1 = 3e8/(22.5e9*10); % Distance from layer 1 to 2
er1 = 1; % Perimittivity of material between layer 2 and 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Layer 2 to 3 %%%%%%%
d2 = 3e8/(22.5e9*3); % Distance from layer 2 to 3
er2 = 1; % Perimittivity of material between layer 2 and 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_v = [d1 d2];
er_v = [er1 er2];
p_v = [px py];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SMatrix = ComputeSLayer(f, path_TM, p_v, d_v, er_v, n_layer);
S11_YY = squeeze(SMatrix(1,1,:));
S21_YY = squeeze(SMatrix(2,1,:));

%SMatrix = ComputeSLayer(f, path_TE, p_v, d_v, er_v, n_layer);
%S11_XX = squeeze(SMatrix(1,1,:));
%S21_XX = squeeze(SMatrix(2,1,:));

figure(3)
subplot(2,2,1)
plot(f/1e9,(imag(S11_YY)),'r','Linewidth',1.2)
xlim([16 32])
grid ON
xlabel('Frequency (GHz)')
ylabel('S11 Im')
set(gca,'FontSize',15)
set(gca,'FontName','Times New Roman')

subplot(2,2,2)
plot(f/1e9,(real(S11_YY)),'r','Linewidth',1.2)
xlim([16 32])
grid ON
xlabel('Frequency (GHz)')
ylabel('S11 Re')
set(gca,'FontSize',15)
set(gca,'FontName','Times New Roman')

subplot(2,2,3)
plot(f/1e9,(imag(S21_YY)),'r','Linewidth',1.2)
xlim([16 32])
grid ON
xlabel('Frequency (GHz)')
ylabel('S21 Im')
set(gca,'FontSize',15)
set(gca,'FontName','Times New Roman')

subplot(2,2,4)
plot(f/1e9,(real(S21_YY)),'r','Linewidth',1.2)
xlim([16 32])
grid ON
xlabel('Frequency (GHz)')
ylabel('S21 Re')
set(gca,'FontSize',15)
set(gca,'FontName','Times New Roman')