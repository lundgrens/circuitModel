
% Implementation 2D-analytical model with field profile imported from CST

c = 299792458;
mu0 = 4*pi*10^(-7);
eps0 = 1/(c^2*mu0);
Z0 = sqrt(mu0/eps0);
Y0 = 1/Z0;

% Simulated field-profiles
%path = 'field_profiles/e-field (f=22.5)_TE.txt'; %TE-polarization rectangular slot
path = 'field_profiles/e-field (f=22.5)_TM.txt'; %TM-polarization rectangular slot
%path = 'field_profiles/e-field (f=9) [1].txt'; %TE-polarization CSRR

px = 10e-3; % Periodicity of the structure in x-dir. [m]
py = 10e-3; % Periodicity of the structure in y-dir. [m]

f = linspace(16,32e9,1002); % Frequency band [Hz]
omega = 2*pi*f;
k0 = omega/c;

N = 50;
n = linspace(1,N,N);
n = [flip(-n) 0 n];

M = 50;
m = linspace(1,M,M);
m = [flip(-m) 0 m];

[N_m,M_m]=meshgrid(n,m);

kn = 2*pi/px*N_m;
km = 2*pi/py*M_m;

[Field_x, Field_y, steps_x, dx] = Read_Field(path);
steps_y = steps_x;
dy = dx;

y_pad = 2^(nextpow2(steps_y)+4);
x_pad = 2^(nextpow2(steps_x)+4);

Fx = fftshift(fft2(Field_x,y_pad,x_pad)); % Compute padded FFT
Fy = fftshift(fft2(Field_y,y_pad,x_pad)); % Compute padded FFT

Fs_x = 2*pi/dx;
dF_x = Fs_x/x_pad;
kv_x = -Fs_x/2:dF_x:Fs_x/2-dF_x;

Fs_y = 2*pi/dy;
dF_y = Fs_y/y_pad;
kv_y = -Fs_y/2:dF_y:Fs_y/2-dF_y;

Eh_x = zeros(size(kn));
Eh_y = zeros(size(kn));

for j = 1:length(kn)
    dif = abs(kv_x-kn(1,j));
    [~,index_x] = min(dif);
    for k = 1:length(km)
        dif = abs(kv_y-km(k,1));
        [~,index_y] = min(dif);  
        Eh_x(k,j) = (Fx(index_y,index_x));
        Eh_y(k,j) = (Fy(index_y,index_x));
    end
end

% Nh = E_a(kth) dot e_h
Nh_TM = Eh_x.*kn./sqrt(kn.^2+km.^2) + 1*Eh_y.*km./sqrt(kn.^2+km.^2);
Nh_TE = Eh_x.*km./sqrt(kn.^2+km.^2) - 1*Eh_y.*kn./sqrt(kn.^2+km.^2);

N0 = Eh_y(N+1,N+1); % Eh_y for inc. y-pol.
Nh_TE(N+1,N+1) = N0;
Nh_TM(N+1,N+1) = N0;

Yeq = zeros(1,length(f));
for i = 1:length(f)
    betah0 = -1i*sqrt((abs(kn).^2 + abs(km).^2) - k0(i)^2);
    betah1 = -1i*sqrt((abs(kn).^2 + abs(km).^2) - k0(i)^2);

    Yh0_TM = omega(i)*eps0./betah0;
    Yh1_TM = omega(i)*eps0./betah1;
    
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
S11 = (abs(N0)^2*Y0 - abs(N0)^2*Y0 - Yeq)./(abs(N0)^2*Y0 + abs(N0)^2*Y0 + Yeq);
S21 = 2*abs(N0)^2*Y0./(abs(N0)^2*Y0 + abs(N0)^2*Y0 + Yeq);

figure(3)
plot(f/1e9,abs(S11),'r','Linewidth',1.2)
xlim([16 32])
grid ON
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
xlabel('Frequency (GHz)')
ylabel('Reflection coefficient')

figure(4)
plot(f/1e9,abs(S21),'r','Linewidth',1.2)
xlim([16 32])
grid ON
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
xlabel('Frequency (GHz)')
ylabel('Transmission coefficient')
