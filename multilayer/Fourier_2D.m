function [Fx, Fy, kv_x, kv_y] = Fourier_2D(Field_x, Field_y, dx, N)

steps_x = N;
steps_y = N;
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
     
end