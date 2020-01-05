function [Field_x, Field_y, N, dx] = Read_Field(path)
% Input: path to .txt file with simulation data
% Return: Field_x = Matrix with x-component of field profile
%         Field_y = Matrix with y-component of field profile
%         N = number of points in one direction
%         dx = step size between two points
% Johan Lundgren 2019

% y-component correspond to TM mode and x-component corrspeond to TE-mode
    A = readtable(path);

    x_v = 1e-3*table2array(A(:,1));
    y_v = 1e-3*table2array(A(:,2));

    Ex_Re = table2array(A(:,4));
    Ey_Re = table2array(A(:,5));
    %Ez_Re = table2array(A(:,6));

    Ex_Im = table2array(A(:,7));
    Ey_Im = table2array(A(:,8));
    %Ez_Im = table2array(A(:,9));

    N=1;
    k=1;
    while k ~= 0
        k = abs(x_v(1))-abs(x_v(N+1));
        N = N + 1;
    end

    Field_x = zeros(N,N);
    Field_y = zeros(N,N);
    for i=1:N
        Field_x(i,1:N) = Ex_Re(1+N*(i-1):N+N*(i-1))...
            + 1i*Ex_Im(1+N*(i-1):N+N*(i-1));
        Field_y(i,1:N) = Ey_Re(1+N*(i-1):N+N*(i-1))...
            + 1i*Ey_Im(1+N*(i-1):N+N*(i-1));
    end

    dx = abs(x_v(1)-x_v(2));

    figure(1)
    subplot(2,2,1)  
    imagesc(x_v/1e-3,y_v/1e-3,real(Field_x))
    colormap jet;
    colorbar
    title('Re x-component')
    xlabel('x-direction [mm]')
    ylabel('y-direction [mm]')
    set(gca,'FontSize',15)

    subplot(2,2,2) 
    imagesc(x_v/1e-3,y_v/1e-3,real(Field_y))
    colormap jet;
    colorbar
    title('Re y-component')
    xlabel('x-direction [mm]')
    ylabel('y-direction [mm]')
    set(gca,'FontSize',15)

    subplot(2,2,3)  
    imagesc(x_v/1e-3,y_v/1e-3,imag(Field_x))
    colormap jet;
    colorbar
    title('Im x-component')
    xlabel('x-direction [mm]')
    ylabel('y-direction [mm]')
    set(gca,'FontSize',15)

    subplot(2,2,4)  
    imagesc(x_v/1e-3,y_v/1e-3,imag(Field_y))
    colormap jet;
    colorbar
    title('Im y-component')
    xlabel('x-direction [mm]')
    ylabel('y-direction [mm]')
    set(gca,'FontSize',15)

    figure(2)
    imagesc(x_v/1e-3,y_v/1e-3,abs(Field_x) + abs(Field_y))
    colormap jet;
    colorbar
    title('absolute transverse field')
    xlabel('x-direction [mm]')
    ylabel('y-direction [mm]')
    set(gca,'FontSize',15)
    end