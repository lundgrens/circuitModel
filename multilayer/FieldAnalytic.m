function [Field_x, Field_y, steps_x, dx] = FieldAnalytic(path)

    px = 9.6e-3;
    py = 9.6e-3;
    wLowTM = 1e-3;
    wLowTE = 1e-3;
    LLowTM = 7.8e-3 + 0.00001i;
    LLowTE = 8e-3;
    
    steps_x = 256;
    lim = px/2;
    dx = 2*lim/steps_x;
    x_v = -lim:dx:lim-dx;
    
    steps_y = 256;
    lim = py/2;
    dy = 2*lim/steps_y;
    y_v = -lim:dy:lim-dy;
    
    [X,Y]=meshgrid(x_v,y_v);
    
    if  strcmp('TM', path)
        Field_x = cos(pi*Y/LLowTM).*(1 - (2*Y/LLowTM).^2).^(-1/2).*rectpuls(Y/LLowTM).*rectpuls(X/wLowTM);
        Field_y = zeros(size(Field_x));
    elseif strcmp('TE', path)
        Field_y = cos(pi*X/LLowTE).*(1 - (2*X/LLowTE).^2).^(-1/2).*rectpuls(X/LLowTE).*rectpuls(Y/wLowTE);
        Field_x = zeros(size(Field_y));
    end

    figure(1)
    imagesc(x_v/1e-3,y_v/1e-3,abs(Field_x) + abs(Field_y))
    colormap jet;
    colorbar
    title('absolute transverse field')
    xlabel('x-direction [mm]')
    ylabel('y-direction [mm]')
    set(gca,'FontSize',15)
    end