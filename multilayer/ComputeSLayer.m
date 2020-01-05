% Johan Lundgren 2019-11-04

function SMatrix = ComputeSLayer(f, path, p_v, d_v, er_v, n_layer)
    
    %%%%% General Constants %%%%%
    c = 299792458;
    mu0 = 4*pi*10^(-7);
    eps0 = 1/(c^2*mu0);
    Z0 = sqrt(mu0/eps0);
    Y0 = 1/Z0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%% Number of modes %%%%%%
    N = 50;
    n = linspace(1,N,N);
    n = [flip(-n) 0 n];

    M = 50;
    m = linspace(1,M,M);
    m = [flip(-m) 0 m];

    [N_m,M_m]=meshgrid(n,m);

    kn = 2*pi/p_v(1)*N_m;
    km = 2*pi/p_v(2)*M_m;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Use either analytical field profile or exported data
    [Field_x, Field_y, nx, dx] = Read_Field(path);
    %[Field_x, Field_y, nx, dx] = FieldAnalytic(path);
    
    [Fx, Fy, kv_x, kv_y] = Fourier_2D(Field_x, Field_y, dx, nx);
    
    [Nh_TM, Nh_TE, N0] = ComputeNh(Fx, Fy, kn, km, kv_x, kv_y, N);
    
    [Ype, Ypi, Ys] = ComputeAdmittance(f, d_v(1), er_v(1), kn, km, Nh_TM, Nh_TE, N0, N);
    
    A_shunt = ShuntABCD(Ype, f);
    A_pi = PiABCD(Ypi, Ypi, Ys, f);
    ABCD = CascadeABCD(A_shunt, A_pi);
    
    for i=3:n_layer
        [Ype, Ypi, Ys] = ComputeAdmittance(f, d_v(i-1), er_v(i-1), kn, km, Nh_TM, Nh_TE, N0, N);
        A_pi = PiABCD(Ypi, Ypi, Ys, f);
        ABCD = CascadeABCD(ABCD,A_pi);
    end
    A_shunt = ShuntABCD(Ype, f);
    ABCD = CascadeABCD(ABCD,A_shunt);
    
    Y_ref = abs(N0)^2*Y0;
    SMatrix = ABCDtoS(ABCD,Y_ref);
    


end