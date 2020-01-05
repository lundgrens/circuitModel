function [Ype, Ypi, Ys] = ComputeAdmittance(f, d, er, kn, km, Nh_TM, Nh_TE, N0, N)
    
    c = 299792458;
    mu0 = 4*pi*10^(-7);
    eps0 = 1/(c^2*mu0);
    Z0 = sqrt(mu0/eps0);
    Y0 = 1/Z0;
    
    omega = 2*pi*f;
    k0 = omega/c;

    Ype = zeros(1,length(f));
    Ypi = zeros(1,length(f));
    Ys = zeros(1,length(f));
    for i = 1:length(f)
        betah0 = -1i*sqrt((abs(kn).^2 + abs(km).^2) - er*k0(i)^2);
        betah1 = -1i*sqrt((abs(kn).^2 + abs(km).^2) - er*k0(i)^2);

        Yh0_TM = omega(i)*eps0*sqrt(er)./betah0;
        Yh1_TM = omega(i)*eps0*sqrt(er)./betah1;

        Yh0_TE = betah0./(omega(i)*mu0);
        Yh1_TE = betah1./(omega(i)*mu0);

        YhL_TM = Yh0_TM;
        YhR_TM = 1i*Yh1_TM.*tan(betah1*d/2);

        YhL_TE = Yh0_TE;
        YhR_TE = 1i*Yh1_TE.*tan(betah1*d/2);

        Ype(i) = sum(sum(abs(Nh_TM).^2.*YhL_TM))...
            + sum(sum(abs(Nh_TE).^2.*YhL_TE)) - 2*abs(N0)^2*Y0;

        Ypi(i) = sum(sum(abs(Nh_TM).^2.*YhR_TM))...
            + sum(sum(abs(Nh_TE).^2.*YhR_TE)) - abs(N0)^2*YhR_TM(N+1,N+1);

        Ys(i) = -1i*sum(sum(abs(Nh_TM).^2.*Yh1_TM.*csc(betah1*d)))...
            -1i*sum(sum(abs(Nh_TE).^2.*Yh1_TE.*csc(betah1*d)))...
            + 1i*abs(N0)^2*Y0*csc(betah1(N+1,N+1)*d);

    end
end