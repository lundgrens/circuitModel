function [Nh_TM, Nh_TE, N0] = ComputeNh(Fx, Fy, kn, km, kv_x, kv_y, N)

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
    
    if abs(Eh_y(N+1,N+1))>abs(Eh_x(N+1,N+1))
        N0 = Eh_y(N+1,N+1);
    else
        N0 = Eh_x(N+1,N+1);
    end
    Nh_TE(N+1,N+1) = N0;
    Nh_TM(N+1,N+1) = N0;
     
end