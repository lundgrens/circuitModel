function ABCD = PiABCD(Y1, Y2, Y3, f)
    
    ABCD = zeros(2,2,length(f));
    ABCD(1,1,:) = 1 + Y2./Y3;
    ABCD(1,2,:) = 1./Y3;
    ABCD(2,1,:) = Y1 + Y2 + Y1.*Y2./Y3;
    ABCD(2,2,:) = 1 + Y1./Y3;
end