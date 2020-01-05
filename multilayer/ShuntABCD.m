function ABCD = ShuntABCD(Y,f)
    ABCD = zeros(2,2,length(f));
    ABCD(1,1,:) = 1;
    ABCD(2,1,:) = Y;
    ABCD(2,2,:) = 1;
end