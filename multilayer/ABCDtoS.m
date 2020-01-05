function SMatrix = ABCDtoS(ABCD,Y_ref)
    SMatrix = zeros(size(ABCD));

    A = ABCD(1,1,:);
    B = ABCD(1,2,:);
    C = ABCD(2,1,:);
    D = ABCD(2,2,:);
    
    SMatrix(1,1,:) = (A + B*Y_ref - C/Y_ref - D)./(A + B*Y_ref + C/Y_ref + D);
    SMatrix(2,1,:) = 2./(A + B*Y_ref + C/Y_ref + D);
end