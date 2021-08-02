function [I] = f_integration(chi, Lp, N, z)
    %trapezoid integration rule
    I = 0;
    dz = (z(2) - z(1))/2;
    for i = 2:N
        I = I + dz * (((chi(i)*Lp(i)) + (chi(i-1)*Lp(i-1))));
    end
end

