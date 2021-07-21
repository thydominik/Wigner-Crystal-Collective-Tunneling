function [S_int] = action(q0, Legendre, Legendre_dif, C, r, alpha, div, z, N)
    for j = 1:N
        q_d(j,:)    = C(j) * (Legendre_dif(j,:) .* (1 - z.^2) - (2 .* z .* Legendre(j,:)));
    end
    q_dif = zeros(1,div);
    
    for i = 1:N
        q_dif   = q_dif + q_d(i,:);
    end 

    S1 = 0;
    for i = 1:div
        S1 = S1 + ((1 - z(i)^2)/r) * q_d(i)^2;
    end
    
    S2 = 0;
    for i = 1:div
        S2 = S2 + r/(1 - z(i)^2) * ((1/4) * q0(i)^4 - (alpha/2 * q0(i)^2));
    end
    S_int = S1 + S2;
end

