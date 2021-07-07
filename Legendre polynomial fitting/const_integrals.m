function [S_21, S_23] = const_integrals(Legendre, Legendre_dif, N, div, z, r, alpha)
    
    S_21 = zeros(N, N);
    for n = 1:N
        Legendre_temp_n = (1-z.^2).*Legendre_dif(n,:) - 2*z.*Legendre(n,:);
        for m = 1:N
            Legendre_temp_m = (1-z.^2).*Legendre_dif(m,:) - 2*z.*Legendre(m,:);
            for i = 1:div
                S_21(n, m) = S_21(n, m) + (((1 - z(i)^2)/r) * Legendre_temp_n(i) * Legendre_temp_m(i));
            end
        end
    end
    
    S_23 = zeros(N, N);
    for n = 1:N
        for m = 1:N
            for i = 1:div
                S_23(n, m) = S_23(n, m) - ((r/2 * alpha * (1 - z(i)^2)) * Legendre(n,i) * Legendre(m,i));
            end
        end
    end
    
end

