function [S_21, S_23] = const_integrals(Legendre, Legendre_dif, N, div, z, r, alpha)
    dz = z(2) - z(1);
    S_21 = zeros(N, N);
    for n = 1:N
        Legendre_temp_n = (1-z.^2).*Legendre_dif(n,:) - 2*z.*Legendre(n,:);
        for m = 1:N
            Legendre_temp_m = (1-z.^2).*Legendre_dif(m,:) - 2*z.*Legendre(m,:);
            for i = 1:(div-1)
                func_val = (((1 - z(i)^2)/r) * Legendre_temp_n(i) * Legendre_temp_m(i)) + (((1 - z(i + 1)^2)/r) * Legendre_temp_n(i + 1) * Legendre_temp_m(i + 1));
                S_21(n, m) = S_21(n, m) + dz/2 * func_val;
            end
        end
    end
    
    S_23 = zeros(N, N);
    for n = 1:N
        for m = 1:N
            for i = 1:(div-1)
                func_val = ((r/2 * alpha * (1 - z(i)^2)) * Legendre(n,i) * Legendre(m,i)) + ((r/2 * alpha * (1 - z(i+1)^2)) * Legendre(n,i+1) * Legendre(m,i+1));
                S_23(n, m) = S_23(n, m) - dz/2 * func_val;
            end
        end
    end
    
end

