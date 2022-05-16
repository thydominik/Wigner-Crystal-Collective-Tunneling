function [S1, S2] = integrals(Legendre, Legendre_dif, C, N, div, z, r, alpha, S_21, S_23, Q_0, dQ)
%here I will calculate the necessary integrals in order to obtain the C_n
%coefficients.

    dz      = z(2) - z(1);
    q_d     = zeros(N, div);
    q_s     = zeros(N, div);
    
    for j = 1:(N)
        %diff of q_0 with only the legendre polynomials
        q_d(j,:)    = C(j) * (Legendre_dif(j,:) .* (1 - z.^2) - (2 .* z .* Legendre(j,:)));
        %original q_0 with only the Legendre polynomials
        q_s(j,:)    = C(j) * Legendre(j,:) .* (1 - z.^2);
    end

    q_dif = zeros(1,div);
    q     = zeros(1,div);

    for i = 1:N
        q_dif   = q_dif + q_d(i,:);
        q       = q + q_s(i,:);
    end
    q = q + Q_0;        %q_0 with the guess curve
    q_dif = q_dif + dQ; %derivative of q_0 with the guess curve
%first order terms:

    S_11 = zeros(N, 1);

    for n = 1:N
        Legendre_temp = (1-z.^2).*Legendre_dif(n,:) - 2*z.*Legendre(n,:);
        for i = 1:(div-1)
            func_val = (((1 - z(i)^2)/r) * q_dif(i) * Legendre_temp(i)) + (((1 - z(i+1)^2)/r) * 2 * q_dif(i+1) * Legendre_temp(i+1));
            S_11(n) = S_11(n) + dz/2 * func_val;
        end
    end

    S_12 = zeros(N, 1);
    for n = 1:N
        for i = 1:(div-1)
            func_val = (r * q(i)^3 * Legendre(n,i)) + (r * q(i)^3 * Legendre(n,i+1));
            S_12(n) = S_12(n) + dz/2 * func_val;
        end
    end

    S_13 = zeros(N, 1);
    for n = 1:N
        for i = 1:(div-1)
            func_val = (r/2 * alpha *  q(i) * Legendre(n,i)) + (r * alpha *  q(i+1) * Legendre(n,i+1));
            S_13(n) = S_13(n) - dz/2 * func_val;
        end
    end

%second order terms:
    S_22 = zeros(N, N);
    for n = 1:N
        for m = 1:N
            for i = 1:(div-1)
                func_val = ((1.5 * r * (1 - z(i)^2)) * q(i)^2 * Legendre(n,i) * Legendre(m,i)) + (1.5 * r * (1 - z(i+1)^2)) * q(i+1)^2 * Legendre(n,i+1) * Legendre(m,i+1);
                S_22(n, m) = S_22(n, m) + dz/2 * func_val;
            end
        end
    end

    S1 = S_11 + S_12 + S_13;
    S2 = S_21 + S_22 + S_23;
end
