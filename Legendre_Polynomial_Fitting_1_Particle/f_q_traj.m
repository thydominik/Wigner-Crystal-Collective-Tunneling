function [q_0] = f_q_traj(z, C, Legendre, N, r, alpha)

    q_0 = sqrt(alpha) * tanh(atanh(z) * sqrt(alpha/2) * r );
    for i = 1:N
        q_0 =  q_0 + (C(i) * Legendre(i,:)) .* (1 - z.^2);
    end
end

