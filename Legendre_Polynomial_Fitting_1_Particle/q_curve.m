function q_0 = q_curve(C, Legendre, N, div, alpha)
    q_0 = zeros(1, div);
    z = linspace(-1,1,div);
    for i = 1:N
        q_0 = q_0 + C(i) * Legendre(i,:) .* (1 - z.^2);
    end
    q_0 = q_0 + sqrt(alpha) * Legendre(1,:);
end

