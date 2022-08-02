function [omega] = OmegaSquared(B, C, NoP, Nt)
%OMEGASQUARED:

temp_mtx = zeros(NoP - 1, NoP - 1, Nt);

for k = 1:Nt
    for alpha_ind = 2:NoP
        for beta_ind = 2:NoP
            product = C(alpha_ind, k) * C(beta_ind, k);
            temp_mtx(alpha_ind - 1, beta_ind - 1, k) = product;
        end
    end
    temp_mtx(:, :, k) = temp_mtx(:, :, k) * 3;
    omega(:, :, k) = B(:, :, k) - temp_mtx(:, :, k);
end

end

