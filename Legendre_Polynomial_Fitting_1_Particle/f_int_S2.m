function [I_S2] = f_int_S2(z, q_0, div, leg, dleg, N_Lp, alpha, r)
    dz = z(2) - z(1);
    
    I_S2 = zeros(N_Lp, N_Lp);
    
    for i = 1:N_Lp
        for j = 1:N_Lp
            L_tilde_i   = (1 - z.^2) .* dleg(i,:) - 2 .* z .* leg(i,:);
            L_tilde_j   = (1 - z.^2) .* dleg(j,:) - 2 .* z .* leg(j,:);
            S2_func_p1  = (1 - z.^2)./(2 * r) .* L_tilde_i .* L_tilde_j;
            S2_func_p2  = (r/2) .* (1 - z.^2) .* leg(i,:) .* leg(j,:) .* (3 .* q_0.^2 - alpha);
            S2_func     = S2_func_p1 + S2_func_p2;
            for k = 2:div
                I_S2(i,j)   = I_S2(i,j) + dz * (S2_func(k) + S2_func(k - 1))/2;                
            end
            
        end
    end    
end

