function E_0 = f_initial_energy(N_p, alpha, eq_pos, eta)
%F_INITIAL_ENERGY 
    E_0 = zeros(N_p,1);

    for i = 1:N_p
        E_0(i) =  (0.5 * (eq_pos(i)^2 + alpha)^2);
        for j = 1:N_p
           if j ~= i
               E_0(i) = E_0(i) + (eta * 2 /(abs(eq_pos(i) - eq_pos(j))));
           end
        end
    end
    E_0 = E_0;
end

