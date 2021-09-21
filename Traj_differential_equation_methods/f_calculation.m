function ans = f_calculation(chi, N_p, alpha, eta, E0)
    for i = 1:N_p
        
        part1 = 0.5 * (chi(i)^2 + alpha)^2;
        
        part2 = 0;
        
        for j = 1:N_p
            if j ~= i
                part2 = part2 + eta /abs(chi(i) - chi(j));
            end
        end
        ans(i) = sqrt(abs(part1 + part2 - E0(i)));
    end
end

