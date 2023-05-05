function action = f_actioncalc(pos, r, a, rs, N, Nt, z, dz, shift)
    action = 0;
    
    %first calculating the interaction term (without the prefactor)
    Q = zeros(Nt, 1);
    for i = 1:Nt
        for j = 1:N
            for k = (j+1):N
                %Q(i) = (rs) * (1/abs(pos(2,i) - pos(1,i)) + 1/abs(pos(3,i) - pos(1,i)) + 1/abs(pos(3,i) - pos(2,i))); 
                Q(i) = Q(i) + rs/abs(pos(j, i) - pos(k, i));
            end
        end
    end
    %second: Potential term (without the prefactor)
    
    V = zeros(N, Nt);
    for i = 1:Nt
        for j = 1:N
            V(j, i) = 0.25 * (pos(j, i)^2 + a)^2;
        end
    end
    
    %integrating the potential and the interaction terms:
    for i = 2:Nt-1
        %if Q(i) + V(1, i) + V(2, i) + V(3, i) - shift == 0% + V(4, i) + V(5, i) - shift == 0
        %    action = action;
        %elseif i > 1 || i < Nt
            pre     = ((1 - z(i)^2)/r);
            func1   = sum(V(:, i)) + Q(i);
            func2   = sum(V(:, i - 1)) + Q(i - 1);
            action  = action + dz * 0.5 * 1/pre * (func1 + func2 - (2 * shift));
        %end
    end
    
    %the derivative part:
    der     = zeros(N, Nt-1);
    z_temp  = z + dz/2;
    for i = 2:Nt
        for j = 1:N
            der(j, i) = (pos(j, i) - pos(j, i - 1))/dz;
        end
    end
    
    %integrating the kinetic part:
    for i = 1:Nt-1
        func = sum(der(:, i).^2);
        action = action + dz * 0.5 * (1 - z_temp(i)^2)/r * (func);
    end
end

