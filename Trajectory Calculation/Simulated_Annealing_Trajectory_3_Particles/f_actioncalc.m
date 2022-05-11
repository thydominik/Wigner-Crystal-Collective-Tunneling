function action = f_actioncalc(pos, r, a, rs, N, z, dz, shift)
    action = 0;
    
    %first: calculating the interaction term (without the prefactor)
    Q = zeros(N, 1);
    for i = 1:N
        Q(i) = (rs) * (1/abs(pos(2,i) - pos(1,i)) + 1/abs(pos(3,i) - pos(1,i)) + 1/abs(pos(3,i) - pos(2,i))); 
    end
    
    %second: Potential term (without the prefactor)
    
    V = zeros(3, N);
    for i = 1:N        
        V(1, i) = 0.25 * (pos(1,i)^2 + a)^2;
        V(2, i) = 0.25 * (pos(2,i)^2 + a)^2;
        V(3, i) = 0.25 * (pos(3,i)^2 + a)^2;
    end
    
    %integrating the potential and the interaction terms:
    for i = 2:N-1
        if Q(i) + V(1, i) + V(2, i) + V(3, i) - shift == 0
            action = action;
        elseif i > 1 || i < N
            pre = ((1 - z(i)^2)/r);
            action = action + dz * 0.5 * 1/pre * (V(1, i-1) + V(1, i) + V(2, i-1) + V(2, i) + V(3, i-1) + V(3, i) + Q(i - 1) + Q(i) - 2 * shift);
        end
    end
    
    %the derivative part:
    der     = zeros(3, N-1);
    z_temp  = z + dz/2;
    for i = 2:N
        der(1, i) = (pos(1, i) - pos(1, i - 1))/dz;
        der(2, i) = (pos(2, i) - pos(2, i - 1))/dz;
        der(3, i) = (pos(3, i) - pos(3, i - 1))/dz;
    end
    
    %integrating the kinetic part:
    for i = 1:N-1
        action = action + dz * 0.5 * (1 - z_temp(i)^2)/r * (der(1, i)^2 + der(2, i)^2 + der(3, i)^2);
    end
end

