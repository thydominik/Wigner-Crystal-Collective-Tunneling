function [I_S1] = f_int_S1(z, q_0, C, div, leg, dleg, N_Lp, alpha, r)
    dz = z(2) - z(1);
    %we need q_0's deriative which is now a tanh(..) + legendreP 
    %first the derivative of c1*tanh(atanh(z)*sqrt(alpha/2)*r)
    for i = 1:div
        prefactor = sech(sqrt(alpha/2) * r *atanh(z(i)))^2;
        if prefactor == 0
            dq_0(i) = 0;
        else
            dq_0(i) = prefactor * (r/(1-z(i)^2))*alpha /sqrt(2);
        end
    end
    %then we need the derivatives of the Legendre polynomials:
    for i = 1:N_Lp
        dq_0 = dq_0 + dleg(i,:) .* C(i);
    end
    
    prefactor   = (1 - z.^2)/r;
    for i = 1:N_Lp       
        S1_func = ((1 - z.^2)./(r) .* dq_0 .* dleg(i, :)) + (r./(1 - z.^2) .* (q_0.^3 - alpha .* q_0) .* leg(i,:));
        %trapezoid integration:
        I_S1(i) = 0;
        for j = 2:div
            I_S1(i) = I_S1(i) + dz * (S1_func(j) + S1_func(j - 1))/2;
        end
    end
    
    
end

