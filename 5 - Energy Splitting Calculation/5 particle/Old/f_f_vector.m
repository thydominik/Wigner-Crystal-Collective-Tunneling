function [f, delta_e, Normalization, Nominator] = f_f_vector(e, NoP)
    
%determining the difference between neighboring unit vectors
    delta_e = zeros(NoP, length(e));
    for i= 1:(length(e) - 1)
        delta_e(:, i) = e(:, i +1) - e(:, i);
    end
    
    f = zeros(NoP, length(e));
    Normalization = f;
    Nominator = f;
    
%creating a normalized and perpendicular f vector to e
    for i  = 1:length(e)
        n1 = norm(delta_e(:, i));
        n2 = e(:, i)' * delta_e(:, i);
        
        Normalization (i) = sqrt(abs(n1^2 + n2^2));
        Nominator(:, i) = delta_e(:, i) - ((e(:, i)' * delta_e(:, i)) * e(:,i));
        if Normalization(i) == 0
            f(:,i) = f(:,i-1);
        else
            f(:, i) = Nominator(:, i)/ Normalization(i);
        end
    end
end

