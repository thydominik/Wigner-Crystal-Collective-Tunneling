function [f, delta_e, Normalization, Nominator] = f_f_vector(e)
    
    delta_e = zeros(3, length(e));
    for i= 1:(length(e) - 1)
        delta_e(:, i) = e(:, i +1) - e(:, i);
    end
    
    f = zeros(3, length(e));
    Normalization = f;
    Nominator = f;
    
    for i  = 1:length(e)
        n1 = norm(delta_e(:, i));
        n2 = e(:, i)' * delta_e(:, i);
        
        Normalization (i) = sqrt(abs(n1^2 + n2^2));
        Nominator(:, i) = delta_e(:, i) - ((e(:, i) * delta_e(:, i)') * e(:,i));
        f(:, i) = Nominator(:, i)/ Normalization(i);
    end
end

