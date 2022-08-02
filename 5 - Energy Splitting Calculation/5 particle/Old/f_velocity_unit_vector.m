function e = f_velocity_unit_vector(v, NoP)
    e = v;
    
    for i = 1:length(v)
        normalization = norm(v(:, i));
        
        %the first vector on paper should be a vector pointing at the
        %instantons direction
        if v(:, i) == 0
            warning('zero velocity :(')
            e(:, i) = zeros(1, NoP);
        else
            e(:, i) = v(:,i) / normalization;
        end
    end
end

