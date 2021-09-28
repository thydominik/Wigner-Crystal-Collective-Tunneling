function e = f_velocity_unit_vector(v)
    e = v;
    
    for i = 1:length(v)
        normalization = norm(v(:, i));
        
        %the first vector on paper should be a vector pointing at the
        %instantons direction
        if v(:, i) == 0
            warning('zero velocity :(')
            e(:, i) = [0, 0, 0];
        else
            e(:, i) = v(:,i) / normalization;
        end
    end
end

