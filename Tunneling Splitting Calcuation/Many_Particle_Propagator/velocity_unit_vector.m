function [unit_vector] = velocity_unit_vector(velocity)   
    unit_vector = velocity;
    
    for i = 1:length(velocity)
        normalization = norm(velocity(:,i));        
        unit_vector(:,i) = velocity(:,i)/normalization;
    end
end

