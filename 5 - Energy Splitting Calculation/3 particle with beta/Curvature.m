function [C] = Curvature(Tauspace, z, dz, v, NoP, Nt)
    %CURVATURE:
    
    C = zeros(NoP, Nt);

    for i = 2:Nt-1
        dt = dz * 2;
        for part_i = 1:NoP
            derivative(part_i, :) = (Tauspace(:, part_i, i + 1) - Tauspace(:, part_i, i - 1))/(dt);
            C(part_i, i) = derivative(part_i, :) * Tauspace(:, part_i, i); 
        end

    end
end

