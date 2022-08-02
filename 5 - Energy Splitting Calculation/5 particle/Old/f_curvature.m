function [C] = f_curvature(e, tauspace, z_time, v, NoP)
    C = zeros(NoP,length(tauspace(1,1,:)));

    %numerically derivating the tauspace vectors
    for i = 2: length(C)-1
        %dt here is 2 times dz
        dt              = z_time(3) - z_time(1);
        for part_i = 1:NoP
            derivative(part_i, :) = (tauspace(:, part_i, i+1) - tauspace(:, part_i, i-1))/(dt);
            C(part_i, i) = derivative(part_i, :) * tauspace(:, part_i, i);
        end       
    end
end
