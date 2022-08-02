function [tauspace] = f_tau_space(R, tau, NoP)

    tauspace = zeros(NoP, NoP, length(R(3,3,:)));
    for i = 1:NoP
        temp_tau(i, :)      = tau(:, i);
        tauspace(:, i, 1)   = temp_tau(i, :);
    end
    

    for i = 2:length(R(3,3,:))
        for part_i = 1:NoP
            temp_tau(i, :) = R(:, :, i) * temp_tau(i - 1, :)';
            tauspace(:, part_i, i) = temp_tau(part_i, :);
        end
    end
    
end

