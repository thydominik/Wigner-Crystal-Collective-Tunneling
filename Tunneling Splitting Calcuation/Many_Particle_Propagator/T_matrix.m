function [T] = T_matrix(tau)
    T = zeros(length(tau(1,:,1)), 2, length(tau(1,1,:)));
    for i = 1: length(tau(1,1,:))
        T(:,:,i) = tau(:,2:3,1);
    end
end

