function [omega] = f_omega_squared(B, v, C, NoP)

temp_mtx = zeros(NoP - 1, NoP - 1, length(C(1,:)));
for k = 1:length(C(1,:))
    %sum for the time
    for i = 2:NoP                                             %length(C(:,1))
        %sum for the \alpha indices
        for j = 2:NoP                                        %length(C(:,1))
            %sum for the \beta indices
            
            product         = C(i,k) * C(j,k);
            temp_mtx(i - 1,j - 1,k) = product;           
        end
    end
    temp_mtx(:,:,k) = temp_mtx(:,:,k) * 3;
    omega(:,:,k)    =  B(:,:,k) - temp_mtx(:,:,k);


end

% omega(:, :, end - 1) = diag(flip(diag(omega(:,:,2))));
% omega(:,:,end) = diag(flip(diag(omega(:,:,1))));

% omega(1,2,end/2 + 1 : end) = -omega(1,2,1:end/2);
% omega(2,1, : ) = omega(1,2,:);
% 
% omega(:,:,end) = omega(:,:,1);
% omega(1,1,end) = omega(2,2,1);
% omega(2,2,end) = omega(1,1,1);
end


