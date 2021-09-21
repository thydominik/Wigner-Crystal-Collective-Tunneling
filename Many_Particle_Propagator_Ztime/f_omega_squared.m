function [omega] = f_omega_squared(T, B, v, C)

temp_mtx = zeros(2, 2, length(C(1,:)));
for k = 1:length(C(1,:))
    %sum for the time
    for i = 2:3                                             %length(C(:,1))
        %sum for the \alpha indices
        for j = 2:3                                         %length(C(:,1))
            %sum for the \beta indices
            
            product         = C(i,k)*C(j,k);
            temp_mtx(i - 1,j - 1,k) = product;           
        end
    end
    temp_mtx(:,:,k) = temp_mtx(:,:,k) * 3;
    omega(:,:,k)    =  B(:,:,k) - temp_mtx(:,:,k);


end

end


