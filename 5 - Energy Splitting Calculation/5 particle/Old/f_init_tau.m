function [ModVec1, EigVal] = f_init_tau(trajectory, a, v, eta, NoP)
    V           = [];
    U_diag      = zeros(NoP, 1);
    Uoffdiag    = zeros(NoP, NoP);

    x0 = trajectory(:, 1);

    for part_i = 1:NoP
        V(part_i) = (3 * x0(part_i)^2 - a);

        for part_j = 1:NoP
            if part_i == part_j
            else
                U_diag(part_i)              = U_diag(part_i) + eta * (2/abs((x0(part_i) - x0(part_j))^3));
                U_offdiag(part_i, part_j)   = - eta * (2/abs((x0(part_i) - x0(part_j))^3));
            end
        end
    end

    for i = 1:NoP
        for j = 1:NoP
            if i == j
                mtx(i, j) = V(i) + U_diag(i);
            else
                mtx(i, j) = U_offdiag(i, j);
            end
        end
    end

%eigenvectors and values of 'mtx'                     
[ModVec1, EigVal] = eig(mtx);

disp('Eigenvectors (columns)')
disp(num2str(ModVec1))
disp('Frequency Eigevalues (dimles)')
disp(num2str(((EigVal))))
end

