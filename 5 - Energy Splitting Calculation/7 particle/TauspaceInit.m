function [Tau0, EigVals, Tauspace] = TauspaceInit(Trajectory, Alpha, Velocity, Eta, NoP, Nt, Rmtx)
    %TAUSPACEINIT:
    
    V           = [];
    U_diag      = zeros(NoP, 1);
    Uoffdiag    = zeros(NoP, NoP);
    
    x0 = Trajectory(:, 1);
    
    for part_i = 1:NoP
        V(part_i) = (3 * x0(part_i)^2 - Alpha);
    
        for part_j = 1:NoP
            if part_i == part_j
            else
                U_diag(part_i)              = U_diag(part_i) + Eta * (2/abs((x0(part_i) - x0(part_j))^3));
                U_offdiag(part_i, part_j)   = - Eta * (2/abs((x0(part_i) - x0(part_j))^3));
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
    [Tau0, EigVals] = eig(mtx);
    disp('Eigenvectors (columns)')
    disp(num2str(Tau0))
    disp('Frequency Eigevalues (dimles)')
    disp(num2str(sqrtm(EigVals)))


    Tauspace = zeros(NoP, NoP, Nt);
    
    for part_i = 1:NoP
        Tauspace(:, part_i, 1) = Tau0(:, part_i);
    end
    
    for i = 2:Nt
        for part_i = 1:NoP
            Tauspace(:, part_i, i) = Rmtx(:, :, i) * Tauspace(:, part_i, i - 1);
        end
    end

end

