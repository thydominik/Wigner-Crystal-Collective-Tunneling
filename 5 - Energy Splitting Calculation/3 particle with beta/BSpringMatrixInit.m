function [BB] = BSpringMatrixInit(Trajectory, Alpha, tauspace, Eta, NoP, beta)
%UNTITLED:

 %first just the potential without the interaction:
    %creating the second derivative diagonal potential elements for every point
    %in the trajectory
    for time_i = 1:length(Trajectory(1, :))
        V           = [];
        U_diag      = zeros(NoP, 1);
        Uoffdiag    = zeros(NoP, NoP);
    
        x0 = Trajectory(:, time_i);
    
        for part_i = 1:NoP
            V(part_i) = (3 * x0(part_i)^2 - Alpha);
    
            for part_j = 1:NoP
                if part_i == part_j
                else
                    U_diag(part_i)              = U_diag(part_i) + Eta * (2/abs((x0(part_i) - x0(part_j) + beta)^3));
                    U_offdiag(part_i, part_j)   = - Eta * (2/abs((x0(part_i) - x0(part_j) + beta)^3));
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
    
        %this is 3 by 3 by N_division size tensor
        B(:, :, time_i) =  mtx;
    end

    BB = zeros(NoP-1, NoP-1, length(B(1,1,:)));
    
    for i = 1:(NoP-1)
        for j = 1:(NoP-1)
            for k = 1:length(tauspace(1,1,:))
                temp        =  B(:, :, k) * tauspace(:,j+1,k);
                temp2       =  tauspace(:,i+1,k)' * temp;
                BB(i,j,k)   =  temp2;
            end
        end
    end


end

