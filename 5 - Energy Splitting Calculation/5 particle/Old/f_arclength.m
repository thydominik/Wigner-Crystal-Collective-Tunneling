function [chiS, LS, VS] = f_arclength(traj, a, eta, zz, N)
    %F_ARCLENGTH
    
    LS(1) = 0;
    for i = 2:length(traj)
        for partIdx = 1:N
            x(partIdx) = (traj(partIdx, i) - traj(partIdx, i-1))^2;
        end
        S(i) = sqrt(sum(x));
        LS(i) = LS(i-1) + S(i);
    end
    
    for i = 1:length(traj)
        VS(i) = 0;
        for partIdx = 1:N
            VS(i) = VS(i) + 0.25 * (traj(partIdx, i )^2 - a)^2;
            for partIdx2 = (partIdx+1):N
                VS(i) = VS(i) + eta * 1/abs(traj(partIdx, i) - traj(partIdx2, i));
            end
        end
    end
    chiS = 0;
end

