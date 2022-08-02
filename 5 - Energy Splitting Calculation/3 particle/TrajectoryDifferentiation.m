function [DiffCurve, VeloUnitVec] = TrajectoryDifferentiation(Trajecotories, NoP, Nt, time, dtime, SmoothSwitch)
    %TRAJECTORYDIFFERENTIATION: Takes a set of trajectories and calculates the
    %" velocity" based on the "time" data provided
    
    DiffCurve = zeros(NoP, Nt);
    
    for part_ind = 1:NoP
        for traj_ind = 1:Nt
            if traj_ind == 1
                DiffCurve(part_ind, traj_ind) = abs(Trajecotories(part_ind, traj_ind + 1) - Trajecotories(part_ind, traj_ind))/dtime;
            elseif traj_ind == Nt
                DiffCurve(part_ind, traj_ind) = abs(Trajecotories(part_ind, traj_ind) - Trajecotories(part_ind, traj_ind - 1))/dtime;
            else
                DiffCurve(part_ind, traj_ind) = abs(Trajecotories(part_ind, traj_ind + 1) - Trajecotories(part_ind, traj_ind - 1))/(2 * dtime);
            end
        end
    end

%     DiffCurve(:, 1)     = 0;
%     DiffCurve(:, Nt)    = 0;

    if SmoothSwitch 
        for i = 1:NoP
            DiffCurve(i, :) = smooth(DiffCurve(i, :));
        end
    end

    %Calculating the velocity Unit vector as well

    for traj_ind = 1:Nt
        Normalization = norm(DiffCurve(:, traj_ind));

        if Normalization == 0
            VeloUnitVec(:, traj_ind) = zeros(1, NoP);
        else
            VeloUnitVec(:, traj_ind) = DiffCurve(:, traj_ind) / Normalization;
        end
    end



end



