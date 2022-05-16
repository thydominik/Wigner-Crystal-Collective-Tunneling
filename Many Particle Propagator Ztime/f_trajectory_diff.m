function diffcurve = f_trajectory_diff(N, traj, z )
%numerical derivative of the curve
    diffcurve = zeros(3, N);
    
    %in Z time, the dz is const.
    dz = z(2) - z(1);
    for i = 1:3
        for j = 1:N
            
            %exactly at the endpoints the particle should have no mom.
            if j == 1
                diffcurve(i, j) = (traj(i, j + 1) - traj(i, j))/ dz;
            elseif j == N
                diffcurve(i, j) = (traj(i, j) - traj(i, j - 1))/ dz;
            else
                diffcurve(i, j) = (traj(i, j + 1) - traj(i, j - 1))/(2 * dz);
            end
        end
    end
    %it is a question that what should we do with the large step at the
    %first derivative
end